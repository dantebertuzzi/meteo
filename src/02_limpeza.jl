# =============================================================================
# 02_limpeza.jl — auditoria de qualidade e agregação horária → diária → mensal → anual
#
# Nenhuma observação é removida silenciosamente. Valores só viram `missing` quando
# (i) ultrapassam um limite físico absoluto (LIMITES em 00_config.jl) ou
# (ii) o período agregado não cumpre o critério de completude. Toda exclusão é
# registrada em tabela.
# =============================================================================

nvalidos(x) = count(!ismissing, x)
agrega_se(f, x, minimo) = nvalidos(x) >= minimo ? f(skipmissing(x)) : missing

"Maior sequência de `true` consecutivos."
function max_consecutivo(v)
    best = run = 0
    for b in v
        run = b ? run + 1 : 0
        best = max(best, run)
    end
    return best
end

# ----------------------------------------------------------------------------- auditoria
"Resumo por variável horária: ausentes, faixa observada e violações de limite físico."
function auditoria_variaveis(d::DataFrame)
    linhas = map(VARS_NUM) do v
        x = d[!, v]; xs = collect(skipmissing(x)); lo, hi = LIMITES[v]
        (variavel = String(v), n = length(x), n_validos = length(xs),
         pct_ausente = round(100 * (1 - length(xs) / length(x)), digits = 2),
         minimo = minimum(xs), p01 = quantile(xs, 0.01), mediana = median(xs),
         p99 = quantile(xs, 0.99), maximo = maximum(xs),
         n_fora_limite = count(y -> y < lo || y > hi, xs))
    end
    DataFrame(linhas)
end

"Checagens de consistência interna entre variáveis na mesma hora."
function checagens_consistencia(d::DataFrame)
    regras = [
        ("Tmín(hora) > Tmáx(hora)",            [:tmin, :tmax], r -> r.tmin > r.tmax),
        ("T > Tmáx(hora) + 0,5 °C",            [:t, :tmax],    r -> r.t > r.tmax + 0.5),
        ("T < Tmín(hora) − 0,5 °C",            [:t, :tmin],    r -> r.t < r.tmin - 0.5),
        ("Td > T + 0,5 °C",                    [:td, :t],      r -> r.td > r.t + 0.5),
        ("URmín > URmáx",                      [:urmin, :urmax], r -> r.urmin > r.urmax),
        ("Pmín > Pmáx",                        [:pmin, :pmax], r -> r.pmin > r.pmax),
        ("Rajada < velocidade média",          [:vraj, :vvel], r -> r.vraj < r.vvel),
        ("Chuva > 50 mm/h (plausível; sinalizada)", [:prec], r -> r.prec > 50),
    ]
    DataFrame(map(regras) do (nome, cols, f)
        sub = dropmissing(d[:, cols])
        (checagem = nome, n_avaliado = nrow(sub), n_violacoes = count(f, eachrow(sub)))
    end)
end

"Número de sequências de ≥ `L` horas com valor idêntico (sensor travado)."
function sequencias_constantes(x, L)
    n = 0; run = 1
    for i in 2:length(x)
        if !ismissing(x[i]) && !ismissing(x[i-1]) && x[i] == x[i-1]
            run += 1
            run == L && (n += 1)
        else
            run = 1
        end
    end
    return n
end

"Ausência (%) por ano × mês para uma variável horária (base da figura de completude)."
function completude_mensal(d::DataFrame, v::Symbol)
    g = combine(groupby(transform(d, :dt_local => ByRow(year) => :ano,
                                     :dt_local => ByRow(month) => :mes), [:ano, :mes]),
                v => (x -> 100 * mean(ismissing.(x))) => :pct_ausente,
                nrow => :horas)
    # meses totalmente fora do registro contam como 100 % ausentes
    todos = DataFrame(ano = repeat(unique(g.ano), inner = 12), mes = repeat(1:12, length(unique(g.ano))))
    g = leftjoin(todos, g, on = [:ano, :mes])
    g.pct_ausente = coalesce.(g.pct_ausente, 100.0)
    sort!(g, [:ano, :mes])
end

"Diagnóstico por ano da radiação global — detecta mudança de escala/sensor."
function auditoria_radiacao(d::DataFrame)
    dh = transform(d, :dt_utc => ByRow(year) => :ano)
    combine(groupby(dh, :ano)) do s
        x = collect(skipmissing(s.rad)); x = x[x .> 0]
        isempty(x) && return (n_pos = 0, p50 = NaN, p95 = NaN, maximo = NaN, n_acima_teto = 0)
        (n_pos = length(x), p50 = quantile(x, 0.5), p95 = quantile(x, 0.95),
         maximo = maximum(x), n_acima_teto = count(>(LIMITES[:rad][2]), x))
    end
end

"""
Aplica limites físicos (valor impossível → missing) e devolve o log de decisões.
A radiação é tratada à parte: é *excluída integralmente* das análises (ver relatório).
"""
function aplica_limites!(d::DataFrame)
    log = DataFrame(variavel = String[], n_anulados = Int[], motivo = String[])
    for v in VARS_NUM
        v == :rad && continue
        lo, hi = LIMITES[v]
        idx = findall(x -> !ismissing(x) && (x < lo || x > hi), d[!, v])
        d[idx, v] .= missing
        push!(log, (String(v), length(idx), "fora do limite físico [$lo, $hi]"))
    end
    n_rad = count(!ismissing, d.rad)
    d.rad .= missing
    push!(log, ("rad", n_rad, "variável excluída: mudança de escala 2005–2008, " *
                "quebra de nível ≈2014 e 2023–24 sem dados (ver tabela de radiação)"))
    return log
end

# ----------------------------------------------------------------------------- agregação
ano_hidro(dt::Date) = month(dt) >= MES_INICIO_HIDRO ? year(dt) : year(dt) - 1
rotulo_hidro(a::Integer) = @sprintf("%d/%02d", a, (a + 1) % 100)

"Agregação diária em hora local, com critério de completude horária."
function agrega_diario(d::DataFrame)
    dh = transform(d, :dt_local => ByRow(Date) => :data)
    dd = combine(groupby(dh, :data),
        nrow => :nh,
        :prec => nvalidos => :nh_prec,
        :prec => (x -> sum(skipmissing(x); init = 0.0)) => :prec_obs,
        :t     => (x -> agrega_se(mean,    x, MIN_HORAS_DIA)) => :tmed,
        :tmax  => (x -> agrega_se(maximum, x, MIN_HORAS_DIA)) => :tmax,
        :tmin  => (x -> agrega_se(minimum, x, MIN_HORAS_DIA)) => :tmin,
        :ur    => (x -> agrega_se(mean,    x, MIN_HORAS_DIA)) => :ur,
        :urmin => (x -> agrega_se(minimum, x, MIN_HORAS_DIA)) => :urmin,
        :td    => (x -> agrega_se(mean,    x, MIN_HORAS_DIA)) => :td,
        :p     => (x -> agrega_se(mean,    x, MIN_HORAS_DIA)) => :p,
        :vvel  => (x -> agrega_se(mean,    x, MIN_HORAS_DIA)) => :vvel)
    dd.prec = [r.nh_prec >= MIN_HORAS_PREC ? r.prec_obs : missing for r in eachrow(dd)]
    dd.dtr  = dd.tmax .- dd.tmin
    dd.ano  = year.(dd.data); dd.mes = month.(dd.data); dd.ano_hidro = ano_hidro.(dd.data)
    return dd
end

const VARS_MENSAIS = [:tmed, :tmax, :tmin, :dtr, :ur, :urmin, :td, :p, :vvel]

"""
Agregação mensal. Variáveis de estado: média dos dias válidos se o mês cumprir a
regra OMM (≤ $MAX_DIAS_FALTA dias ausentes e ≤ $MAX_DIAS_CONSEC consecutivos).
Chuva: total mensal só se ≤ 2 % das horas estiverem ausentes (totais não toleram falha).
"""
function agrega_mensal(dd::DataFrame)
    linhas = NamedTuple[]
    for s in groupby(dd, [:ano, :mes])
        a, m = s.ano[1], s.mes[1]
        nd = daysinmonth(Date(a, m))
        dias = Date(a, m, 1):Day(1):Date(a, m, nd)
        pos  = Dict(s.data .=> 1:nrow(s))
        vals = Dict{Symbol, Any}()
        for v in VARS_MENSAIS
            falta = [!haskey(pos, dia) || ismissing(s[pos[dia], v]) for dia in dias]
            ok = count(falta) <= MAX_DIAS_FALTA && max_consecutivo(falta) <= MAX_DIAS_CONSEC
            vals[v] = ok ? mean(skipmissing(s[!, v])) : missing
        end
        frac = 1 - sum(s.nh_prec) / (24nd)
        push!(linhas, (; ano = a, mes = m, frac_ausente_prec = frac, prec_obs = sum(s.prec_obs),
                        prec = frac <= MAX_HORAS_FALTA_MES_PREC ? sum(s.prec_obs) : missing,
                        (v => vals[v] for v in VARS_MENSAIS)...))
    end
    dm = DataFrame(linhas)
    dm.ano_hidro = [m >= MES_INICIO_HIDRO ? a : a - 1 for (a, m) in zip(dm.ano, dm.mes)]
    dm.tempo = dm.ano .+ (dm.mes .- 0.5) ./ 12
    return dm
end

"Climatologia mensal (média dos meses válidos) e anomalias `v_anom` = v − clim(mês)."
function anomalias_mensais!(dm::DataFrame, vars)
    clims = Dict{Symbol, Dict{Int, Float64}}()
    for v in vars
        c = combine(groupby(dropmissing(dm[:, [:mes, v]]), :mes), v => mean => :c)
        clims[v] = Dict(c.mes .=> c.c)
        dm[!, Symbol(v, :_anom)] = [ismissing(x) ? missing : x - clims[v][m] for (x, m) in zip(dm[!, v], dm.mes)]
    end
    return clims
end

"""
Série anual por ano civil (temperatura, UR etc.). O valor anual é
clim_anual + média das anomalias mensais válidas, exigindo ≥ $MIN_MESES_ANO meses.
Usar anomalias evita o viés sazonal quando faltam meses (ex.: faltar o inverno
'aquece' uma média bruta).
"""
function agrega_anual(dm::DataFrame, clims, vars; chave = :ano)
    anos = sort(unique(dm[!, chave]))
    out = DataFrame(chave => anos)
    for v in vars
        clim_anual = mean(values(clims[v]))
        va = Symbol(v, :_anom)
        nm = Int[]; val = Union{Missing, Float64}[]
        for a in anos
            x = collect(skipmissing(dm[dm[!, chave] .== a, va]))
            push!(nm, length(x))
            push!(val, length(x) >= MIN_MESES_ANO ? clim_anual + mean(x) : missing)
        end
        out[!, v] = val
        out[!, Symbol(:nmeses_, v)] = nm
    end
    return out
end

"""
Índices anuais de extremos de temperatura (ano civil), padrão ETCCDI simplificado:
TXx, TNn, TX90p e TN10p (percentis mês-específicos sobre todo o período).
Só calculados nos anos em que Tmáx/Tmín têm ≥ $MIN_MESES_ANO meses válidos.
"""
function indices_temperatura(dd::DataFrame, anual::DataFrame)
    p90 = Dict(r.mes => r.q for r in eachrow(combine(groupby(dropmissing(dd[:, [:mes, :tmax]]), :mes), :tmax => (x -> quantile(x, 0.9)) => :q)))
    p10 = Dict(r.mes => r.q for r in eachrow(combine(groupby(dropmissing(dd[:, [:mes, :tmin]]), :mes), :tmin => (x -> quantile(x, 0.1)) => :q)))
    linhas = map(anual.ano) do a
        s = dd[dd.ano .== a, :]
        okx = anual.nmeses_tmax[anual.ano .== a][1] >= MIN_MESES_ANO
        okn = anual.nmeses_tmin[anual.ano .== a][1] >= MIN_MESES_ANO
        sx = dropmissing(s[:, [:mes, :tmax]]); sn = dropmissing(s[:, [:mes, :tmin]])
        (ano = a,
         txx   = okx ? maximum(sx.tmax) : missing,
         tnn   = okn ? minimum(sn.tmin) : missing,
         tx90p_n = okx ? count(r -> r.tmax > p90[r.mes], eachrow(sx)) : missing,
         tx_dias = okx ? nrow(sx) : missing,
         tn10p_n = okn ? count(r -> r.tmin < p10[r.mes], eachrow(sn)) : missing,
         tn_dias = okn ? nrow(sn) : missing)
    end
    DataFrame(linhas)
end

"""
Série por ano hidrológico (jul–jun) para precipitação e índices diários.

Critério de validade — *fração climatológica ausente* (FCA):
    FCA = Σ_dias (h_ausentes/24)·c̄(mês) / P̄_anual
onde c̄(mês) é a chuva diária climatológica do mês (dos meses completos). A FCA
estima que fração da chuva *esperada* do ano caiu em horas sem registro. Um ano
com julho ausente (≈1 % da chuva) permanece; um ano sem março–abril (≈30 %) sai.
"""
function agrega_hidro(dd::DataFrame, dm::DataFrame, clims, vars)
    cm = combine(groupby(dropmissing(dm[:, [:mes, :prec]]), :mes), :prec => mean => :p)
    clim_mes = Dict(cm.mes .=> cm.p)
    clim_dia = Dict(m => clim_mes[m] / (m == 2 ? 28.25 : daysinmonth(Date(2001, m))) for m in 1:12)
    clim_anual = sum(values(clim_mes))
    umidos = collect(skipmissing(filter(x -> !ismissing(x) && x >= 1, dd.prec)))
    limiar95 = quantile(umidos, 0.95)

    anos = sort(unique(dd.ano_hidro))
    linhas = NamedTuple[]
    for a in anos
        s = dd[dd.ano_hidro .== a, :]
        inicio, fim = Date(a, MES_INICIO_HIDRO, 1), Date(a + 1, MES_INICIO_HIDRO, 1) - Day(1)
        cobertura = (minimum(s.data) == inicio) && (maximum(s.data) == fim)
        esperado_aus = sum((1 .- s.nh_prec ./ 24) .* [clim_dia[m] for m in s.mes])
        esperado_aus += sum(clim_dia[month(x)] for x in inicio:Day(1):fim if !(x in s.data); init = 0.0)
        fca = esperado_aus / clim_anual
        valido = cobertura && fca <= MAX_FCA
        p = collect(skipmissing(s.prec)); u = p[p .>= 1]
        push!(linhas, (; ano_hidro = a, rotulo = rotulo_hidro(a), prec_obs = sum(s.prec_obs),
            fca = fca, valido = valido,
            prec = valido ? sum(s.prec_obs) : missing,
            dias_validos = length(p),
            r1mm   = valido ? length(u) : missing,
            rx1day = valido ? maximum(p) : missing,
            sdii   = valido && !isempty(u) ? mean(u) : missing,
            r95_n  = valido ? count(>(limiar95), p) : missing))
    end
    dh = DataFrame(linhas)
    anual = agrega_anual(dm, clims, vars; chave = :ano_hidro)
    dh = leftjoin(dh, anual, on = :ano_hidro)
    sort!(dh, :ano_hidro)
    return dh, (clim_mes = clim_mes, clim_anual = clim_anual, limiar95 = limiar95)
end
