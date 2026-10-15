# =============================================================================
# 02b_preenchimento.jl — preenchimento das lacunas diárias da A307 com a estação
# convencional INMET 82983 (Petrolina), obtida no BDMEP.
#
# A A307 continua sendo a série principal: um dia só recebe o valor da 82983
# quando o valor diário da A307 é `missing` (falha do registrador ou dia
# incompleto). Cada variável é tratada isoladamente e marcada em `<v>_pre`.
#
# As duas estações usam instrumentos e protocolos diferentes (observador 3×/dia,
# psicrômetro, pluviômetro lido às 12 UTC), então os valores da 82983 são
# ajustados antes de entrar na série:
#   • tmed, tmax, tmin, ur, urmin: viés aditivo por mês do ano, Δ(mês) = média(A307 − 82983)
#   • vvel: fator multiplicativo por mês (anemômetros em alturas/instrumentos distintos)
#   • prec: sem ajuste — o total da 82983 difere < 2 % no período comum, e um
#     fator único distorceria os dias de chuva forte.
# Os ajustes são estimados só nos dias em que as duas estações têm valor válido.
# =============================================================================

const PREENCHER_COM_CONVENCIONAL = true
const ARQ_CONV = joinpath(RAIZ, "data", "raw", "bdmep", "dados_82983_D_2003-01-01_2025-06-30.csv")
const ESTACAO_CONV = "82983"

const MAPA_CONV = [
    "PRECIPITACAO TOTAL"             => :prec,
    "TEMPERATURA MAXIMA"             => :tmax,
    "TEMPERATURA MEDIA COMPENSADA"   => :tmed,
    "TEMPERATURA MINIMA"             => :tmin,
    "UMIDADE RELATIVA DO AR, MEDIA"  => :ur,
    "UMIDADE RELATIVA DO AR, MINIMA" => :urmin,
    "VENTO, VELOCIDADE MEDIA"        => :vvel,
]
const VARS_PREENCHE = last.(MAPA_CONV)
# Limites para dados DIÁRIOS. A chuva não pode reusar LIMITES[:prec] (150 mm, um limite
# HORÁRIO): totais diários acima de 150 mm são raros no semiárido, mas ocorrem (82983,
# 01/01/1978: 151,3 mm; Cabrobó tem 5 dias assim). 400 mm fica acima de qualquer registro
# diário da região e só elimina erro grosseiro de digitação.
const PREC_DIARIA_MAX = 400.0
const LIMITE_CONV = Dict(:prec => (0.0, PREC_DIARIA_MAX), :tmax => LIMITES[:tmax], :tmed => LIMITES[:t],
                         :tmin => LIMITES[:tmin], :ur => LIMITES[:ur], :urmin => LIMITES[:urmin],
                         :vvel => LIMITES[:vvel])

"""
    le_bdmep_diario(caminho) -> DataFrame

Lê o CSV diário do BDMEP: 10 linhas de metadados, separador ';', vírgula decimal,
ausentes como `null`. Valores fora do limite físico viram `missing`.
"""
function le_bdmep_diario(caminho::AbstractString)
    linhas = readlines(caminho)
    i = findfirst(l -> startswith(l, "Data Medicao"), linhas)
    isnothing(i) && error("Cabeçalho 'Data Medicao' não encontrado em $(basename(caminho))")
    df = CSV.read(IOBuffer(join(linhas[i:end], '\n')), DataFrame; delim = ';', header = 1,
                  types = String, missingstring = ["", "null"], silencewarnings = true)
    out = DataFrame(data = Date.(df[!, 1]))
    for nome in names(df)[2:end]
        k = findfirst(p -> startswith(uppercase(nome), p.first), MAPA_CONV)
        isnothing(k) && continue
        v = MAPA_CONV[k].second
        lo, hi = LIMITE_CONV[v]
        out[!, v] = [ismissing(x) || x < lo || x > hi ? missing : x for x in parse_num.(df[!, nome])]
    end
    faltando = setdiff(VARS_PREENCHE, propertynames(out))
    isempty(faltando) || error("Colunas não encontradas em $(basename(caminho)): $faltando")
    return out
end

const AJUSTE_MULT = Set([:vvel])

"""
    calibra_convencional(dd, conv) -> (ajuste, resumo)

Estima o ajuste mês a mês da 82983 para a escala da A307 nos dias com as duas
estações válidas. `ajuste[v][mês]` é aditivo (°C, %) ou multiplicativo (vvel).
`resumo` traz, por variável, n, correlação diária, viés médio e RMSE antes e
depois do ajuste.
"""
function calibra_convencional(dd::DataFrame, conv::DataFrame)
    j = innerjoin(dd[:, [:data, :mes, VARS_PREENCHE...]],
                  rename(conv, (v => Symbol(v, :_c) for v in VARS_PREENCHE)...), on = :data)
    ajuste = Dict{Symbol, Dict{Int, Float64}}()
    linhas = NamedTuple[]
    for v in VARS_PREENCHE
        vc = Symbol(v, :_c)
        s = dropmissing(j[:, [:mes, v, vc]])
        mult = v in AJUSTE_MULT
        ajuste[v] = Dict(m => begin
                g = s[s.mes .== m, :]
                v == :prec ? (mult ? 1.0 : 0.0) :
                mult ? mean(g[!, v]) / mean(g[!, vc]) : mean(g[!, v] .- g[!, vc])
            end for m in 1:12)
        aplicado = [aplica_ajuste(v, x, m, ajuste) for (x, m) in zip(s[!, vc], s.mes)]
        push!(linhas, (variavel = String(v), tipo_ajuste = v == :prec ? "nenhum" : mult ? "multiplicativo" : "aditivo",
                       n_dias = nrow(s), r_diario = cor(s[!, v], s[!, vc]),
                       vies_antes = mean(s[!, vc] .- s[!, v]), rmse_antes = sqrt(mean((s[!, vc] .- s[!, v]) .^ 2)),
                       vies_depois = mean(aplicado .- s[!, v]), rmse_depois = sqrt(mean((aplicado .- s[!, v]) .^ 2)),
                       total_a307 = v == :prec ? sum(s[!, v]) : missing,
                       total_82983 = v == :prec ? sum(s[!, vc]) : missing))
    end
    return ajuste, DataFrame(linhas)
end

aplica_ajuste(v, x, m, ajuste) = ismissing(x) ? missing :
    v in AJUSTE_MULT ? x * ajuste[v][m] : x + ajuste[v][m]

"""
    preenche_diario(dd, conv, ajuste) -> (dd_preenchido, log_por_ano)

Estende a série diária até o fim do registro da 82983 e substitui cada valor
diário ausente da A307 pelo valor ajustado da 82983. Na chuva, o dia preenchido
passa a contar como 24 h observadas (`nh_prec = 24`), de modo que os critérios
de completude mensal e a FCA o tratam como dia completo.
"""
function preenche_diario(dd::DataFrame, conv::DataFrame, ajuste)
    todas = minimum(dd.data):Day(1):max(maximum(dd.data), maximum(conv.data))
    novas = setdiff(todas, dd.data)
    if !isempty(novas)
        extra = DataFrame(data = novas, nh = 0, nh_prec = 0, prec_obs = 0.0)
        dd = vcat(dd, extra; cols = :union)
    end
    dd = sort(dd, :data)
    dd.ano = year.(dd.data); dd.mes = month.(dd.data); dd.ano_hidro = ano_hidro.(dd.data)

    pos = Dict(conv.data .=> 1:nrow(conv))
    for v in VARS_PREENCHE
        dd[!, v] = convert(Vector{Union{Missing, Float64}}, dd[!, v])
        flag = falses(nrow(dd))
        for i in 1:nrow(dd)
            ismissing(dd[i, v]) || continue
            k = get(pos, dd.data[i], nothing)
            isnothing(k) && continue
            x = aplica_ajuste(v, conv[k, v], dd.mes[i], ajuste)
            ismissing(x) && continue
            dd[i, v] = x; flag[i] = true
            if v == :prec
                dd.prec_obs[i] = x; dd.nh_prec[i] = 24
            end
        end
        dd[!, Symbol(v, :_pre)] = flag
    end
    dd.dtr = dd.tmax .- dd.tmin

    log = combine(groupby(dd, :ano_hidro),
        nrow => :dias,
        [Symbol(v, :_pre) => sum => Symbol(:pre_, v) for v in VARS_PREENCHE]...,
        :prec => (x -> count(ismissing, x)) => :dias_sem_chuva_registrada)
    log.rotulo = rotulo_hidro.(log.ano_hidro)
    select!(log, :rotulo, Not(:rotulo))
    return dd, log
end
