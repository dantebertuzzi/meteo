# =============================================================================
# dados_longos.jl — série longa da estação convencional 82983 (1961–2025)
#
# Decisão D1: a chuva do artigo vem só da 82983 (um instrumento, sem
# preenchimento). A A307 entra como série de validação/homogeneidade.
#
# Critérios de completude (mesma lógica da análise principal):
# • Mês de chuva válido se a chuva *esperada* nos dias ausentes, pela
#   climatologia diária do mês, for ≤ max(1 mm, 5 % da normal do mês). Na
#   estação seca, alguns dias ausentes não comprometem o total; em março, sim.
# • Ano hidrológico (jul–jun) válido se a fração climatológica ausente FCA ≤ 5 %.
# • Temperatura mensal: regra da OMM (≤ 10 dias ausentes e ≤ 4 consecutivos).
# =============================================================================

const ARQ_LONGO = joinpath(RAIZ, "data", "raw", "bdmep", "dados_82983_D_1961-01-01_2025-06-30.csv")
const LAT_82983 = -9.3886111
const LON_82983 = -40.52333332
const ALT_82983 = 372.54

"""
Temperatura, UR e vento de 1962–1965 diferem 3–5 °C de Cabrobó e de Remanso
(contra poucos décimos após 1971): segmento inomogêneo, anterior à lacuna de
1966–1970. Essas variáveis só são usadas a partir de `INICIO_TEMP`; a chuva
desse segmento é compatível com as vizinhas e é mantida.
"""
const INICIO_TEMP = 1971

function carrega_longo(caminho = ARQ_LONGO)
    d = le_bdmep_diario(caminho)
    d.ano = year.(d.data); d.mes = month.(d.data); d.ano_hidro = ano_hidro.(d.data)
    for v in (:tmax, :tmin, :tmed, :ur, :urmin, :vvel)
        d[d.ano .< INICIO_TEMP, v] .= missing
    end
    return d
end

"Climatologia diária da chuva por mês (mm/dia), dos dias válidos de `ref`."
function clim_diaria(d; ref = (1991, 2020))
    s = d[(ref[1] .<= d.ano .<= ref[2]) .& .!ismissing.(d.prec), :]
    Dict(m => mean(Float64[x for x in s.prec[s.mes .== m]]) for m in 1:12)    # NaN se não houver dados
end

"Série mensal contínua (todas as combinações ano × mês do registro)."
function mensal_longo(d; ref = (1991, 2020))
    cd = clim_diaria(d; ref)
    linhas = NamedTuple[]
    for g in groupby(d, [:ano, :mes]; sort = true)
        a, m = g.ano[1], g.mes[1]; nd = daysinmonth(Date(a, m))
        nmiss = count(ismissing, g.prec) + (nd - nrow(g))
        esperado = nmiss * cd[m]; normal = cd[m] * nd
        okp = nmiss < nd && esperado <= max(1.0, 0.05normal)
        falta(v) = begin
            f = [ismissing(x) for x in g[!, v]]
            count(f) + (nd - nrow(g)) <= MAX_DIAS_FALTA && max_consecutivo(f) <= MAX_DIAS_CONSEC
        end
        med(v) = falta(v) ? mean(skipmissing(g[!, v])) : missing
        push!(linhas, (ano = a, mes = m, dias_aus_prec = nmiss, prec_esperada_aus = esperado,
                       prec = okp ? sum(skipmissing(g.prec); init = 0.0) : missing,
                       tmax = med(:tmax), tmin = med(:tmin), tmed = med(:tmed),
                       ur = med(:ur), vvel = med(:vvel)))
    end
    dm = DataFrame(linhas)
    dm.data = Date.(dm.ano, dm.mes, 1)
    dm.ano_hidro = [m >= MES_INICIO_HIDRO ? a : a - 1 for (a, m) in zip(dm.ano, dm.mes)]
    return dm
end

"Totais por ano hidrológico com FCA (fração da chuva climatológica caída em dias sem registro)."
function hidro_longo(d; ref = (1991, 2020))
    cd = clim_diaria(d; ref)
    anual = sum(cd[m] * daysinmonth(Date(2001, m)) for m in 1:12)
    out = combine(groupby(d, :ano_hidro; sort = true)) do s
        a = s.ano_hidro[1]
        dias = Date(a, MES_INICIO_HIDRO, 1):Day(1):Date(a + 1, MES_INICIO_HIDRO, 1) - Day(1)
        pos = Dict(s.data .=> s.prec)
        aus = sum(cd[month(x)] for x in dias if ismissing(get(pos, x, missing)); init = 0.0)
        fca = aus / anual
        tot = sum(skipmissing(s.prec); init = 0.0)
        (rotulo = rotulo_hidro(a), prec_obs = tot, dias_aus = count(x -> ismissing(get(pos, x, missing)), dias),
         fca = fca, valido = fca <= MAX_FCA, prec = fca <= MAX_FCA ? tot : missing)
    end
    return out, anual
end
