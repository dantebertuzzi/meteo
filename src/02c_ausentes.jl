# =============================================================================
# 02c_ausentes.jl — diagnóstico dos valores ausentes com MissingPatterns.jl
#
# Três perguntas:
#   1. Onde e quanto falta? (mapa de calor por ano, horário e diário)
#   2. As variáveis falham juntas? (padrões únicos e ϕ entre máscaras) — se a
#      falha é do registrador, todas somem ao mesmo tempo e a ausência depende
#      do tempo, não do valor medido.
#   3. O que o preenchimento com a 82983 resolveu? (diff antes/depois)
# Saídas: tabelas em results/tabelas/02_ausentes_*.csv e mapas HTML em figures/.
# =============================================================================
using MissingPatterns

const VARS_AUSENTES_DIA = [:prec, :tmed, :tmax, :tmin, :ur, :urmin, :td, :p, :vvel]

"Padrões únicos como DataFrame plano: uma coluna Bool por variável (true = ausente)."
padroes_df(tbl) = DataFrame([(; r.pattern..., n_ausentes = r.nmissing, n = r.n, pct = r.pct)
                             for r in missingpatternstats(tbl)])

function analisa_ausentes(d::DataFrame, dd_antes::DataFrame, dd::DataFrame)
    # --- horário (A307 bruta, após limites físicos; rad fica de fora: foi excluída)
    vh = filter(!=(:rad), VARS_NUM)
    H = select(d, :dt_local, vh)
    println("\n=== Ausentes — A307 horária (resumo por variável) ===")
    missingsummary(select(H, vh); color = :never)
    println("\n=== Ausentes — A307 horária (padrões entre variáveis) ===")
    missingpatterns(select(H, vh); max_patterns = 10)
    missinghtml(joinpath(DIR_FIG, "ausentes_horario_a307.html"), H; by = :dt_local, period = :year,
                title = "A307 — ausentes por ano (dados horários)", emphasis = :missing)

    # --- diário antes × depois do preenchimento, na mesma grade de datas
    A = leftjoin(DataFrame(data = dd.data), select(dd_antes, :data, VARS_AUSENTES_DIA), on = :data)
    sort!(A, :data)
    B = select(dd, :data, VARS_AUSENTES_DIA)
    println("\n=== Ausentes — série diária ANTES do preenchimento ===")
    missingsummary(select(A, VARS_AUSENTES_DIA); color = :never)
    println("\n=== Ausentes — série diária DEPOIS do preenchimento ===")
    missingsummary(select(B, VARS_AUSENTES_DIA); color = :never)
    missingpatterns(select(B, VARS_AUSENTES_DIA); max_patterns = 10)
    println("\n=== Diferença antes → depois (− = resolvido, + = introduzido) ===")
    plotmissingdiff(select(A, VARS_AUSENTES_DIA), select(B, VARS_AUSENTES_DIA); color = :never)
    println("\n=== Mapa de ausentes por ano — série diária final ===")
    plotmissing(B; by = :data, period = :year, color = :never, layout = :classic)

    for (tbl, nome, tit) in ((A, "antes", "antes do preenchimento"), (B, "depois", "após preenchimento com a $ESTACAO_CONV"))
        missinghtml(joinpath(DIR_FIG, "ausentes_diario_$nome.html"), tbl; by = :data, period = :year,
                    title = "Série diária — ausentes por ano, $tit", emphasis = :missing)
    end

    # --- tabelas
    st = leftjoin(rename(DataFrame(missingstats(select(A, VARS_AUSENTES_DIA)))[:, [:column, :nmissing, :pct]],
                         :nmissing => :ausentes_antes, :pct => :pct_antes),
                  rename(DataFrame(missingstats(select(B, VARS_AUSENTES_DIA)))[:, [:column, :nmissing, :pct]],
                         :nmissing => :ausentes_depois, :pct => :pct_depois), on = :column)
    rename!(st, :column => :variavel)
    st.variavel = String.(st.variavel)
    st.resolvidos = st.ausentes_antes .- st.ausentes_depois
    pares = DataFrame(missingpairstats(select(B, VARS_AUSENTES_DIA)))
    pares.a = String.(pares.a); pares.b = String.(pares.b)
    salva_tabela(st, "02_ausentes_por_variavel")
    salva_tabela(padroes_df(select(A, VARS_AUSENTES_DIA)), "02_ausentes_padroes_antes")
    salva_tabela(padroes_df(select(B, VARS_AUSENTES_DIA)), "02_ausentes_padroes_depois")
    salva_tabela(sort(pares, :phi, rev = true), "02_ausentes_coocorrencia")
    salva_tabela(padroes_df(select(H, vh)), "02_ausentes_padroes_horario")
    return st
end
