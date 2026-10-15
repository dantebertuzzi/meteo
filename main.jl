# =============================================================================
# main.jl — Análise climatológica da estação INMET A307 (Petrolina-PE), 2003–2024
#
# Uso:  julia --project=. main.jl
# Saídas: data/processed/*.csv, results/tabelas/*.csv, figures/*.{png,pdf}
# =============================================================================
include(joinpath(@__DIR__, "src", "00_config.jl"))
for f in ["01_importacao", "02_limpeza", "02b_preenchimento", "02c_ausentes", "03_descritiva", "04_temporal", "05_distribuicoes",
          "06_extremos", "07_correlacoes", "08_modelagem", "09_figuras"]
    include(joinpath(@__DIR__, "src", f * ".jl"))
end

function executa_analise()
    R = Dict{Symbol, Any}()

    # ------------------------------------------------------------------ 1. importação
    d, inventario = carrega_estacao(DIR_BRUTO, ESTACAO)
    mostra(inventario; titulo = "Inventário de arquivos")
    salva_tabela(inventario, "01_inventario_arquivos")
    @printf("\nRegistros horários: %d | %s → %s | timestamps duplicados: %d | horas faltantes na grade: %d\n",
            nrow(d), minimum(d.dt_utc), maximum(d.dt_utc), nrow(d) - length(unique(d.dt_utc)),
            length(minimum(d.dt_utc):Hour(1):maximum(d.dt_utc)) - length(unique(d.dt_utc)))

    # ------------------------------------------------------------------ 2. auditoria
    aud = auditoria_variaveis(d);        mostra(aud; titulo = "Auditoria por variável horária (antes do tratamento)")
    cons = checagens_consistencia(d);    mostra(cons; titulo = "Checagens de consistência interna")
    rad = auditoria_radiacao(d);         mostra(rad; titulo = "Radiação global por ano (kJ/m², horas com rad > 0)")
    travados = DataFrame(variavel = ["t", "ur", "p"],
                         seq_12h_constantes = [sequencias_constantes(d[!, v], 12) for v in (:t, :ur, :p)])
    mostra(travados; titulo = "Sequências ≥ 12 h com valor idêntico")
    R[:completude_prec] = completude_mensal(d, :prec)
    R[:completude_t]    = completude_mensal(d, :t)
    R[:radiacao] = rad
    foreach(((t, n),) -> salva_tabela(t, n), [(aud, "02_auditoria_variaveis"), (cons, "02_consistencia"),
            (rad, "02_radiacao_por_ano"), (R[:completude_prec], "02_completude_prec_mensal")])
    log_limites = aplica_limites!(d); mostra(log_limites; titulo = "Decisões de tratamento")
    salva_tabela(log_limites, "02_decisoes_tratamento")

    # ------------------------------------------------------------------ agregação
    dd = agrega_diario(d)
    if PREENCHER_COM_CONVENCIONAL
        conv = le_bdmep_diario(ARQ_CONV)
        ajuste, calib = calibra_convencional(dd, conv)
        mostra(calib; titulo = "Calibração A307 × convencional $ESTACAO_CONV (dias com as duas válidas)")
        salva_tabela(calib, "02_calibracao_convencional")
        aj = DataFrame(mes = 1:12; (v => [ajuste[v][m] for m in 1:12] for v in VARS_PREENCHE if v != :prec)...)
        mostra(aj; titulo = "Ajuste mensal aplicado à $ESTACAO_CONV (aditivo; vvel multiplicativo)")
        salva_tabela(aj, "02_ajuste_mensal_convencional")
        dd_antes = dd
        dd, log_pre = preenche_diario(dd, conv, ajuste)
        mostra(log_pre; titulo = "Dias preenchidos com a $ESTACAO_CONV por ano hidrológico")
        salva_tabela(log_pre, "02_preenchimento_por_ano")
        R[:calib] = calib; R[:preench] = log_pre
    else
        dd_antes = dd
    end
    R[:ausentes] = analisa_ausentes(d, dd_antes, dd)
    dm = agrega_mensal(dd)
    clims = anomalias_mensais!(dm, [VARS_MENSAIS; :prec])
    anual = agrega_anual(dm, clims, VARS_MENSAIS)
    itemp = indices_temperatura(dd, anual)
    anual = leftjoin(anual, itemp, on = :ano); sort!(anual, :ano)
    dh, climp = agrega_hidro(dd, dm, clims, VARS_MENSAIS)
    CSV.write(joinpath(DIR_PROC, "diario.csv"), dd)
    CSV.write(joinpath(DIR_PROC, "mensal.csv"), dm)
    CSV.write(joinpath(DIR_PROC, "anual_civil.csv"), anual)
    CSV.write(joinpath(DIR_PROC, "anual_hidrologico.csv"), dh)
    R[:d] = d; R[:dd] = dd; R[:dm] = dm; R[:anual] = anual; R[:dh] = dh; R[:climp] = climp; R[:clims] = clims

    val = dh[:, [:rotulo, :prec_obs, :fca, :valido, :dias_validos]]
    mostra(val; titulo = "Validade dos anos hidrológicos (FCA ≤ $(100MAX_FCA) %)")
    salva_tabela(val, "02_validade_anos_hidrologicos")
    mostra(anual[:, [:ano, :nmeses_tmed, :nmeses_tmax, :nmeses_tmin, :nmeses_ur]];
           titulo = "Meses válidos por ano civil (temperatura/UR)")
    @printf("\nClimatologia mensal de chuva (mm): %s\nTotal anual climatológico: %.1f mm | limiar P95 dia úmido: %.1f mm\n",
            join([@sprintf("%d:%.1f", m, climp.clim_mes[m]) for m in 1:12], "  "), climp.clim_anual, climp.limiar95)

    # ------------------------------------------------------------------ 3. descritiva
    desc_h = tabela_descritiva(dh, [:prec, :r1mm, :rx1day, :sdii, :r95_n])
    desc_a = tabela_descritiva(anual, [:tmed, :tmax, :tmin, :dtr, :ur, :urmin, :td, :p, :vvel, :txx, :tnn];
                               razao = Dict(v => false for v in (:tmed, :tmax, :tmin, :dtr, :td, :txx, :tnn)))
    desc_d = tabela_descritiva(dd, [:prec, :tmax, :tmin, :tmed, :ur];
                               razao = Dict(v => false for v in (:tmax, :tmin, :tmed)))
    desc_d.variavel .= desc_d.variavel .* "_diario"
    mostra(desc_h; titulo = "Descritiva — ano hidrológico (anos válidos)")
    mostra(desc_a; titulo = "Descritiva — ano civil (anos válidos); CV omitido p/ escala intervalar (°C)")
    mostra(desc_d; titulo = "Descritiva — escala diária (todos os dias válidos)")
    salva_tabela(vcat(desc_h, desc_a, desc_d), "03_descritiva")
    R[:desc] = vcat(desc_h, desc_a)

    # ------------------------------------------------------------------ 4. tendência
    series = [(dh.ano_hidro, dh.prec, "prec (mm/ano hidrológico)"), (dh.ano_hidro, dh.r1mm, "r1mm (dias ≥1 mm)"),
              (dh.ano_hidro, dh.rx1day, "rx1day (mm)"), (dh.ano_hidro, dh.sdii, "sdii (mm/dia)"),
              (anual.ano, anual.tmed, "tmed (°C)"), (anual.ano, anual.tmax, "tmax (°C)"),
              (anual.ano, anual.tmin, "tmin (°C)"), (anual.ano, anual.dtr, "dtr (°C)"),
              (anual.ano, anual.ur, "ur (%)"), (anual.ano, anual.td, "td (°C)"),
              (anual.ano, anual.p, "p (hPa)"), (anual.ano, anual.vvel, "vvel (m/s)"),
              (anual.ano, anual.txx, "txx (°C)"), (anual.ano, anual.tnn, "tnn (°C)")]
    tend = DataFrame([analise_tendencia(t, y; nome = n) for (t, y, n) in series])
    tend.q_mk = bh_ajuste(tend.p_mk)
    mostra(tend; titulo = "Tendência: Mann-Kendall, TFPW, Theil-Sen (por década) e Pettitt")
    salva_tabela(tend, "04_tendencias")
    R[:tend] = tend

    # anomalias mensais: MQO com EP HAC (≈ 12× mais pontos que a série anual)
    tend_m = DataFrame(map([:tmed, :tmax, :tmin, :dtr, :ur, :td, :vvel]) do v
        va = Symbol(v, :_anom)
        rg = regressao(dm, va, [:tempo]; hac = true, nome = String(va))
        c = rg.coefs[2, :]
        (variavel = String(va), n_meses = rg.n, tend_decada = 10c.estimativa,
         ic_inf = 10c.ic_inf, ic_sup = 10c.ic_sup, p_hac = c.p, R2 = rg.r2)
    end)
    tend_m.q_bh = bh_ajuste(tend_m.p_hac)
    mostra(tend_m; titulo = "Tendência em anomalias mensais — MQO com EP Newey-West")
    salva_tabela(tend_m, "04_tendencias_mensais_hac")
    R[:tend_m] = tend_m

    # ------------------------------------------------------------------ 5. distribuições
    dist = compara_distribuicoes(dh.prec)
    mostra(select(dist, Not(:modelo)); titulo = "Chuva anual (ano hidrológico): ajuste de distribuições")
    salva_tabela(select(dist, Not(:modelo)), "05_distribuicoes_prec")
    norm = normalidade(anual, [:tmed, :tmax, :tmin, :ur, :dtr])
    append!(norm, normalidade(dh, [:prec, :r1mm, :rx1day]))
    norm.p_bh = bh_ajuste(norm.p)
    mostra(norm; titulo = "Shapiro-Wilk (H₀: normalidade)")
    salva_tabela(norm, "05_normalidade")
    R[:dist] = dist

    # ------------------------------------------------------------------ 6. extremos
    dh.spi, gspi = spi_anual(dh.prec)
    dh.classe_spi = classe_spi.(dh.spi)
    anual.z_tmed = anomalia_padronizada(anual.tmed)
    anual.z_tmax = anomalia_padronizada(anual.tmax)
    anual.z_tmin = anomalia_padronizada(anual.tmin)
    lo, hi = cercas_tukey(dh.prec)
    mostra(dh[:, [:rotulo, :prec, :spi, :classe_spi, :r1mm, :rx1day]]; titulo = "SPI anual (Gamma: $(round.(params(gspi), sigdigits=4))); cercas de Tukey: [$(round(lo, digits=0)), $(round(hi, digits=0))] mm")
    mostra(anual[:, [:ano, :tmed, :z_tmed, :tmax, :z_tmax, :tmin, :z_tmin, :txx, :tnn]]; titulo = "Anomalias padronizadas de temperatura")
    salva_tabela(dh, "06_anos_hidrologicos_spi"); salva_tabela(anual, "06_anos_civis_anomalias")

    # POT / GPD na chuva diária
    dv = dropmissing(dd[:, [:data, :mes, :prec]])
    umid = dv.prec[dv.prec .>= 1]
    limiares = quantile(umid, 0.70:0.01:0.97)
    diag_u = diagnostico_limiar(dv.data, dv.prec, limiares)
    mostra(diag_u; titulo = "Diagnóstico de limiar (POT)")
    salva_tabela(diag_u, "06_pot_diagnostico_limiar")
    anos_ef = anos_efetivos(dd, climp.clim_mes, climp.clim_anual)
    u = quantile(umid, Q_LIMIAR_POT)
    picos, dpicos = desagrupa(dv.data, dv.prec, u)
    y = picos .- u; λ = length(y) / anos_ef
    gpd = ajusta_gpd(y)
    @printf("\nPOT: u = %.1f mm (P%d dos dias úmidos), %d picos desagrupados, anos efetivos = %.2f, λ = %.2f/ano\n",
            u, round(Int, 100Q_LIMIAR_POT), length(y), anos_ef, λ)
    @printf("GPD: σ = %.2f (EP %.2f), ξ = %.3f (EP %.3f)\n", gpd.σ, gpd.ep_σ, gpd.ξ, gpd.ep_ξ)
    rl = niveis_retorno(y, u, λ, anos_ef, [2, 5, 10, 20, 50])
    mostra(rl; titulo = "Níveis de retorno da chuva diária (IC 95 % bootstrap)")
    salva_tabela(rl, "06_niveis_retorno")
    top = sort(DataFrame(data = dpicos, prec = picos), :prec, rev = true)[1:min(10, end), :]
    mostra(top; titulo = "Maiores eventos diários desagrupados")
    R[:pot] = (u = u, y = y, λ = λ, gpd = gpd, rl = rl, anos = anos_ef, diag = diag_u, picos = picos, datas = dpicos)

    cont = DataFrame([
        tendencia_contagem(dh.ano_hidro, dh.r95_n, dh.dias_validos; familia = :poisson, nome = "R95 (dias > P95 úmido)"),
        tendencia_contagem(dh.ano_hidro, dh.r1mm, dh.dias_validos; familia = :poisson, nome = "R1mm (dias ≥ 1 mm)"),
        tendencia_contagem(anual.ano, anual.tx90p_n, anual.tx_dias; familia = :binomial, nome = "TX90p (dias Tmáx > P90)"),
        tendencia_contagem(anual.ano, anual.tn10p_n, anual.tn_dias; familia = :binomial, nome = "TN10p (dias Tmín < P10)")])
    cont.q_bh = bh_ajuste(cont.p)
    mostra(select(cont, Not(:modelo)); titulo = "Tendência em contagens de extremos (quase-Poisson/quase-binomial)")
    salva_tabela(select(cont, Not(:modelo)), "06_tendencia_contagens")
    R[:cont] = cont

    # ------------------------------------------------------------------ 7. correlações
    vars_h = [:prec, :r1mm, :rx1day, :tmed, :tmax, :tmin, :dtr, :ur, :td, :vvel]
    tc = tabela_correlacoes(dh, vars_h)
    mostra(sort(tc, :q_spearman)[1:min(20, end), :]; titulo = "Correlações anuais (ano hidrológico) — 20 menores q (Spearman)")
    salva_tabela(tc, "07_correlacoes_anuais")
    R[:cor_R], R[:cor_N] = matriz_correlacao(dh, vars_h); R[:cor_vars] = vars_h
    R[:cor_P] = [i == j ? 0.0 : p_correlacao(R[:cor_R][i, j], R[:cor_N][i, j]) for i in eachindex(vars_h), j in eachindex(vars_h)]
    cm = DataFrame([correlacao_mensal(dm, :tmax_anom, :tmin_anom), correlacao_mensal(dm, :tmax_anom, :ur_anom),
                    correlacao_mensal(dm, :tmin_anom, :td_anom), correlacao_mensal(dm, :ur_anom, :prec_anom; metodo = :spearman),
                    correlacao_mensal(dm, :tmax_anom, :prec_anom; metodo = :spearman), correlacao_mensal(dm, :dtr_anom, :ur_anom)])
    cm.q_bh = bh_ajuste(cm.p_nef)
    mostra(cm; titulo = "Correlações entre anomalias mensais (n efetivo de Bretherton)")
    salva_tabela(cm, "07_correlacoes_mensais")
    R[:cor_mensal] = cm

    regs = [regressao(dm, :tmax_anom, [:ur_anom]; hac = true, nome = "Tmáx ~ UR (mensal)"),
            regressao(dm, :tmin_anom, [:tmax_anom]; hac = true, nome = "Tmín ~ Tmáx (mensal)"),
            regressao(dh, :ur, [:prec]; nome = "UR ~ chuva (anual)"),
            regressao(dh, :tmax, [:prec]; nome = "Tmáx ~ chuva (anual)"),
            regressao(dm, :tmax_anom, [:ur_anom, :prec_anom]; hac = true, nome = "Tmáx ~ UR + chuva (mensal)")]
    tr = vcat([r.coefs for r in regs]...)
    tr_res = DataFrame(modelo = [r.coefs.modelo[1] for r in regs], n = [r.n for r in regs],
                       R2 = [r.r2 for r in regs], R2aj = [r.r2aj for r in regs], AIC = [r.aic for r in regs])
    mostra(tr; titulo = "Regressões entre variáveis — coeficientes")
    mostra(tr_res; titulo = "Regressões entre variáveis — ajuste")
    salva_tabela(tr, "07_regressoes_coef"); salva_tabela(tr_res, "07_regressoes_ajuste")
    R[:regs] = regs

    # ------------------------------------------------------------------ 8/9. modelagem e períodos
    mods = Dict{String, Any}()
    for (t, y, n) in [(dh.ano_hidro, dh.prec, "prec"), (anual.ano, anual.tmed, "tmed"), (anual.ano, anual.tmax, "tmax"),
                      (anual.ano, anual.tmin, "tmin"), (anual.ano, anual.ur, "ur"), (anual.ano, anual.dtr, "dtr")]
        mods[n] = modelo_tendencia(t, y; nome = n)
    end
    tm = DataFrame([mods[k].resumo for k in ["prec", "tmed", "tmax", "tmin", "ur", "dtr"]])
    tm.q_bh = bh_ajuste(tm.p_beta1)
    mostra(tm; titulo = "Modelos de tendência MQO (+ MQG-AR(1)) com diagnóstico")
    salva_tabela(tm, "08_modelos_tendencia")
    gam = tendencia_gamma(dh.ano_hidro, dh.prec)
    @printf("\nGLM Gamma (log) p/ chuva: %.1f %%/década [%.1f; %.1f], p = %.3f, AIC = %.1f (MQO: %.1f)\n",
            gam.pct_decada, gam.ic_inf, gam.ic_sup, gam.p, gam.AIC, mods["prec"].resumo.AIC)
    R[:mods] = mods; R[:tm] = tm; R[:gamma] = gam

    corte_h = 2012; corte_c = 2013
    per = DataFrame([
        compara_grupos(dh.prec[dh.ano_hidro .<= corte_h], dh.prec[dh.ano_hidro .> corte_h]; nome = "prec", rotulo_a = "2003/04–2012/13", rotulo_b = "2013/14–2023/24"),
        compara_grupos(dh.r1mm[dh.ano_hidro .<= corte_h], dh.r1mm[dh.ano_hidro .> corte_h]; nome = "r1mm", rotulo_a = "2003/04–2012/13", rotulo_b = "2013/14–2023/24"),
        [compara_grupos(anual[anual.ano .<= corte_c, v], anual[anual.ano .> corte_c, v]; nome = String(v), rotulo_a = "2003–2013", rotulo_b = "2014–2024")
         for v in (:tmed, :tmax, :tmin, :dtr, :ur)]...,
        compara_grupos(dh.tmax[.!ismissing.(dh.spi) .& (coalesce.(dh.spi, 0.0) .>= 0)],
                       dh.tmax[.!ismissing.(dh.spi) .& (coalesce.(dh.spi, 0.0) .< 0)];
                       nome = "tmax por SPI", rotulo_a = "SPI ≥ 0 (úmido)", rotulo_b = "SPI < 0 (seco)")])
    per.q_welch = bh_ajuste(per.p_welch)
    mostra(per; titulo = "Comparação entre períodos / grupos")
    salva_tabela(per, "08_comparacao_periodos")
    R[:per] = per; R[:cortes] = (corte_h, corte_c)
    R[:dh] = dh; R[:anual] = anual
    return R
end

R = executa_analise()
gera_figuras(R)
println("\nConcluído. Tabelas em results/tabelas, figuras em figures/.")
