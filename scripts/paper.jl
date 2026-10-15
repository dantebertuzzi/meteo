# =============================================================================
# scripts/paper.jl — regenera todos os números (paper/resultados.tex), tabelas
# (results/paper/) e figuras (paper/figuras/) do artigo.
#
# Pré-requisito: `julia --project=. main.jl` (gera data/processed/).
# Uso:          julia --project=. -t auto scripts/paper.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "paper", "carrega.jl"))

const REF = (1991, 2020)            # período de referência (normal OMM) para SPI/SPEI
const HORIZONTES = (5, 10, 20, 30)

titulo(s) = println("\n", "="^length(s), "\n", s, "\n", "="^length(s))

# ============================================================================= dados
function bloco_dados()
    titulo("Série longa 82983")
    d = carrega_longo()
    dm = mensal_longo(d; ref = REF)
    dh, clim_anual = hidro_longo(d; ref = REF)
    v = dh[dh.valido, :]
    @printf("Registro diário: %s → %s | anos hidrológicos válidos: %d de %d | normal 1991–2020: %.1f mm\n",
            minimum(d.data), maximum(d.data), nrow(v), nrow(dh), clim_anual)
    registra!("NanosValidos", nrow(v)); registra!("NanosHidro", nrow(dh))
    registra!("PrimeiroAnoValido", v.rotulo[1]); registra!("UltimoAnoValido", v.rotulo[end])
    registra!("NormalChuva", clim_anual; dig = 0)
    registra!("MediaChuvaValidos", mean(v.prec); dig = 0); registra!("DPChuvaValidos", std(v.prec); dig = 0)
    salva_ptab(dh, "anos_hidrologicos_82983")
    salva_ptab(dm, "mensal_82983")
    return (d = d, dm = dm, dh = dh, clim_anual = clim_anual)
end

# ============================================================================= homogeneidade
"""
A307 *sem* os dias preenchidos × 82983, só nos dias com as duas válidas.
Chuva: razão anual Σ A307 / Σ 82983 nos dias comuns. Temperatura: diferença
média anual. Anos com < 300 dias comuns (chuva) ou < 200 (temperatura) saem.
"""
function bloco_homogeneidade(L)
    titulo("Homogeneidade relativa A307 × 82983")
    a = CSV.read(joinpath(DIR_PROC, "diario.csv"), DataFrame)
    for v in (:prec, :tmax, :tmin)
        pre = Symbol(v, :_pre)
        a[!, v] = [p ? missing : x for (x, p) in zip(a[!, v], a[!, pre])]
    end
    j = innerjoin(select(a, :data, :prec => :prec_a, :tmax => :tmax_a, :tmin => :tmin_a),
                  select(L.d, :data, :ano_hidro, :ano, :prec, :tmax, :tmin), on = :data)
    linhas = NamedTuple[]
    for g in groupby(j, :ano_hidro; sort = true)
        sp = dropmissing(g[:, [:prec_a, :prec]])
        st = dropmissing(g[:, [:tmax_a, :tmax, :tmin_a, :tmin]])
        push!(linhas, (ano_hidro = g.ano_hidro[1], dias_prec = nrow(sp),
                       razao_prec = nrow(sp) >= 300 && sum(sp.prec) > 0 ? sum(sp.prec_a) / sum(sp.prec) : missing,
                       prec_a307 = nrow(sp) >= 300 ? sum(sp.prec_a) : missing,
                       prec_82983 = nrow(sp) >= 300 ? sum(sp.prec) : missing,
                       dias_t = nrow(st),
                       dif_tmax = nrow(st) >= 200 ? mean(st.tmax_a .- st.tmax) : missing,
                       dif_tmin = nrow(st) >= 200 ? mean(st.tmin_a .- st.tmin) : missing))
    end
    h = DataFrame(linhas)
    res = DataFrame(map([(:razao_prec, "chuva (razão A307/82983)"), (:dif_tmax, "Tmáx (A307 − 82983, °C)"),
                         (:dif_tmin, "Tmín (A307 − 82983, °C)")]) do (v, nome)
        s = dropmissing(h[:, [:ano_hidro, v]]); t = Float64.(s.ano_hidro); y = Float64.(s[!, v])
        sn = snht_teste(t, y); pt = pettitt(t, y)
        (serie = nome, n = length(y), media = mean(y), dp = std(y), min = minimum(y), max = maximum(y),
         snht_T0 = sn.T0, snht_T95 = sn.T95, p_snht = sn.p, quebra_snht_apos = rotulo_hidro(Int(sn.ultimo_antes)),
         p_pettitt = pt.p, quebra_pettitt_apos = rotulo_hidro(Int(pt.ultimo_antes)),
         media_antes = sn.media_antes, media_depois = sn.media_depois)
    end)
    res.q_bh = bh_ajuste(vcat(res.p_snht))
    mostra(h; titulo = "Série anual de comparação (dias comuns)")
    mostra(res; titulo = "SNHT (p por Monte Carlo) e Pettitt")
    salva_ptab(h, "homogeneidade_anual"); salva_ptab(res, "homogeneidade_testes")
    r = res[1, :]
    registra!("RazaoChuvaMedia", r.media; dig = 3); registra!("RazaoChuvaMin", r.min; dig = 2); registra!("RazaoChuvaMax", r.max; dig = 2)
    registra_p!("PSnhtChuva", r.p_snht); registra!("QuebraSnhtChuva", r.quebra_snht_apos)
    registra_p!("PPettittChuva", r.p_pettitt)
    registra!("RazaoChuvaAntes", r.media_antes; dig = 2); registra!("RazaoChuvaDepois", r.media_depois; dig = 2)
    for (i, nm) in ((2, "Tmax"), (3, "Tmin"))
        registra!("DifMedia" * nm, res.media[i]; dig = 2); registra_p!("PSnht" * nm, res.p_snht[i])
        registra_p!("PPettitt" * nm, res.p_pettitt[i])
    end
    i1213 = findfirst(==(2012), h.ano_hidro)
    !isnothing(i1213) && registra!("RazaoChuvaDozeTreze", h.razao_prec[i1213]; dig = 2)
    return h, res
end

# ============================================================================= homogeneidade regional
const VIZINHAS = Dict("82886" => "Cabrobó", "82979" => "Remanso", "82986" => "Paulo Afonso", "82753" => "Ouricuri",
                      "82987" => "Petrolândia", "82985" => "Caraíba", "83088" => "Senhor do Bonfim", "82887" => "Floresta")

"Homogeneidade absoluta-relativa da 82983 contra referência ponderada das convencionais vizinhas."
function bloco_vizinhas(L)
    titulo("Homogeneidade da 82983 com estações vizinhas")
    dir = joinpath(RAIZ, "data", "raw", "bdmep")
    arqs = Dict(c => joinpath(dir, f) for f in readdir(dir) for c in keys(VIZINHAS) if startswith(f, "dados_$(c)_D_"))
    if isempty(arqs)
        println("Arquivos das vizinhas ausentes em $dir — bloco pulado.")
        return nothing
    end
    viz = Dict{String, Dict{Int, Float64}}(); resumo = NamedTuple[]
    for (c, f) in arqs
        d = le_bdmep_diario(f); d.ano = year.(d.data); d.mes = month.(d.data); d.ano_hidro = ano_hidro.(d.data)
        # estações fechadas antes de 1991 não têm dados na normal 1991–2020: a FCA usa a climatologia do próprio registro
        nref = count(i -> !ismissing(d.prec[i]) && REF[1] <= d.ano[i] <= REF[2], 1:nrow(d))
        refv = nref >= 3650 ? REF : (minimum(d.ano), maximum(d.ano))
        h, normal = hidro_longo(d; ref = refv)
        hv = h[h.valido, :]
        viz[VIZINHAS[c]] = Dict(zip(hv.ano_hidro, Float64.(hv.prec)))
        tv = nrow(hv) >= 15 ? analise_tendencia(hv.ano_hidro, hv.prec; nome = VIZINHAS[c]) : nothing
        push!(resumo, (estacao = VIZINHAS[c], codigo = c, anos_validos = nrow(hv), normal = normal,
                       sen_decada = isnothing(tv) ? missing : tv.sen_decada, p_mk = isnothing(tv) ? missing : tv.p_mk,
                       p_tfpw = isnothing(tv) ? missing : tv.p_mk_tfpw,
                       quebra_pettitt = isnothing(tv) ? missing : tv.mudanca_apos, p_pettitt = isnothing(tv) ? missing : tv.p_pettitt))
    end
    rs = DataFrame(resumo); rs.q_mk = bh_ajuste(coalesce.(rs.p_mk, 1.0))
    mostra(rs; titulo = "Tendência da chuva anual nas vizinhas (mesmo critério de validade)")
    salva_ptab(rs, "vizinhas_tendencia")
    cand = Dict(zip(L.dh.ano_hidro[L.dh.valido], Float64.(L.dh.prec[L.dh.valido])))
    Q, info = serie_referencia(cand, viz)
    mostra(info; titulo = "Vizinhas: correlação das primeiras diferenças com a 82983")
    salva_ptab(info, "vizinhas_correlacao"); salva_ptab(Q, "vizinhas_serie_Q")
    sn = snht_teste(Float64.(Q.ano), Q.Q); pt = pettitt(Float64.(Q.ano), Q.Q)
    tq = analise_tendencia(Q.ano, Q.Q; nome = "Q = 82983/referência")
    @printf("Série Q (n = %d): SNHT T₀ = %.2f (T95 = %.2f, p = %.3f, quebra após %s, Q %.2f → %.2f) | Pettitt p = %.3f | tendência de Q: %.3f/década (p = %.3f)\n",
            nrow(Q), sn.T0, sn.T95, sn.p, rotulo_hidro(Int(sn.ultimo_antes)), sn.media_antes, sn.media_depois, pt.p, tq.sen_decada, tq.p_mk)
    registra!("NvizinhasUsadas", count(info.usada)); registra!("NanosQ", nrow(Q))
    registra_p!("PSnhtQ", sn.p); registra!("QuebraQ", rotulo_hidro(Int(sn.ultimo_antes)))
    registra!("QAntes", sn.media_antes; dig = 2); registra!("QDepois", sn.media_depois; dig = 2)
    registra_p!("PPettittQ", pt.p); registra_p!("PTendQ", tq.p_mk)
    registra!("NvizinhasTendNeg", count(r -> !ismissing(r.sen_decada) && r.sen_decada < 0, eachrow(rs)))
    registra!("NvizinhasTendSig", count(<(0.05), rs.q_mk))
    # tendência da referência regional (sinal climático comum)
    tr = analise_tendencia(Q.ano, Q.referencia; nome = "referência regional")
    registra_ic!("SenReferencia", tr.sen_decada, tr.sen_ic_inf, tr.sen_ic_sup; dig = 0); registra_p!("PMKReferencia", tr.p_mk)
    mostra(DataFrame([tq, tr]); titulo = "Tendência: razão Q (inomogeneidade) × referência regional (clima)")
    registra_ic!("SenQ", 100tq.sen_decada, 100tq.sen_ic_inf, 100tq.sen_ic_sup; dig = 1)

    # par a par: log(82983/vizinha) — a referência composta muda de composição quando
    # estações fecham (S. do Bonfim 2011, P. Afonso 2015); pares fixos evitam esse artefato.
    pares = DataFrame(map([r for r in eachrow(info) if r.usada]) do r
        V = viz[r.estacao]; an = sort([a for a in keys(cand) if haskey(V, a)])
        lq = [log(cand[a] / V[a]) for a in an]
        t = analise_tendencia(an, lq; nome = r.estacao); s = snht_teste(Float64.(an), lq)
        (vizinha = r.estacao, rho_dif = r.rho_dif, n = length(an), inicio = rotulo_hidro(an[1]), fim = rotulo_hidro(an[end]),
         sen_pct_decada = 100(exp(t.sen_decada) - 1), ic_inf = 100(exp(t.sen_ic_inf) - 1), ic_sup = 100(exp(t.sen_ic_sup) - 1),
         p_mk = t.p_mk, p_snht = s.p, quebra_snht = rotulo_hidro(Int(s.ultimo_antes)))
    end)
    pares.q_snht = bh_ajuste(pares.p_snht)
    mostra(pares; titulo = "Razões par a par log(82983/vizinha): tendência relativa e quebra")
    salva_ptab(pares, "vizinhas_pares")
    registra!("NparesQuebraSig", count(<(0.05), pares.q_snht)); registra!("Npares", nrow(pares))
    registra!("PSnhtParesMin", minimum(pares.p_snht); dig = 2)
    return (Q = Q, info = info, rs = rs, snht = sn, pares = pares)
end

# ============================================================================= P1: evento de referência
function bloco_referencia(L)
    titulo("Evento de referência")
    v = sort(L.dh[L.dh.valido, :], :prec)
    mostra(first(v[:, [:rotulo, :prec, :fca]], 8); titulo = "Anos hidrológicos mais secos (82983, válidos)")
    x1213 = only(L.dh.prec[L.dh.ano_hidro .== 2012])
    xmin = v.prec[1]; rot_min = v.rotulo[1]
    registra!("ChuvaDozeTreze", x1213; dig = 1); registra!("ChuvaMinima", xmin; dig = 1); registra!("AnoMinimo", rot_min)
    for (i, r) in enumerate(eachrow(first(v, 5)))
        registra!("RankAno" * string(i), r.rotulo); registra!("RankChuva" * string(i), r.prec; dig = 1)
    end
    rank1213 = findfirst(==("2012/13"), v.rotulo)
    registra!("RankDozeTreze", rank1213)
    registra!("RankDozeTrezeEn", ("driest", "second driest", "third driest", "fourth driest", "fifth driest")[min(rank1213, 5)])
    registra!("PctNormalDozeTreze", 100x1213 / L.clim_anual; dig = 0)
    return (x1213 = x1213, xmin = xmin, rot_min = rot_min, v = v)
end

# ============================================================================= P1: distribuições e período de retorno
function bloco_retorno(L, E)
    titulo("Ajuste de distribuições e período de retorno (82983)")
    x = Float64.(E.v.prec); n = length(x)
    cmp = compara_distribuicoes(x)
    cmp.p_x1213 = [cdf(m, E.x1213) for m in cmp.modelo]
    cmp.T_x1213 = 1 ./ cmp.p_x1213
    cmp.p_xmin = [cdf(m, E.xmin) for m in cmp.modelo]
    cmp.T_xmin = 1 ./ cmp.p_xmin
    mostra(select(cmp, Not(:modelo)); titulo = "Candidatas (n = $n)")
    salva_ptab(select(cmp, Not(:modelo)), "distribuicoes_chuva_anual")
    p_med = sum(cmp.peso_akaike .* cmp.p_x1213)
    registra!("TMediaAkaike", 1 / p_med; sig = 2)
    # IC da média ponderada: reamostra os anos, reajusta as 4 candidatas e recalcula os pesos
    rng = MersenneTwister(SEMENTE); pm = Float64[]
    for _ in 1:2000
        xb = x[rand(rng, 1:n, n)]
        fits = [try ajusta_mle(D, xb) catch; nothing end for D in CANDIDATAS]
        any(isnothing, fits) && continue
        aic = [2length(params(f)) - 2loglikelihood(f, xb) for f in fits]
        w = exp.(-(aic .- minimum(aic)) ./ 2); w ./= sum(w)
        push!(pm, sum(w .* [cdf(f, E.x1213) for f in fits]))
    end
    # estimativa pontual = amostra completa; o bootstrap fornece só o intervalo (percentis)
    registra_ic!("TMediaAkaikeB", 1 / p_med, 1 / quantile(pm, 0.975), 1 / quantile(pm, 0.025); sig = 2)
    for nh in HORIZONTES
        registra_ic!("ProbMediaAkaikeN" * string(nh), 100prob_ao_menos_um(p_med, nh),
                     100quantile(prob_ao_menos_um.(pm, nh), 0.025), 100quantile(prob_ao_menos_um.(pm, nh), 0.975); dig = 0)
    end
    registra!("TModelosMin", minimum(cmp.T_x1213); sig = 2); registra!("TModelosMax", maximum(cmp.T_x1213); sig = 2)
    # faixa só entre as candidatas não rejeitadas pelo Anderson–Darling (p ≥ 0,05)
    aceitas = cmp.p_AD_boot .>= 0.05
    registra!("TModelosMaxAceitos", maximum(cmp.T_x1213[aceitas]); sig = 2)
    registra!("DistRejeitadasAD", join(lowercase.(cmp.distribuicao[.!aceitas]), ", "))
    registra!("NDistRejeitadasAD", count(!, aceitas))
    g = only(cmp.modelo[cmp.distribuicao .== "Gamma"]); ig = findfirst(==("Gamma"), cmp.distribuicao)
    registra!("GammaForma", shape(g); dig = 2); registra!("GammaEscala", scale(g); dig = 1)
    registra!("GammaPesoAkaike", cmp.peso_akaike[ig]; dig = 2); registra!("GammaPAD", cmp.p_AD_boot[ig]; dig = 2)
    registra!("MelhorAIC", cmp.distribuicao[1])
    for r in eachrow(cmp)
        registra!("T" * r.distribuicao, r.T_x1213; sig = 2)
    end

    # empírico
    k1213 = count(<=(E.x1213), x); kmin = count(<=(E.xmin), x)
    registra!("NanosAteDozeTreze", k1213)
    registra!("TEmpiricoDozeTreze", (n + 1) / k1213; sig = 2)

    out = Dict{String, Any}()
    for (rot, x0) in (("DozeTreze", E.x1213), ("Minimo", E.xmin))
        pr = perfil_p_gamma(x, x0)
        bs = bootstrap_p(Gamma, x, x0; B = 10_000)
        rb = resumo_replicas(bs)
        @printf("\n[%s] x₀ = %.1f mm | Gamma: p̂ = %.4f, T̂ = %.1f anos\n", rot, x0, pr.p, 1 / pr.p)
        @printf("  perfilada: p ∈ [%.4f, %.4f] → T ∈ [%.1f, %.1f]\n", pr.inf, pr.sup, 1 / pr.sup, 1 / pr.inf)
        @printf("  bootstrap: p ∈ [%.4f, %.4f] → T ∈ [%.1f, %.1f]\n", rb.p[1], rb.p[3], rb.T[1], rb.T[3])
        registra!("P" * rot, pr.p; dig = 3)
        registra_ic!("T" * rot, 1 / pr.p, 1 / pr.sup, 1 / pr.inf; sig = 2)          # perfilada
        registra_ic!("TBoot" * rot, rb.T[2], rb.T[1], rb.T[3]; sig = 2)
        for nh in HORIZONTES
            pn = prob_ao_menos_um(pr.p, nh); q = rb.prob_n[nh]
            registra_ic!("Prob" * rot * "N" * string(nh), 100pn, 100q[1], 100q[3]; dig = 0)
        end
        out[rot] = (perfil = pr, boot = bs, resumo = rb)
    end
    return (cmp = cmp, n = n, out = out)
end

# ============================================================================= P1: bayes
function bloco_bayes(L, E)
    titulo("Abordagem bayesiana (Turing/NUTS)")
    x = Float64.(E.v.prec)
    res = DataFrame()
    draws = Dict{String, Any}()
    for (cen, sk, sm) in (("base", 1.0, 1.0), ("vaga", 2.0, 2.0))
        ps, ch = posterior_p(x, E.x1213; sd_k = sk, sd_mu = sm)
        rk = diagnostico_cadeias(ch, :logk); rm = diagnostico_cadeias(ch, :logmu)
        rb = resumo_replicas(ps)
        push!(res, (cenario = cen, rhat_logk = rk, rhat_logmu = rm, p_med = rb.p[2], T_med = rb.T[2],
                     T_inf = rb.T[1], T_sup = rb.T[3],
                     p30_med = rb.prob_n[30][2], p30_inf = rb.prob_n[30][1], p30_sup = rb.prob_n[30][3]))
        draws[cen] = ps
        if cen == "base"
            registra_ic!("TBayes", rb.T[2], rb.T[1], rb.T[3]; sig = 2)
            for nh in HORIZONTES
                q = rb.prob_n[nh]; registra_ic!("ProbBayesN" * string(nh), 100q[2], 100q[1], 100q[3]; dig = 0)
            end
            registra!("RhatMax", max(rk, rm); dig = 3)
            registra!("ESSMin", round(Int, min(ess(ch, :logk), ess(ch, :logmu))))
            lk = vec(Array(ch[:logk])); lm = vec(Array(ch[:logmu]))
            ppc, _ = ppc_minimo(s -> (k = exp(lk[s]); [Gamma(k, exp(lm[s]) / k) for _ in 1:length(x)]),
                                (ndraws = length(lk),), minimum(x))
            registra!("PPCMinEst", ppc; dig = 2)
            # preditiva posterior: P(próximo ano ≤ x₀) = E_post[p]
            registra!("PPreditiva", mean(ps); dig = 3)
            registra!("TPreditiva", 1 / mean(ps); sig = 2)
        else
            registra_ic!("TBayesVaga", rb.T[2], rb.T[1], rb.T[3]; sig = 2)
        end
    end
    mostra(res; titulo = "Posterior de T (evento ≤ 2012/13) — sensibilidade à priori")
    salva_ptab(res, "bayes_periodo_retorno")
    return (res = res, draws = draws)
end

# ============================================================================= P1: não estacionário
function oni_estacao_chuvosa()
    linhas = readlines(joinpath(RAIZ, "data", "raw", "indices", "oni.data"))
    tab = Dict{Tuple{Int, Int}, Float64}()
    for l in linhas[2:end]
        c = split(strip(l)); length(c) == 13 || continue
        a = tryparse(Int, c[1]); isnothing(a) && continue
        for m in 1:12
            v = parse(Float64, c[m+1]); v > -99 && (tab[(a, m)] = v)
        end
    end
    # média nov–abr do ano hidrológico a (nov/dez de a, jan–abr de a+1)
    Dict(a => mean([tab[(a, 11)], tab[(a, 12)], tab[(a + 1, 1)], tab[(a + 1, 2)], tab[(a + 1, 3)], tab[(a + 1, 4)]])
         for a in 1950:2025 if all(haskey(tab, k) for k in ((a, 11), (a, 12), (a + 1, 1), (a + 1, 2), (a + 1, 3), (a + 1, 4))))
end

function bloco_nao_estacionario(L, E, V)
    titulo("Modelos não estacionários")
    v = E.v; sort!(v, :ano_hidro)
    y = Float64.(v.prec); t = (v.ano_hidro .- mean(v.ano_hidro)) ./ 10
    tend = analise_tendencia(v.ano_hidro, v.prec; nome = "prec 82983")
    mostra(DataFrame([tend]); titulo = "Tendência na chuva anual (82983)")
    registra_ic!("SenChuva", tend.sen_decada, tend.sen_ic_inf, tend.sen_ic_sup; dig = 0)
    registra_p!("PMKChuva", tend.p_mk); registra_p!("PTFPWChuva", tend.p_mk_tfpw)
    registra_p!("PPettittChuvaLonga", tend.p_pettitt)
    nt = gamma_nao_estacionaria(y, t)
    oni = oni_estacao_chuvosa()
    ok = [haskey(oni, a) for a in v.ano_hidro]
    xo = [oni[a] for a in v.ano_hidro[ok]]
    ne = gamma_nao_estacionaria(y[ok], xo)
    tab = DataFrame([(modelo = "log μ ~ tempo (década)", n = length(y), b = nt.b, ep_b = nt.ep_b,
                      efeito_pct = 100(exp(nt.b) - 1), LR = nt.LR, p_LR = nt.p_LR, dAIC = nt.AIC - nt.AIC0),
                     (modelo = "log μ ~ ONI (nov–abr)", n = length(xo), b = ne.b, ep_b = ne.ep_b,
                      efeito_pct = 100(exp(ne.b) - 1), LR = ne.LR, p_LR = ne.p_LR, dAIC = ne.AIC - ne.AIC0)])
    mostra(tab; titulo = "Gamma não estacionária vs. estacionária (dAIC < 0 favorece a não estacionária)")
    salva_ptab(tab, "gamma_nao_estacionaria")
    registra!("NSTempoEfeito", tab.efeito_pct[1]; dig = 1); registra_p!("NSTempoP", tab.p_LR[1])
    registra!("NSTempoDAIC", tab.dAIC[1]; dig = 1)
    registra!("NSOniEfeito", tab.efeito_pct[2]; dig = 1); registra_p!("NSOniP", tab.p_LR[2])
    registra!("NSOniDAIC", tab.dAIC[2]; dig = 1)
    # p do evento de referência sob ONI = +1 (El Niño moderado), 0 e −1
    for (rot, o) in (("Nino", 1.0), ("Neutro", 0.0), ("Nina", -1.0))
        p = ne.p_em(E.x1213, o); registra!("TOni" * rot, 1 / p; sig = 2)
    end
    tc(a) = (a - mean(v.ano_hidro)) / 10
    registra!("TTempoDoisMilVinteCinco", 1 / nt.p_em(E.x1213, tc(2025)); sig = 2)

    # --- sensibilidade: inclinação fixada na taxa regional (referência das vizinhas)
    b_reg = missing
    if !isnothing(V)
        Q = V.Q
        nr = gamma_nao_estacionaria(Float64.(Q.referencia), tc.(Q.ano))
        b_reg = nr.b
        registra!("NSRegionalEfeito", 100(exp(b_reg) - 1); dig = 1)
    end
    nf = ismissing(b_reg) ? nothing : gamma_inclinacao_fixa(y, t, b_reg)

    # --- resultado principal: média por pesos de Akaike das candidatas não estacionárias,
    #     no mesmo esquema do caso estacionário (a Gamma sozinha fica como sensibilidade).
    #     Cada família entra com forma fixa (DP constante na normal, CV constante nas demais)
    #     e com forma log-linear no tempo, que libera a relação entre variabilidade e média.
    ult = maximum(v.ano_hidro); registra!("UltimoAnoNS", rotulo_hidro(ult))
    ma = ns_media_akaike(y, t)
    ma.tab.T_ult = [1 / f.F(E.x1213, tc(ult)) for f in ma.fits]
    ma.tab.T_2025 = [1 / f.F(E.x1213, tc(2025)) for f in ma.fits]
    sufixo(f) = f == "fixa" ? "" : "Var"
    for r in eachrow(ma.tab)
        registra!("TNS" * r.distribuicao * sufixo(r.forma), r.T_ult; sig = 2)
        registra!("PesoNS" * r.distribuicao * sufixo(r.forma), r.peso_akaike; dig = 2)
    end
    # razão de verossimilhança forma variável × fixa, por família (χ²₁)
    for D in CANDIDATAS
        nome = string(nameof(D))
        ffix = only(f for f in ma.fits if f.D === D && !f.forma_var)
        fvar = only(f for f in ma.fits if f.D === D && f.forma_var)
        registra_p!("PFormaVar" * nome, ccdf(Chisq(1), max(2(fvar.ll - ffix.ll), 0.0)))
        registra!("EfeitoForma" * nome, 100(exp(fvar.c) - 1); dig = 1)
    end
    fixas = ma.tab.forma .== "fixa"
    registra!("TNSModelosMin", minimum(ma.tab.T_ult[fixas]); sig = 2); registra!("TNSModelosMax", maximum(ma.tab.T_ult[fixas]); sig = 2)
    registra!("TNSVarMin", minimum(ma.tab.T_ult[.!fixas]); sig = 2); registra!("TNSVarMax", maximum(ma.tab.T_ult[.!fixas]); sig = 2)
    registra!("PesoNSVarTotal", sum(ma.tab.peso_akaike[.!fixas]); dig = 2)
    mf = ns_media_akaike(y, t; formas = (false,))          # sensibilidade: só forma fixa
    registra!("TNSMediaFixaUlt", 1 / mf.p_em(E.x1213, tc(ult)); sig = 2)
    sort!(ma.tab, :AIC)
    mostra(ma.tab; titulo = "Candidatas não estacionárias (posição log-linear no tempo)")
    salva_ptab(ma.tab, "nao_estacionarias_candidatas")
    registra!("MelhorAICNS", lowercase(ma.tab.distribuicao[1]))
    registra!("MelhorAICNSForma", ma.tab.forma[1])
    pbm = ns_media_akaike_boot(y, t, E.x1213, [tc(ult), tc(2025)])
    registra!("NBootNSMedia", size(pbm, 1))
    for (j, (rot, a)) in enumerate((("Ult", ult), ("DoisMilVinteCinco", 2025)))
        p0 = ma.p_em(E.x1213, tc(a))
        registra_ic!("TNSMedia" * rot, 1 / p0, 1 / quantile(pbm[:, j], 0.975), 1 / quantile(pbm[:, j], 0.025); sig = 2)
    end
    pu_m = ma.p_em(E.x1213, tc(ult))
    for nh in HORIZONTES
        registra_ic!("ProbNSMediaUltN" * string(nh), 100prob_ao_menos_um(pu_m, nh),
                     100quantile(prob_ao_menos_um.(pbm[:, 1], nh), 0.025),
                     100quantile(prob_ao_menos_um.(pbm[:, 1], nh), 0.975); dig = 0)
    end
    # inclinação regional, também como média entre as candidatas
    if !ismissing(b_reg)
        mr = ns_media_akaike(y, t; b = b_reg)
        pbr = ns_media_akaike_boot(y, t, E.x1213, [tc(2025)]; b = b_reg)
        p0 = mr.p_em(E.x1213, tc(2025))
        registra_ic!("TRegionalMedia", 1 / p0, 1 / quantile(pbr[:, 1], 0.975), 1 / quantile(pbr[:, 1], 0.025); sig = 2)
        for nh in HORIZONTES
            registra_ic!("ProbRegMediaN" * string(nh), 100prob_ao_menos_um(p0, nh),
                         100quantile(prob_ao_menos_um.(pbr[:, 1], nh), 0.025),
                         100quantile(prob_ao_menos_um.(pbr[:, 1], nh), 0.975); dig = 0)
        end
    end

    # --- bayesiano não estacionário
    po = posterior_ns(y, t)
    rh = maximum(diagnostico_cadeias(po.ch, s) for s in (:logk, :a, :b))
    registra!("RhatMaxNS", rh; dig = 3)
    registra_ic!("BayesNSEfeito", 100(exp(median(po.b)) - 1), 100(exp(quantile(po.b, 0.025)) - 1),
                 100(exp(quantile(po.b, 0.975)) - 1); dig = 1)
    registra!("ProbBNegativo", mean(po.b .< 0); dig = 3)
    pdraw(i, a) = cdf(Gamma(exp(po.logk[i]), exp(po.a[i] + po.b[i] * tc(a)) / exp(po.logk[i])), E.x1213)
    nd = length(po.b)
    p25 = [pdraw(i, 2025) for i in 1:nd]
    registra_ic!("TBayesNS", median(1 ./ p25), quantile(1 ./ p25, 0.025), quantile(1 ./ p25, 0.975); sig = 2)
    registra!("ESSMinNS", round(Int, minimum(ess(po.ch, s) for s in (:logk, :a, :b))))
    ppcns, _ = ppc_minimo(s -> (k = exp(po.logk[s]); [Gamma(k, exp(po.a[s] + po.b[s] * ti) / k) for ti in t]),
                          (ndraws = nd,), minimum(y))
    registra!("PPCMinNS", ppcns; dig = 2)
    # último ano hidrológico válido (sensibilidade Gamma bayesiana), sem extrapolar além dos dados
    pu =[pdraw(i, ult) for i in 1:nd]
    registra_ic!("TBayesNSUlt", median(1 ./ pu), quantile(1 ./ pu, 0.025), quantile(1 ./ pu, 0.975); sig = 2)
    registra!("TTempoUlt", 1 / nt.p_em(E.x1213, tc(ult)); sig = 2)
    for nh in HORIZONTES
        qq = quantile([1 - (1 - pu[i])^nh for i in 1:nd], [0.025, 0.5, 0.975])
        registra_ic!("ProbNSUltN" * string(nh), 100qq[2], 100qq[1], 100qq[3]; dig = 0)
    end

    # sensibilidade: critério FCA ≤ 7,5 % (inclui anos recentes quase completos, p.ex. 2021/22 e 2023/24)
    d75 = sort(L.dh[coalesce.(L.dh.fca .<= 0.075, false) .& (L.dh.prec_obs .> 0), :], :ano_hidro)
    novos = setdiff(d75.ano_hidro, v.ano_hidro)
    registra!("AnosFCASete", join(rotulo_hidro.(novos), ", ")); registra!("NanosFCASete", nrow(d75))
    y75 = Float64.(d75.prec_obs); t75 = tc.(d75.ano_hidro)
    n75 = gamma_nao_estacionaria(y75, t75)
    tr75 = analise_tendencia(d75.ano_hidro, y75; nome = "prec FCA ≤ 7,5 %")
    registra_ic!("SenChuvaFCASete", tr75.sen_decada, tr75.sen_ic_inf, tr75.sen_ic_sup; dig = 0)
    registra!("NSFCASeteEfeito", 100(exp(n75.b) - 1); dig = 1)
    registra!("TFCASeteUlt", 1 / n75.p_em(E.x1213, tc(maximum(d75.ano_hidro))); sig = 2)
    registra!("UltimoAnoFCASete", rotulo_hidro(maximum(d75.ano_hidro)))

    # --- probabilidade de ≥ 1 evento em n anos, por cenário
    cen = NamedTuple[]
    g = fit_mle(Gamma, y); p_est = cdf(g, E.x1213)
    for nh in HORIZONTES
        congel = [1 - (1 - p25[i])^nh for i in 1:nd]
        continua = [prob_n_variavel(j -> pdraw(i, 2025 + j), nh) for i in 1:nd]
        q(vv) = quantile(vv, [0.025, 0.5, 0.975])
        qc = q(congel); qt = q(continua)
        push!(cen, (n_anos = nh, estacionario_MV = prob_ao_menos_um(p_est, nh),
                    ns_MV_congelado = prob_ao_menos_um(nt.p_em(E.x1213, tc(2025)), nh),
                    ns_MV_continua = prob_n_variavel(j -> nt.p_em(E.x1213, tc(2025 + j)), nh),
                    regional_congelado = isnothing(nf) ? missing : prob_ao_menos_um(nf.p_em(E.x1213, tc(2025)), nh),
                    regional_continua = isnothing(nf) ? missing : prob_n_variavel(j -> nf.p_em(E.x1213, tc(2025 + j)), nh),
                    bayes_congelado = qc[2], bayes_congelado_inf = qc[1], bayes_congelado_sup = qc[3],
                    bayes_continua = qt[2], bayes_continua_inf = qt[1], bayes_continua_sup = qt[3]))
        s = string(nh)
        registra_ic!("ProbNSCongN" * s, 100qc[2], 100qc[1], 100qc[3]; dig = 0)
        registra_ic!("ProbNSContN" * s, 100qt[2], 100qt[1], 100qt[3]; dig = 0)
        registra!("ProbEstN" * s, 100prob_ao_menos_um(p_est, nh); dig = 0)
        isnothing(nf) || registra!("ProbRegN" * s, 100prob_ao_menos_um(nf.p_em(E.x1213, tc(2025)), nh); dig = 0)
    end
    tcen = DataFrame(cen)
    mostra(tcen; titulo = "P(≥ 1 ano com chuva ≤ 2012/13 em n anos, a partir de 2025) por cenário")
    salva_ptab(tcen, "probabilidade_cenarios")
    tT = DataFrame(modelo = ["estacionário (Gamma, MV)", "não estacionário local (MV)", "não estacionário, inclinação regional (MV)",
                             "não estacionário local (bayes, mediana)"],
                   T_2025 = [1 / p_est, 1 / nt.p_em(E.x1213, tc(2025)),
                             isnothing(nf) ? missing : 1 / nf.p_em(E.x1213, tc(2025)), median(1 ./ p25)],
                   efeito_pct_decada = [0.0, 100(exp(nt.b) - 1), ismissing(b_reg) ? missing : 100(exp(b_reg) - 1),
                                        100(exp(median(po.b)) - 1)])
    mostra(tT; titulo = "Período de retorno em 2025 do evento ≤ 2012/13")
    salva_ptab(tT, "periodo_retorno_2025")
    isnothing(nf) || registra!("TRegionalDoisMilVinteCinco", tT.T_2025[3]; sig = 2)
    # T do evento ao longo do tempo (para figura): MV local e regional
    anos = 1965:2025
    salva_ptab(DataFrame(ano = collect(anos), T_local = [1 / nt.p_em(E.x1213, tc(a)) for a in anos],
                         T_regional = [isnothing(nf) ? missing : 1 / nf.p_em(E.x1213, tc(a)) for a in anos],
                         T_bayes_med = [median(1 ./ [pdraw(i, a) for i in 1:20:nd]) for a in anos]), "periodo_retorno_no_tempo")
    return (tend = tend, tab = tab, oni = oni, nt = nt, nf = nf, po = po, cen = tcen, tT = tT, ma = ma)
end

# ============================================================================= SPI / SPEI / corridas
function bloco_indices(L, E)
    titulo("SPI, SPEI e teoria de corridas")
    dm = sort(L.dm, :data)
    # garante grade mensal contínua
    grade = DataFrame(data = minimum(dm.data):Month(1):maximum(dm.data))
    dm = leftjoin(grade, dm, on = :data); sort!(dm, :data)
    dm.ano = year.(dm.data); dm.mes = month.(dm.data)
    for k in (3, 6, 12)
        dm[!, Symbol(:spi, k)], _ = spi(dm.prec, dm.mes, dm.ano; k = k, ref = REF)
    end
    dm.eto = [eto_mensal(a, m, tx, tn, LAT_82983) for (a, m, tx, tn) in zip(dm.ano, dm.mes, dm.tmax, dm.tmin)]
    dm.balanco = [ismissing(p) || ismissing(e) ? missing : p - e for (p, e) in zip(dm.prec, dm.eto)]
    dm.spei12, ajll = spei(dm.balanco, dm.mes, dm.ano; k = 12, ref = REF, dist = :loglogistica)
    dm.spei12_gev, ajgev = spei(dm.balanco, dm.mes, dm.ano; k = 12, ref = REF, dist = :gev)

    # Stagge et al. (2015): qualidade do ajuste pela normalidade do índice resultante
    sw = DataFrame(map([(:spi12, "SPI-12 (Gamma)"), (:spei12, "SPEI-12 (log-logística)"), (:spei12_gev, "SPEI-12 (GEV)")]) do (v, nm)
        ref = [z for (z, a) in zip(dm[!, v], dm.ano) if !ismissing(z) && REF[1] <= a <= REF[2]]
        t = ShapiroWilkTest(Float64.(ref[1:12:end]))      # 1 valor por ano (evita sobreposição)
        (indice = nm, n_meses = count(!ismissing, dm[!, v]), media_ref = mean(ref), dp_ref = std(ref),
         W = t.W, p_shapiro = pvalue(t))
    end)
    mostra(sw; titulo = "Normalidade dos índices no período de referência")
    salva_ptab(sw, "indices_normalidade")

    # índice no fim do ano hidrológico (junho) = SPI-12/SPEI-12 do ano hidrológico
    jun = dm[dm.mes .== 6, [:ano, :spi12, :spei12, :spei12_gev]]
    jun.ano_hidro = jun.ano .- 1
    for (rot, a) in (("DozeTreze", 2012), ("DezesseisDezessete", 2016), ("OnzeDoze", 2011))
        r = jun[jun.ano_hidro .== a, :]
        nrow(r) == 1 || continue
        registra!("SPI" * rot, r.spi12[1]; dig = 2); registra!("SPEI" * rot, r.spei12[1]; dig = 2)
    end
    salva_ptab(select(dm, :data, :ano, :mes, :prec, :tmax, :tmin, :eto, :balanco, :spi3, :spi6, :spi12, :spei12, :spei12_gev),
               "indices_mensais")

    # corridas no SPI-12 e SPEI-12
    ev = Dict{Symbol, DataFrame}()
    for v in (:spi12, :spei12)
        e = corridas(dm.data, dm[!, v]; z0 = -1.0)
        mostra(e; titulo = "Eventos de seca ($(v) < −1)")
        salva_ptab(e, "eventos_" * String(v)); ev[v] = e
    end
    e = ev[:spi12]
    iref = findfirst(r -> r.inicio <= Date(2013, 6, 1) <= r.fim, eachrow(e))
    anos_cobertos = count(!ismissing, dm.spi12) / 12
    registra!("NeventosSPI", nrow(e)); registra!("AnosCobertosSPI", anos_cobertos; dig = 0)
    if !isnothing(iref)
        r = e[iref, :]
        registra!("EvRefInicio", Dates.format(r.inicio, "mm/yyyy")); registra!("EvRefFim", Dates.format(r.fim, "mm/yyyy"))
        registra!("EvRefDuracao", r.duracao); registra!("EvRefSeveridade", r.severidade; dig = 1)
        registra!("EvRefIntensidade", r.intensidade; dig = 2); registra!("EvRefPico", r.pico; dig = 2)
        rank_sev = count(>=(r.severidade), e.severidade)
        registra!("EvRefRankSeveridade", rank_sev)
        λ = nrow(e) / anos_cobertos
        # severidade ~ Exponencial (ajuste simples; n pequeno) → T = 1/(λ·P(S ≥ s_ref))
        es = fit_mle(Exponential, e.severidade)
        Tsev = 1 / (λ * ccdf(es, r.severidade))
        registra!("TaxaEventos", λ * 10; dig = 1)    # eventos por década
        registra!("TSeveridadeRef", Tsev; sig = 2)
        registra!("TSeveridadeEmp", anos_cobertos / rank_sev; sig = 2)
        @printf("Evento de referência (SPI-12): %s a %s, D = %d meses, S = %.1f; %d evento(s) com S ≥ S_ref em %.0f anos; T(exp) ≈ %.0f anos\n",
                r.inicio, r.fim, r.duracao, r.severidade, rank_sev, anos_cobertos, Tsev)
    end
    return (dm = dm, jun = jun, ev = ev, sw = sw)
end

# ============================================================================= Markov
function bloco_markov(L, I)
    titulo("Persistência: cadeia de Markov")
    v = L.dh[L.dh.valido, :]
    corte = quantile(Float64.(v.prec), 1 / 3)
    est = [ismissing(p) ? missing : (p <= corte ? "S" : "N") for p in L.dh.prec]
    C = transicoes(L.dh.ano_hidro, est)
    m = markov_posterior(C)
    @printf("Corte (tercil inferior): %.0f mm | contagens S→S %d, S→N %d, N→S %d, N→N %d\n", corte, C[1, 1], C[1, 2], C[2, 1], C[2, 2])
    @printf("P(S|S) = %.2f [%.2f; %.2f] | P(S|N) = %.2f [%.2f; %.2f] | P(P(S|S) > P(S|N)) = %.2f\n",
            m.pSS[2], m.pSS[1], m.pSS[3], m.pNS[2], m.pNS[1], m.pNS[3], m.prob_dif_pos)
    ft = FisherExactTest(C[1, 1], C[1, 2], C[2, 1], C[2, 2])
    registra!("CorteSeco", corte; dig = 0); registra!("NtransicoesMarkov", sum(C))
    registra!("NSS", C[1, 1]); registra!("NSN", C[1, 2]); registra!("NNS", C[2, 1]); registra!("NNN", C[2, 2])
    registra_ic!("PSS", m.pSS[2], m.pSS[1], m.pSS[3]; dig = 2)
    registra_ic!("PNS", m.pNS[2], m.pNS[1], m.pNS[3]; dig = 2)
    registra_ic!("DifPersistencia", m.dif[2], m.dif[1], m.dif[3]; dig = 2)
    registra!("ProbPersistencia", m.prob_dif_pos; dig = 2); registra_p!("PFisher", pvalue(ft))
    registra_ic!("DuracaoSeca", m.duracao_seca[2], m.duracao_seca[1], m.duracao_seca[3]; dig = 1)
    registra_ic!("Recorrencia", m.recorrencia[2], m.recorrencia[1], m.recorrencia[3]; dig = 1)
    # sensibilidade: seco = SPI anual ≤ −1
    z, _ = spi_anual(L.dh.prec)
    est2 = [ismissing(q) ? missing : (q <= -1 ? "S" : "N") for q in z]
    C2 = transicoes(L.dh.ano_hidro, est2); m2 = markov_posterior(C2)
    registra_ic!("PSSspi", m2.pSS[2], m2.pSS[1], m2.pSS[3]; dig = 2)
    registra_ic!("PNSspi", m2.pNS[2], m2.pNS[1], m2.pNS[3]; dig = 2)
    # sem tendência: estado pela PIT uᵢ = Fᵢ(yᵢ) sob a Gamma não estacionária (média
    # log-linear no tempo). Seco ⇔ uᵢ ≤ 1/3, isto é, tercil inferior *do seu próprio ano*.
    # Separa persistência real de anos secos agrupados só porque a média caiu.
    tc(a) = (a - mean(v.ano_hidro)) / 10
    nt = gamma_nao_estacionaria(Float64.(v.prec), tc.(v.ano_hidro))
    est3 = [ismissing(p) ? missing : (nt.p_em(p, tc(a)) <= 1 / 3 ? "S" : "N") for (p, a) in zip(L.dh.prec, L.dh.ano_hidro)]
    C3 = transicoes(L.dh.ano_hidro, est3); m3 = markov_posterior(C3)
    registra_ic!("PSSdet", m3.pSS[2], m3.pSS[1], m3.pSS[3]; dig = 2)
    registra_ic!("PNSdet", m3.pNS[2], m3.pNS[1], m3.pNS[3]; dig = 2)
    registra!("ProbPersistenciaDet", m3.prob_dif_pos; dig = 2)
    Cs = (C, C2, C3); ms = (m, m2, m3)
    tab = DataFrame(definicao = ["tercil inferior", "SPI ≤ −1", "tercil sem tendência (PIT)"],
                    SS = [c[1, 1] for c in Cs], SN = [c[1, 2] for c in Cs], NS = [c[2, 1] for c in Cs], NN = [c[2, 2] for c in Cs],
                    pSS = [x.pSS[2] for x in ms], pSS_inf = [x.pSS[1] for x in ms], pSS_sup = [x.pSS[3] for x in ms],
                    pNS = [x.pNS[2] for x in ms], pNS_inf = [x.pNS[1] for x in ms], pNS_sup = [x.pNS[3] for x in ms],
                    prob_persist = [x.prob_dif_pos for x in ms])
    mostra(tab; titulo = "Matriz de transição (posterior Dirichlet, priori de Jeffreys)")
    salva_ptab(tab, "markov")
    return (C = C, m = m, tab = tab, corte = corte)
end

# ============================================================================= sensibilidade (série combinada)
function bloco_sensibilidade(L, E, R)
    titulo("Sensibilidade: série combinada A307 + preenchimento (2003–2024)")
    c = CSV.read(joinpath(DIR_PROC, "anual_hidrologico.csv"), DataFrame)
    c = dropmissing(c[c.valido .== true, [:ano_hidro, :rotulo, :prec]])
    l = L.dh[L.dh.valido .& (L.dh.ano_hidro .>= 2003), [:ano_hidro, :prec]]
    j = innerjoin(rename(c, :prec => :prec_comb), rename(l, :prec => :prec_82983), on = :ano_hidro)
    mostra(j; titulo = "Totais anuais nos anos comuns")
    salva_ptab(j, "sensibilidade_series")
    # p do evento (≤ 2012/13 da 82983) usando a Gamma da série longa, mas o x₀ da série combinada
    x = Float64.(E.v.prec)
    x0c = only(c.prec[c.ano_hidro .== 2012])
    pr = perfil_p_gamma(x, x0c)
    registra!("ChuvaDozeTrezeComb", x0c; dig = 1)
    registra_ic!("TComb", 1 / pr.p, 1 / pr.sup, 1 / pr.inf; sig = 2)
    # p com a Gamma ajustada só aos anos da era da A307 (2003 em diante), na própria 82983
    cc = L.dh[L.dh.valido .& (L.dh.ano_hidro .>= 2003), :]
    xc = Float64.(cc.prec); pc = perfil_p_gamma(xc, E.x1213)
    registra!("NanosCurto", length(xc)); registra!("PeriodoCurto", rotulo_hidro(minimum(cc.ano_hidro)) * "–" * rotulo_hidro(maximum(cc.ano_hidro)))
    registra_ic!("TCurto", 1 / pc.p, 1 / pc.sup, 1 / pc.inf; sig = 2)
    @printf("x₀ combinado = %.1f → T = %.0f [%.0f; %.0f] | série curta (n = %d): T = %.0f [%.0f; %.0f]\n",
            x0c, 1 / pr.p, 1 / pr.sup, 1 / pr.inf, length(xc), 1 / pc.p, 1 / pc.sup, 1 / pc.inf)
    return j
end

# ============================================================================= P2: sazonalidade
"Série diária do ano hidrológico `a` (jul–jun), com `missing` nos dias sem registro."
function diaria_hidro(d, a)
    pos = Dict(d.data .=> d.prec)
    [get(pos, x, missing) for x in Date(a, MES_INICIO_HIDRO, 1):Day(1):Date(a + 1, MES_INICIO_HIDRO, 1) - Day(1)]
end

function bloco_p2_sazonalidade(L)
    titulo("P2 — Sazonalidade (82983)")
    v = L.dh[L.dh.valido, :]
    Rbar = L.clim_anual / 365.25
    linhas = map(eachrow(v)) do r
        p = diaria_hidro(L.d, r.ano_hidro)
        p0 = [ismissing(x) ? 0.0 : x for x in p]
        mes = [month(x) for x in Date(r.ano_hidro, MES_INICIO_HIDRO, 1):Day(1):Date(r.ano_hidro + 1, MES_INICIO_HIDRO, 1) - Day(1)]
        x12 = [sum(p0[mes .== m]) for m in 1:12]
        c = circular_chuva(p0); lb = liebmann(p, Rbar)
        # estação indefinida: a anomalia acumulada não atinge um máximo dentro do ano
        # (fim nos 3 últimos dias) ou a "estação" dura < 30 dias — típico de anos muito secos
        indef = lb.fim >= length(p) - 2 || lb.duracao < 30
        (ano_hidro = r.ano_hidro, rotulo = r.rotulo, prec = r.prec, SI = indice_sazonalidade(x12),
         centroide = c.centroide, concentracao = c.concentracao, estacao_definida = !indef,
         inicio = indef ? missing : lb.inicio, fim = indef ? missing : lb.fim, duracao = indef ? missing : lb.duracao,
         frac_nov_abr = sum(x12[[11, 12, 1, 2, 3, 4]]) / sum(x12))
    end
    s = DataFrame(linhas)
    salva_ptab(s, "sazonalidade_anual")
    dia(n) = Dates.format(Date(2001, 7, 1) + Day(round(Int, n) - 1), "dd/mm")
    tend = DataFrame([tendencia_completa(s.ano_hidro, s[!, c]; nome = String(c))
                      for c in (:SI, :centroide, :concentracao, :inicio, :fim, :duracao, :frac_nov_abr)])
    mostra(tend; titulo = "Tendência dos indicadores de sazonalidade (por década)")
    salva_ptab(tend, "sazonalidade_tendencias")
    registra!("SIMedio", mean(s.SI); dig = 2); registra!("ConcMedia", mean(s.concentracao); dig = 2)
    registra!("CentroideMedio", dia(mean(s.centroide)))
    registra!("InicioMediano", dia(median(skipmissing(s.inicio)))); registra!("FimMediano", dia(median(skipmissing(s.fim))))
    registra!("DuracaoMediana", round(Int, median(skipmissing(s.duracao))))
    registra!("NanosSemEstacao", count(.!s.estacao_definida))
    s.sem_estacao = Float64.(.!s.estacao_definida); s.decada = (s.ano_hidro .- mean(s.ano_hidro)) ./ 10
    lg = glm(@formula(sem_estacao ~ decada), s, Binomial(), LogitLink())
    lg0 = glm(@formula(sem_estacao ~ 1), s, Binomial(), LogitLink())
    LRl = 2(loglikelihood(lg) - loglikelihood(lg0)); bl = coef(lg)[2]; el = stderror(lg)[2]
    registra_ic!("ORSemEstacao", exp(bl), exp(bl - 1.96el), exp(bl + 1.96el); dig = 2)
    registra_p!("PSemEstacao", ccdf(Chisq(1), LRl))
    @printf("Anos sem estação definida: OR por década = %.2f [%.2f; %.2f], LR p = %.3f\n", exp(bl), exp(bl - 1.96el), exp(bl + 1.96el), ccdf(Chisq(1), LRl))
    registra!("AnosSemEstacao", join(s.rotulo[.!s.estacao_definida], ", "))
    registra!("FracNovAbr", 100mean(s.frac_nov_abr); dig = 0)
    for (i, nm) in enumerate(("SI", "Centroide", "Conc", "Inicio", "Fim", "Duracao", "FracNovAbr"))
        r = tend[i, :]; dg = nm in ("SI", "Conc", "FracNovAbr") ? 3 : 1
        registra_ic!("Sen" * nm, r.sen_decada, r.sen_ic_inf, r.sen_ic_sup; dig = dg)
        registra_p!("PMK" * nm, r.p_mk); registra_p!("PTFPW" * nm, r.p_mk_tfpw); registra_p!("PHR" * nm, r.p_hamed_rao)
        registra!("MDT" * nm, r.mdt_decada; dig = dg)
    end

    # climatologia mensal: 1961–1990 × 1991–2020 (meses válidos)
    dm = L.dm
    clim = combine(groupby(dm, :mes; sort = true)) do g
        a = collect(skipmissing(g.prec[1961 .<= g.ano .<= 1990])); b = collect(skipmissing(g.prec[1991 .<= g.ano .<= 2020]))
        ta = collect(skipmissing(g.tmax[1991 .<= g.ano .<= 2020])); tn = collect(skipmissing(g.tmin[1991 .<= g.ano .<= 2020]))
        (prec_6190 = isempty(a) ? missing : mean(a), n_6190 = length(a), prec_9120 = mean(b), n_9120 = length(b),
         p10_9120 = quantile(b, 0.1), p90_9120 = quantile(b, 0.9), tmax_9120 = mean(ta), tmin_9120 = mean(tn))
    end
    mostra(clim; titulo = "Climatologia mensal da 82983 (meses válidos)")
    salva_ptab(clim, "climatologia_mensal")
    registra!("MesMaisChuvoso", ("jan", "fev", "mar", "abr", "mai", "jun", "jul", "ago", "set", "out", "nov", "dez")[argmax(clim.prec_9120)])
    registra!("ChuvaMarNormal", clim.prec_9120[3]; dig = 0)
    return (s = s, tend = tend, clim = clim)
end

# ============================================================================= P2: ETCCDI e médias anuais
function bloco_p2_indices(L)
    titulo("P2 — Índices ETCCDI e médias anuais (82983)")
    d = L.d
    umidos = [x for (x, a) in zip(d.prec, d.ano) if !ismissing(x) && x >= 1 && REF[1] <= a <= REF[2]]
    p95 = quantile(umidos, 0.95); registra!("PNoventaCincoUmido", p95; dig = 1)
    v = L.dh[L.dh.valido, :]
    ich = DataFrame(map(eachrow(v)) do r
        merge((ano_hidro = r.ano_hidro,), indices_chuva(diaria_hidro(d, r.ano_hidro), p95))
    end)
    # CDD por ano civil (a estação seca, mai–out, fica inteira dentro do ano civil); ≤ 5 dias ausentes
    cdd = combine(groupby(d, :ano; sort = true)) do g
        nd = daysinyear(g.ano[1]); aus = count(ismissing, g.prec) + nd - nrow(g)
        (CDD_civil = aus <= 5 ? maior_corrida(g.prec, <(1.0)) : missing, ausentes = aus)
    end
    # temperatura e outras médias: ano civil, ≥ 10 meses válidos, via anomalias mensais
    dm = L.dm
    anu = DataFrame(ano = sort(unique(dm.ano)))
    for vv in (:tmax, :tmin, :tmed, :ur, :vvel)
        cl = Dict(m => mean(skipmissing(dm[(dm.mes .== m) .& (REF[1] .<= dm.ano .<= REF[2]), vv])) for m in 1:12)
        anu[!, vv] = map(anu.ano) do a
            g = dm[dm.ano .== a, :]
            an = [x - cl[m] for (x, m) in zip(g[!, vv], g.mes) if !ismissing(x)]
            length(an) >= MIN_MESES_ANO ? mean(values(cl)) + mean(an) : missing
        end
    end
    anu.dtr = anu.tmax .- anu.tmin
    # extremos de temperatura (percentis por dia do calendário, base 1991–2020)
    p90 = percentil_calendario(d.data, d.tmax, 0.9); p10 = percentil_calendario(d.data, d.tmin, 0.1)
    ext = combine(groupby(d, :ano; sort = true)) do g
        a = g.ano[1]; okx = !ismissing(anu.tmax[anu.ano .== a][1]); okn = !ismissing(anu.tmin[anu.ano .== a][1])
        sx = [(dt, x) for (dt, x) in zip(g.data, g.tmax) if !ismissing(x)]
        sn = [(dt, x) for (dt, x) in zip(g.data, g.tmin) if !ismissing(x)]
        (TXx = okx ? maximum(last.(sx)) : missing, TNn = okn ? minimum(last.(sn)) : missing,
         TX90p = okx ? 100count(((dt, x),) -> x > p90[dayofyear(dt)], sx) / length(sx) : missing,
         TN10p = okn ? 100count(((dt, x),) -> x < p10[dayofyear(dt)], sn) / length(sn) : missing)
    end
    anu = leftjoin(leftjoin(anu, ext, on = :ano), cdd, on = :ano); sort!(anu, :ano)
    salva_ptab(ich, "etccdi_chuva_hidro"); salva_ptab(anu, "anual_civil_82983")

    series = [(ich.ano_hidro, ich.R1mm, "R1mm (dias)"), (ich.ano_hidro, ich.Rx1day, "Rx1day (mm)"),
              (ich.ano_hidro, ich.SDII, "SDII (mm/dia)"), (ich.ano_hidro, ich.R95p, "R95p (mm)"),
              (ich.ano_hidro, ich.CWD, "CWD (dias)"), (anu.ano, anu.CDD_civil, "CDD (dias)"),
              (anu.ano, anu.tmax, "Tmáx (°C)"), (anu.ano, anu.tmin, "Tmín (°C)"), (anu.ano, anu.tmed, "Tméd (°C)"),
              (anu.ano, anu.dtr, "DTR (°C)"), (anu.ano, anu.ur, "UR (%)"), (anu.ano, anu.vvel, "Vento (m/s)"),
              (anu.ano, anu.TXx, "TXx (°C)"), (anu.ano, anu.TNn, "TNn (°C)"),
              (anu.ano, anu.TX90p, "TX90p (% dias)"), (anu.ano, anu.TN10p, "TN10p (% dias)")]
    tend = DataFrame([tendencia_completa(t, y; nome = n) for (t, y, n) in series])
    tend.q_tfpw = bh_ajuste(tend.p_mk_tfpw); tend.q_hr = bh_ajuste(tend.p_hamed_rao)
    mostra(tend; titulo = "Tendências P2 (Sen por década; q de BH sobre TFPW e Hamed–Rao)")
    salva_ptab(tend, "p2_tendencias")
    chaves = ("Rummm", "Rxumday", "SDII", "Rnovecincop", "CWD", "CDD", "Tmax", "Tmin", "Tmed", "DTR", "UR", "Vento",
              "TXx", "TNn", "TXnovezerop", "TNumzerop")
    casas = Dict("Rnovecincop" => 0, "Rummm" => 1, "Rxumday" => 1, "CWD" => 1, "CDD" => 1, "UR" => 1,
                 "TXnovezerop" => 1, "TNumzerop" => 1)
    for (k, r) in zip(chaves, eachrow(tend))
        k = nome_macro(k); dg = get(casas, k, 2)
        registra_ic!("Sen" * k, r.sen_decada, r.sen_ic_inf, r.sen_ic_sup; dig = dg)
        registra_p!("QTFPW" * k, r.q_tfpw); registra_p!("QHR" * k, r.q_hr)
        registra!("MDT" * k, r.mdt_decada; dig = dg); registra!("N" * k, r.n)
        registra!("Periodo" * k, r.periodo)
    end
    registra!("NtendSig", count(<(0.05), tend.q_tfpw)); registra!("NtendTotal", nrow(tend))
    # sensibilidade à quebra A307 × 82983 (≈ 2009): tendência só até 2008
    for (vv, k) in ((:tmin, "Tmin"), (:tmax, "Tmax"))
        sub = anu[anu.ano .<= 2008, :]
        tt = tendencia_completa(sub.ano, sub[!, vv]; nome = k * " até 2008")
        registra_ic!("Sen" * k * "AteDoisMilOito", tt.sen_decada, tt.sen_ic_inf, tt.sen_ic_sup; dig = 2)
        registra_p!("PTFPW" * k * "AteDoisMilOito", tt.p_mk_tfpw)
    end
    return (ich = ich, anu = anu, tend = tend)
end

# ============================================================================= P2: homogeneidade da temperatura
"Diferença anual de Tmáx/Tmín 82983 − vizinha (Cabrobó, Remanso): SNHT + Pettitt; sensibilidade da Tmáx às quebras."
function bloco_p2_homog_temp(L, X)
    titulo("P2 — Homogeneidade da temperatura da 82983 (diferenças com vizinhas)")
    dir = joinpath(RAIZ, "data", "raw", "bdmep")
    # segmento inicial (antes da lacuna de 1966–70), lido do arquivo bruto, sem o corte de INICIO_TEMP
    bruto = le_bdmep_diario(ARQ_LONGO)
    for (c, nm) in (("82886", "Cabrobo"), ("82979", "Remanso"))
        f = filter(x -> startswith(x, "dados_$(c)_D_"), readdir(dir)); isempty(f) && continue
        v = le_bdmep_diario(joinpath(dir, f[1]))
        j = innerjoin(select(bruto, :data, :tmax, :tmin), select(v, :data, :tmax => :tx, :tmin => :tn), on = :data)
        for (a, b, rot) in ((:tmax, :tx, "Tmax"), (:tmin, :tn, "Tmin"))
            dd = [(year(dt), x - y) for (dt, x, y) in zip(j.data, j[!, a], j[!, b]) if !ismissing(x) && !ismissing(y)]
            ini = [d for (an, d) in dd if an <= 1965]; pos = [d for (an, d) in dd if an >= INICIO_TEMP]
            registra!("DifIni" * rot * nm, mean(ini); dig = 1); registra!("DifPos" * rot * nm, mean(pos); dig = 1)
            @printf("%s × %s: diferença média 1963–65 = %.2f °C; ≥ %d = %.2f °C\n", rot, nm, mean(ini), INICIO_TEMP, mean(pos))
        end
    end
    linhas = NamedTuple[]
    for (c, nm) in (("82886", "Cabrobó"), ("82979", "Remanso"))
        f = filter(x -> startswith(x, "dados_$(c)_D_"), readdir(dir)); isempty(f) && continue
        v = le_bdmep_diario(joinpath(dir, f[1]))
        j = innerjoin(select(L.d, :data, :tmax, :tmin), select(v, :data, :tmax => :tmax_v, :tmin => :tmin_v), on = :data)
        j.ano = year.(j.data)
        for (a, b, rot) in ((:tmax, :tmax_v, "Tmáx"), (:tmin, :tmin_v, "Tmín"))
            g = combine(groupby(dropmissing(j[:, [:ano, a, b]]), :ano), [a, b] => ((x, y) -> length(x) >= 300 ? mean(x .- y) : missing) => :dif)
            g = dropmissing(g); nrow(g) < 10 && continue
            sn = snht_teste(Float64.(g.ano), g.dif); pt = pettitt(Float64.(g.ano), g.dif)
            push!(linhas, (variavel = rot, vizinha = nm, n = nrow(g), inicio = minimum(g.ano), fim = maximum(g.ano),
                           p_snht = sn.p, quebra_snht = Int(sn.ultimo_antes), salto = sn.media_depois - sn.media_antes,
                           p_pettitt = pt.p))
        end
    end
    h = DataFrame(linhas); h.q = bh_ajuste(h.p_snht)
    mostra(h; titulo = "SNHT nas diferenças anuais de temperatura")
    salva_ptab(h, "homogeneidade_temperatura")
    registra!("NquebrasTemp", count(<(0.05), h.q)); registra!("NtestesTemp", nrow(h))
    for r in eachrow(h)
        k = (r.variavel == "Tmáx" ? "Tmax" : "Tmin") * (r.vizinha == "Cabrobó" ? "Cabrobo" : "Remanso")
        registra!("Quebra" * k, string(r.quebra_snht + 1))            # primeiro ano depois da quebra (texto: sem separador de milhar)
        registra!("Salto" * k, r.salto; dig = 2); registra_p!("QQuebra" * k, r.q)
    end

    # --- sensibilidade da tendência da Tmáx às quebras frente às vizinhas. A comparação não
    #     diz qual estação mudou, então se usam três hipóteses complementares:
    #     (a) a 82983 mudou: subtrai-se o salto da diferença com Cabrobó a partir da quebra;
    #     (b) agnóstica: MQO com tendência + degrau na data da quebra (absorve qualquer salto,
    #         inclusive climático, e por isso é conservadora);
    #     (c) começo da série: a Tmáx de 1973–74 difere das duas vizinhas; usa-se 1975 em diante.
    anu = X.anu; s = dropmissing(anu[:, [:ano, :tmax]])
    hc = h[(h.variavel .== "Tmáx") .& (h.vizinha .== "Cabrobó"), :]
    if nrow(hc) == 1
        a0 = hc.quebra_snht[1] + 1; sal = hc.salto[1]
        aj = s.tmax .- sal .* (s.ano .>= a0)
        ta = tendencia_completa(s.ano, aj; nome = "Tmáx ajustada")
        registra_ic!("SenTmaxAjustada", ta.sen_decada, ta.sen_ic_inf, ta.sen_ic_sup; dig = 2)
        registra_p!("PTFPWTmaxAjustada", ta.p_mk_tfpw)
        dd = DataFrame(y = s.tmax, t = (s.ano .- mean(s.ano)) ./ 10, d = Float64.(s.ano .>= a0))
        m = lm(@formula(y ~ t + d), dd); ic = confint(m)
        registra_ic!("BTmaxDegrau", coef(m)[2], ic[2, 1], ic[2, 2]; dig = 2)
        registra!("DegrauTmax", coef(m)[3]; dig = 2)
        @printf("Tmáx: Sen ajustada = %.2f [%.2f; %.2f] | MQO com degrau em %d: b = %.2f [%.2f; %.2f], degrau = %.2f\n",
                ta.sen_decada, ta.sen_ic_inf, ta.sen_ic_sup, a0, coef(m)[2], ic[2, 1], ic[2, 2], coef(m)[3])
    end
    for (vv, k) in ((:tmax, "Tmax"), (:tmin, "Tmin"))
        ss = dropmissing(anu[anu.ano .>= 1975, [:ano, vv]])
        tt = tendencia_completa(ss.ano, ss[!, vv]; nome = k * " desde 1975")
        registra_ic!("Sen" * k * "DesdeSetentaCinco", tt.sen_decada, tt.sen_ic_inf, tt.sen_ic_sup; dig = 2)
    end
    return h
end

# ============================================================================= P2: quantílica e STL
function bloco_p2_quantilica(L)
    titulo("P2 — Regressão quantílica da temperatura diária (82983)")
    q = vcat([insertcols!(quantilica_tempo(L.d.data, L.d[!, v]), 1, :variavel => String(v)) for v in (:tmax, :tmin)]...)
    mostra(q; titulo = "β₁ (°C/década) por quantil, IC por bootstrap em blocos anuais")
    salva_ptab(q, "quantilica_temperatura")
    for r in eachrow(q)
        k = (r.variavel == "tmax" ? "QRTmax" : "QRTmin") * ("Q" * string(round(Int, 100r.tau)))
        registra_ic!(k, r.beta_decada, r.ic_inf, r.ic_sup; dig = 2)
    end
    return q
end

function bloco_p2_stl(L)
    titulo("P2 — Decomposição STL (mensal)")
    dm = sort(L.dm, :data)
    out = DataFrame(data = dm.data)
    info = NamedTuple[]
    for v in (:tmax, :tmin, :prec)
        # período contínuo: do 1º ao último mês com dado, e só a partir de 1971 (lacuna de 1965–70)
        sel = dm.ano .>= 1971
        x = dm[sel, v]; cl = Dict(m => mean(skipmissing(dm[dm.mes .== m, v])) for m in 1:12)
        xf = [ismissing(a) ? cl[m] : a for (a, m) in zip(x, dm.mes[sel])]
        y = v == :prec ? log1p.(Float64.(xf)) : Float64.(xf)
        s = stl_simples(y)
        push!(info, (variavel = String(v), meses = length(x), preenchidos = count(ismissing, x),
                     tendencia_ini = s.tendencia[1], tendencia_fim = s.tendencia[end], var_sazonal = var(s.sazonal) / var(y)))
        out[!, Symbol(v, :_tend)] = vcat(fill(missing, count(.!sel)), s.tendencia)
        out[!, Symbol(v, :_saz)] = vcat(fill(missing, count(.!sel)), s.sazonal)
    end
    inf = DataFrame(info); mostra(inf; titulo = "STL: meses preenchidos e variação da tendência")
    salva_ptab(out, "stl_componentes"); salva_ptab(inf, "stl_resumo")
    for r in eachrow(inf)
        r.variavel == "prec" && continue
        k = r.variavel == "tmax" ? "Tmax" : "Tmin"
        registra!("STLDelta" * k, r.tendencia_fim - r.tendencia_ini; dig = 1)
        registra!("STLPreench" * k, r.preenchidos); registra!("STLMeses" * k, r.meses)
    end
    return inf
end

# ============================================================================= teleconexões
"Índice mensal PSL (ano + 12 colunas; −99,9x = ausente) → Dict((ano, mês) => valor)."
function le_indice_psl(nome)
    tab = Dict{Tuple{Int, Int}, Float64}()
    for l in readlines(joinpath(RAIZ, "data", "raw", "indices", nome * ".data"))[2:end]
        c = split(strip(l)); length(c) == 13 || continue
        a = tryparse(Int, c[1]); isnothing(a) && continue
        for m in 1:12
            v = parse(Float64, c[m+1]); v > -99 && (tab[(a, m)] = v)
        end
    end
    tab
end

"Média nov–abr do ano hidrológico a (nov/dez de a, jan–abr de a + 1)."
media_novabr(tab, a) = (ks = ((a, 11), (a, 12), (a + 1, 1), (a + 1, 2), (a + 1, 3), (a + 1, 4));
                        all(haskey(tab, k) for k in ks) ? mean(tab[k] for k in ks) : missing)

"""
Correlação da chuva anual com ONI e com o gradiente meridional do Atlântico
(TNA − TSA), ambos médios em nov–abr. Como a chuva e a TSM têm tendência, a
correlação é calculada também na série sem tendência (PIT sob a Gamma não
estacionária, levada à escala normal): uma correlação que só existe com a
tendência seria espúria. Também: P(ano seco | fase do ENSO), com seco = tercil
inferior e fases ONI ≥ 0,5 (El Niño), ≤ −0,5 (La Niña), demais (neutro).
"""
function bloco_teleconexoes(L, E, N)
    titulo("Teleconexões: ONI e gradiente TNA − TSA")
    oni = le_indice_psl("oni"); tna = le_indice_psl("tna"); tsa = le_indice_psl("tsa")
    v = sort(E.v, :ano_hidro)
    tc(a) = (a - mean(v.ano_hidro)) / 10
    z = [quantile(Normal(), clamp(N.nt.p_em(p, tc(a)), 1e-6, 1 - 1e-6)) for (p, a) in zip(v.prec, v.ano_hidro)]
    df = DataFrame(ano_hidro = v.ano_hidro, prec = Float64.(v.prec), z_sem_tend = z,
                   oni = [media_novabr(oni, a) for a in v.ano_hidro],
                   grad = [(t = media_novabr(tna, a); u = media_novabr(tsa, a); ismissing(t) || ismissing(u) ? missing : t - u) for a in v.ano_hidro])
    linhas = NamedTuple[]
    for (ind, nome) in ((:oni, "ONI"), (:grad, "TNA − TSA"))
        for (resp, rot) in ((:prec, "chuva (com tendência)"), (:z_sem_tend, "chuva sem tendência (PIT)"))
            s = dropmissing(df[:, [ind, resp]])
            r = cor(s[!, ind], s[!, resp]); n = nrow(s)
            lo, hi = ic_pearson(r, n); tt = r * sqrt((n - 2) / (1 - r^2))
            ρs = corspearman(s[!, ind], s[!, resp])
            push!(linhas, (indice = nome, resposta = rot, n = n, r_pearson = r, ic_inf = lo, ic_sup = hi,
                           p = 2ccdf(TDist(n - 2), abs(tt)), rho_spearman = ρs))
        end
    end
    tab = DataFrame(linhas); tab.q_bh = bh_ajuste(tab.p)
    mostra(tab; titulo = "Correlação da chuva anual com ONI e TNA − TSA (nov–abr)")
    salva_ptab(tab, "teleconexoes")
    for r in eachrow(tab)
        k = (r.indice == "ONI" ? "Oni" : "Grad") * (startswith(r.resposta, "chuva sem") ? "Det" : "Bruta")
        registra_ic!("R" * k, r.r_pearson, r.ic_inf, r.ic_sup; dig = 2); registra_p!("Q" * k, r.q_bh)
    end
    # tendência do gradiente (TNA aquece mais que TSA?) — ajuda a interpretar a diferença bruta × sem tendência
    g = dropmissing(df[:, [:ano_hidro, :grad]])
    tg = analise_tendencia(g.ano_hidro, g.grad; nome = "TNA − TSA")
    registra_ic!("SenGrad", tg.sen_decada, tg.sen_ic_inf, tg.sen_ic_sup; dig = 2); registra_p!("PSenGrad", tg.p_mk_tfpw)

    # P(seco | fase do ENSO): seco = tercil inferior da chuva (bruto) e da PIT (sem tendência)
    corte = quantile(df.prec, 1 / 3)
    fase(o) = ismissing(o) ? missing : o >= 0.5 ? "El Niño" : o <= -0.5 ? "La Niña" : "neutro"
    df.fase = fase.(df.oni)
    cond = NamedTuple[]
    for f in ("El Niño", "neutro", "La Niña"), (crit, sel) in (("bruto", df.prec .<= corte), ("sem tendência", df.z_sem_tend .<= quantile(Normal(), 1 / 3)))
        idx = coalesce.(df.fase .== f, false)
        n = count(idx); k = count(sel[idx])
        b = Beta(0.5 + k, 0.5 + n - k)
        push!(cond, (fase = f, criterio = crit, n = n, secos = k, p_med = median(b), ic_inf = quantile(b, 0.025), ic_sup = quantile(b, 0.975)))
    end
    tcd = DataFrame(cond)
    mostra(tcd; titulo = "P(ano seco | fase do ENSO), posterior Beta (priori de Jeffreys)")
    salva_ptab(tcd, "seco_por_fase_enso")
    for r in eachrow(tcd)
        k = (r.fase == "El Niño" ? "Nino" : r.fase == "La Niña" ? "Nina" : "Neutro") * (r.criterio == "bruto" ? "" : "Det")
        registra_ic!("PSeco" * k, r.p_med, r.ic_inf, r.ic_sup; dig = 2); registra!("NFase" * k, r.n); registra!("NSecos" * k, r.secos)
    end
    salva_ptab(df, "teleconexoes_anual")
    return (tab = tab, cond = tcd, df = df)
end

# ============================================================================= execução
function executa_paper()
    L = bloco_dados()
    H = bloco_homogeneidade(L)
    V = bloco_vizinhas(L)
    E = bloco_referencia(L)
    R = bloco_retorno(L, E)
    B = bloco_bayes(L, E)
    N = bloco_nao_estacionario(L, E, V)
    I = bloco_indices(L, E)
    M = bloco_markov(L, I)
    S = bloco_sensibilidade(L, E, R)
    Z = bloco_p2_sazonalidade(L)
    X = bloco_p2_indices(L)
    HT = bloco_p2_homog_temp(L, X)
    QR = bloco_p2_quantilica(L)
    ST = bloco_p2_stl(L)
    T = bloco_teleconexoes(L, E, N)
    println("\nMacros: ", escreve_macros())
    return (L = L, H = H, V = V, E = E, R = R, B = B, N = N, I = I, M = M, S = S, Z = Z, X = X, HT = HT, QR = QR, ST = ST, T = T)
end

using Serialization
const CACHE = joinpath(DIR_PTAB, "resultados.jls")

# `julia scripts/paper.jl`            → análises + figuras + resultados.tex
# `julia scripts/paper.jl figuras`    → só as figuras, a partir do cache da última execução
if abspath(PROGRAM_FILE) == @__FILE__
    if "figuras" in ARGS && isfile(CACHE)
        P = deserialize(CACHE)
    else
        P = executa_paper()
        serialize(CACHE, P)
    end
    gera_figuras_paper(P)
end
