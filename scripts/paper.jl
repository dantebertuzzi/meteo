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
const NBOOT_VALID = 500             # réplicas do p preditivo da validação e do IC da média em 2020/21

titulo(s) = println("\n", "="^length(s), "\n", s, "\n", "="^length(s))

# ============================================================================= dados
function bloco_dados()
    titulo("Série longa 82983")
    d = carrega_longo(); d_bruto = copy(d)              # d_bruto: sem o controle de qualidade da temperatura
    qc = qc_temperatura!(d)
    @printf("CQ da temperatura: segmento da Tmín removido (%d dias); > 4 DP: Tmáx %d, Tmín %d, Tméd %d dias\n",
            qc.segmento_tmin, qc.tmax, qc.tmin, qc.tmed)
    registra!("NDiasSegTmin", qc.segmento_tmin); registra!("NQCTmax", qc.tmax); registra!("NQCTmin", qc.tmin)
    dm = mensal_longo(d; ref = REF)
    dh, clim_anual = hidro_longo(d; ref = REF)
    v = dh[dh.valido, :]
    @printf("Registro diário: %s → %s | anos hidrológicos válidos: %d de %d | normal 1991–2020: %.1f mm\n",
            minimum(d.data), maximum(d.data), nrow(v), nrow(dh), clim_anual)
    registra!("NanosValidos", nrow(v)); registra!("NanosHidro", nrow(dh))
    registra!("PrimeiroAnoValido", v.rotulo[1]); registra!("UltimoAnoValido", v.rotulo[end])
    registra!("NormalChuva", clim_anual; dig = 0)
    registra!("MediaChuvaValidos", mean(v.prec); dig = 0); registra!("DPChuvaValidos", std(v.prec); dig = 0)
    # mesmo período da climatologia de Bebedouro (Embrapa Semiárido, 1963–2009; Teixeira, 2010)
    v09 = v[v.ano_hidro .<= 2008, :]
    registra!("MediaChuvaAteDoisMilNove", mean(v09.prec); dig = 0); registra!("NChuvaAteDoisMilNove", nrow(v09))
    salva_ptab(dh, "anos_hidrologicos_82983")
    salva_ptab(dm, "mensal_82983")
    return (d = d, d_bruto = d_bruto, dm = dm, dh = dh, clim_anual = clim_anual)
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
        registra!("QuebraSnht" * nm, res.quebra_snht_apos[i]); registra!("QuebraPettitt" * nm, res.quebra_pettitt_apos[i])
        registra_p!("PPettitt" * nm, res.p_pettitt[i])
    end
    i1213 = findfirst(==(2012), h.ano_hidro)
    !isnothing(i1213) && registra!("RazaoChuvaDozeTreze", h.razao_prec[i1213]; dig = 2)
    i2021 = findfirst(==(2021), h.ano_hidro)
    !isnothing(i2021) && registra!("RazaoChuvaVinteUm", h.razao_prec[i2021]; dig = 2)
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
    rs = DataFrame(resumo); rs.q_mk = bh_ajuste(rs.p_mk); rs.q_teste = bh_ajuste(rs.p_tfpw)
    pa = only(eachrow(rs[rs.estacao .== "Paulo Afonso", :]))
    registra!("SenVizPauloAfonso", pa.sen_decada; dig = 0); registra_p!("PVizPauloAfonso", pa.p_tfpw)
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
    registra!("NvizinhasTendSig", count(<(0.05), rs.q_teste))       # teste principal com BH
    # tendência da referência regional (sinal climático comum)
    tr = analise_tendencia(Q.ano, Q.referencia; nome = "referência regional")
    registra_ic!("SenReferencia", tr.sen_decada, tr.sen_ic_inf, tr.sen_ic_sup; dig = 0); registra_p!("PMKReferencia", tr.p_mk)
    mostra(DataFrame([tq, tr]); titulo = "Tendência: razão Q (inomogeneidade) × referência regional (clima)")
    registra_ic!("SenQ", 100tq.sen_decada, 100tq.sen_ic_inf, 100tq.sen_ic_sup; dig = 1)
    te = analise_tendencia(Q.ano, Q.candidata; nome = "82983 nos anos da referência")
    registra_ic!("SenEstAnosRef", te.sen_decada, te.sen_ic_inf, te.sen_ic_sup; dig = 0)
    registra!("PctRefEst", 100tr.sen_decada / te.sen_decada; dig = 0)
    a8 = Q.ano .<= 2008
    registra!("SenRefAteOito", analise_tendencia(Q.ano[a8], Q.referencia[a8]).sen_decada; dig = 0)
    registra!("SenEstAteOito", analise_tendencia(Q.ano[a8], Q.candidata[a8]).sen_decada; dig = 0)
    registra_p!("PTesteQ", tq.p_mk_tfpw); registra_p!("PTesteReferencia", tr.p_mk_tfpw)

    # par a par: log(82983/vizinha) — a referência composta muda de composição quando
    # estações fecham (S. do Bonfim 2011, P. Afonso 2015); pares fixos evitam esse artefato.
    pares = DataFrame(map([r for r in eachrow(info) if r.usada]) do r
        V = viz[r.estacao]; an = sort([a for a in keys(cand) if haskey(V, a)])
        lq = [log(cand[a] / V[a]) for a in an]
        t = analise_tendencia(an, lq; nome = r.estacao); s = snht_teste(Float64.(an), lq)
        (vizinha = r.estacao, rho_dif = r.rho_dif, n = length(an), inicio = rotulo_hidro(an[1]), fim = rotulo_hidro(an[end]),
         sen_pct_decada = 100(exp(t.sen_decada) - 1), ic_inf = 100(exp(t.sen_ic_inf) - 1), ic_sup = 100(exp(t.sen_ic_sup) - 1),
         p_mk = t.p_mk, p_teste = t.p_mk_tfpw, p_snht = s.p, quebra_snht = rotulo_hidro(Int(s.ultimo_antes)))
    end)
    pares.q_snht = bh_ajuste(pares.p_snht)
    mostra(pares; titulo = "Razões par a par log(82983/vizinha): tendência relativa e quebra")
    salva_ptab(pares, "vizinhas_pares")
    registra!("NparesQuebraSig", count(<(0.05), pares.q_snht)); registra!("Npares", nrow(pares))
    # tendência relativa nos pares, com o teste principal e BH na família dos pares
    pares.q_teste = bh_ajuste(pares.p_teste)
    registra!("NparesTendSig", count(<(0.05), pares.q_teste))
    neg = pares[pares.sen_pct_decada .< -3, :]                       # pares com queda relativa apreciável
    registra!("NparesQuedaRel", nrow(neg)); registra!("ParesQuedaRel", join(neg.vizinha, " e "))
    registra!("SenParesQuedaMin", minimum(neg.sen_pct_decada); dig = 1); registra!("SenParesQuedaMax", maximum(neg.sen_pct_decada); dig = 1)
    registra_p!("QParesQuedaMin", minimum(neg.q_teste))
    chave_viz = Dict("Cabrobó" => "Cabrobo", "Remanso" => "Remanso", "Paulo Afonso" => "PauloAfonso",
                     "Senhor do Bonfim" => "Bonfim", "Ouricuri" => "Ouricuri", "Petrolândia" => "Petrolandia",
                     "Caraíba" => "Caraiba", "Floresta" => "Floresta")
    for r in eachrow(pares)
        k = chave_viz[r.vizinha]
        registra!("SenPar" * k, r.sen_pct_decada; dig = 1); registra_p!("QPar" * k, r.q_teste)
        registra!("RhoPar" * k, r.rho_dif; dig = 2)
    end
    # composição da referência: vizinhas que encerram antes do último ano da candidata e anos com uma só
    ultc = maximum(keys(cand))
    fim_viz = Dict(r.estacao => maximum(keys(viz[r.estacao])) for r in eachrow(info) if r.usada)
    registra!("NvizEncerram", count(<(ultc), values(fim_viz)))
    um = Q.ano[Q.n_viz .== 1]
    registra!("NanosRefUmaViz", length(um)); registra!("AnosRefUmaViz", join(rotulo_hidro.(um), ", "))
    semref = sort([a for a in keys(cand) if !(a in Q.ano)])
    registra!("NanosSemRef", length(semref)); registra!("AnosSemRef", join(rotulo_hidro.(semref), ", "))
    # o SNHT também rejeita sob tendência gradual sem quebra: simula a reta ajustada a Q com resíduos
    # normais e conta quantas vezes T₀ supera o valor crítico de 5 % da série i.i.d.
    tq_ = Float64.(Q.ano); X = hcat(ones(length(tq_)), tq_); βq = X \ Q.Q; σq = std(Q.Q .- X * βq)
    rngq = MersenneTwister(SEMENTE)
    rej = mean(snht(X * βq .+ σq .* randn(rngq, length(tq_))).T0 > sn.T95 for _ in 1:2000)
    registra!("PctSnhtSobTendencia", 100rej; dig = 0)
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
    registra!("NBootEstMedia", length(pm))                 # réplicas válidas de 2000
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

    registra!("TMinimoAkaike", 1 / sum(cmp.peso_akaike .* cmp.p_xmin); sig = 2)      # E1′ pelo método principal
    # empírico
    k1213 = count(<=(E.x1213), x); kmin = count(<=(E.xmin), x)
    # sensibilidade ao empate: 2011/12 fica 0,1 mm acima de 2012/13
    x1112 = only(E.v.prec[E.v.rotulo .== "2011/12"])
    registra!("ChuvaOnzeDoze", x1112; dig = 1); registra!("NanosAteOnzeDoze", count(<=(x1112), x))
    registra!("TEmpiricoOnzeDoze", (n + 1) / count(<=(x1112), x); sig = 2)
    registra!("NanosAteDozeTreze", k1213)
    registra!("TEmpiricoDozeTreze", (n + 1) / k1213; sig = 2)

    out = Dict{String, Any}()
    for (rot, x0) in (("DozeTreze", E.x1213), ("Minimo", E.xmin))
        pr = perfil_p_gamma(x, x0)
        bs = bootstrap_p(Gamma, x, x0; B = 10_000)
        rb = resumo_replicas(bs)
        registra!("NBootGamma" * rot, length(bs))           # réplicas válidas de 10 000
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
    return (cmp = cmp, n = n, out = out, p_med = p_med, pm = pm)
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
    for (k, a0, a1) in (("Lacerda", 1965, 2006), ("Teixeira", 1978, 2007), ("Gois", 1975, 2016))
        m_ = a0 .<= v.ano_hidro .<= a1
        registra!("SenEst" * k, analise_tendencia(v.ano_hidro[m_], v.prec[m_]).sen_decada; dig = 0)
    end
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
    plr = Dict{String, Float64}()
    for D in CANDIDATAS
        nome = string(nameof(D))
        ffix = only(f for f in ma.fits if f.D === D && !f.forma_var)
        fvar = only(f for f in ma.fits if f.D === D && f.forma_var)
        plr[nome] = ccdf(Chisq(1), max(2(fvar.ll - ffix.ll), 0.0))
        registra_p!("PFormaVar" * nome, plr[nome])
        registra!("EfeitoForma" * nome, 100(exp(fvar.c) - 1); dig = 1)
    end
    nomes_lr = collect(keys(plr)); qlr = bh_ajuste([plr[k] for k in nomes_lr])
    for (k, q) in zip(nomes_lr, qlr)
        registra_p!("QFormaVar" * k, q)
    end
    nfx = only(f for f in ma.fits if f.D === Normal && !f.forma_var)
    registra!("PctNormalNeg", 100nfx.F(0.0, tc(ult)) / nfx.F(E.x1213, tc(ult)); dig = 0)
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
        # avaliada no último ano válido (comparável ao cenário local) e, como extrapolação, em 2025
        pbr = ns_media_akaike_boot(y, t, E.x1213, [tc(ult), tc(2025)]; b = b_reg)
        p0 = mr.p_em(E.x1213, tc(2025)); pu0 = mr.p_em(E.x1213, tc(ult))
        registra_ic!("TRegionalMedia", 1 / p0, 1 / quantile(pbr[:, 2], 0.975), 1 / quantile(pbr[:, 2], 0.025); sig = 2)
        registra!("NBootRegionalB", size(pbr, 1))
        registra_ic!("TRegionalMediaUltB", 1 / pu0, 1 / quantile(pbr[:, 1], 0.975), 1 / quantile(pbr[:, 1], 0.025); sig = 2)
        for nh in HORIZONTES
            registra_ic!("ProbRegMediaN" * string(nh), 100prob_ao_menos_um(pu0, nh),
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
    # estimador principal (média de Akaike das oito candidatas) com o critério FCA ≤ 7,5 %: inclui os
    # anos chuvosos depois de 2020/21, onde o T do fim da série é mais sensível
    m75 = ns_media_akaike(y75, t75); u75 = tc(maximum(d75.ano_hidro))
    pb75 = ns_media_akaike_boot(y75, t75, E.x1213, [u75]); registra!("NBootFCASete", size(pb75, 1))
    registra_ic!("TNSMediaFCASeteUlt", 1 / m75.p_em(E.x1213, u75), 1 / quantile(pb75[:, 1], 0.975), 1 / quantile(pb75[:, 1], 0.025); sig = 2)
    registra!("TNSMediaFCASeteUltStd", 1 / m75.p_em(E.x1213, tc(ult)); sig = 2)          # mesmo modelo, avaliado em 2020/21
    registra!("UltimoAnoFCASete", rotulo_hidro(maximum(d75.ano_hidro)))
    for r in eachrow(d75[in.(d75.ano_hidro, Ref(novos)), :])
        registra!("FCA" * nome_macro(rotulo_hidro(r.ano_hidro)), 100r.fca; dig = 1)   # FCA (%) dos anos acrescentados
    end
    # --- faixa entre cenários para o fim da série: o último ano é o ponto de maior alavancagem da
    #     reta, então o resultado sai como faixa, e não como um ponto em 2020/21
    pdec = mean(ma.p_em(E.x1213, tc(a)) for a in 2011:ult)
    registra!("TNSMediaDecada", 1 / pdec; sig = 2)
    # todos os cenários avaliados no MESMO ano (o último válido): misturar 2020/21, 2023/24 e 2025
    # numa faixa compararia climas de anos diferentes de uma reta com tendência
    cenT = [1 / ma.p_em(E.x1213, tc(ult)), 1 / pdec, 1 / m75.p_em(E.x1213, tc(ult))]
    if !ismissing(b_reg)
        tru = 1 / mr.p_em(E.x1213, tc(ult)); push!(cenT, tru)
        registra!("TRegionalMediaUlt", tru; sig = 2)
    end
    registra!("TNSFaixaMin", minimum(cenT); sig = 2); registra!("TNSFaixaMax", maximum(cenT); sig = 2)

    # --- sensibilidade da tendência da chuva ao período: a série termina na década seca de 2010,
    #     então testam-se a remoção da seca de 2011/12–2016/17, o corte antes dela e o início em 1971
    sens = [("Todos os anos válidos", "Todos", trues(nrow(v))),
            ("Sem 2011/12–2016/17", "SemSeca", .!(2011 .<= v.ano_hidro .<= 2016)),
            ("Até 2010/11", "AteDezOnze", v.ano_hidro .<= 2010),
            ("A partir de 1971/72", "DesdeSetentaUm", v.ano_hidro .>= 1971)]
    tsens = DataFrame(map(sens) do (rot, k, m)
        tr = analise_tendencia(v.ano_hidro[m], v.prec[m]; nome = rot)
        if k != "Todos"
            registra_ic!("SenChuva" * k, tr.sen_decada, tr.sen_ic_inf, tr.sen_ic_sup; dig = 0)
            registra_p!("PTFPWChuva" * k, tr.p_mk_tfpw); registra!("NChuva" * k, tr.n)
        end
        (amostra = rot, n = tr.n, sen = tr.sen_decada, ic_inf = tr.sen_ic_inf, ic_sup = tr.sen_ic_sup,
         p_mk = tr.p_mk, p_tfpw = tr.p_mk_tfpw)
    end)
    registra_p!("PTFPWChuvaFCASete", tr75.p_mk_tfpw); registra!("NChuvaFCASete", tr75.n)
    push!(tsens, (amostra = "FCA ≤ 7,5 %", n = tr75.n, sen = tr75.sen_decada, ic_inf = tr75.sen_ic_inf,
                  ic_sup = tr75.sen_ic_sup, p_mk = tr75.p_mk, p_tfpw = tr75.p_mk_tfpw))
    # anos recentes: FCA ≤ 7,5 % e também 2022/23, cujo total observado, mesmo com falhas, já é
    # chuvoso (limite inferior do total real), de modo que incluí-lo subestima esse ano
    drec = sort(vcat(d75, L.dh[L.dh.ano_hidro .== 2022, :]), :ano_hidro)
    trr = analise_tendencia(drec.ano_hidro, Float64.(drec.prec_obs); nome = "FCA ≤ 7,5 % + 2022/23")
    registra_ic!("SenChuvaRecentes", trr.sen_decada, trr.sen_ic_inf, trr.sen_ic_sup; dig = 0)
    registra_p!("PTFPWChuvaRecentes", trr.p_mk_tfpw); registra!("NChuvaRecentes", trr.n)
    r22 = only(eachrow(L.dh[L.dh.ano_hidro .== 2022, :]))
    registra!("ChuvaObsVinteDois", r22.prec_obs; dig = 0); registra!("DiasAusVinteDois", r22.dias_aus)
    for a in (2021, 2023)
        registra!("ChuvaObs" * (a == 2021 ? "VinteUm" : "VinteTres"), only(L.dh.prec_obs[L.dh.ano_hidro .== a]); dig = 0)
    end
    registra!("ChuvaObsOitentaUm", only(L.dh.prec_obs[L.dh.ano_hidro .== 1981]); dig = 0)
    push!(tsens, (amostra = "FCA ≤ 7,5 % + 2022/23", n = trr.n, sen = trr.sen_decada, ic_inf = trr.sen_ic_inf,
                  ic_sup = trr.sen_ic_sup, p_mk = trr.p_mk, p_tfpw = trr.p_mk_tfpw))
    mostra(tsens; titulo = "Sensibilidade da tendência da chuva anual ao período (Sen por década)")
    salva_ptab(tsens, "sensibilidade_tendencia_chuva")
    # forma da tendência: reta × degrau único × episódio seco (duas datas), com erro normal, por AIC.
    # Os cortes são escolhidos nos dados, então o AIC subestima a penalidade: o p vem de bootstrap
    # paramétrico sob a reta, repetindo a busca dos cortes em cada réplica.
    ftabs = NamedTuple[]
    for (rot, k, aa, yy) in (("critério principal", "", v.ano_hidro, y),
                             ("FCA ≤ 7,5 % + 2022/23", "Rec", drec.ano_hidro, Float64.(drec.prec_obs)))
        f = forma_tendencia(aa, yy)
        registra!("AICChuvaLinear" * k, f.lin; dig = 1)
        registra!("AICChuvaDegrau" * k, f.deg; dig = 1); registra!("AnoDegrauChuva" * k, rotulo_hidro(round(Int, f.c_deg)))
        registra!("AICEpisodio" * k, f.ep; dig = 1)
        registra!("AnosEpisodio" * k, rotulo_hidro(round(Int, f.c_ep[1])) * "--" * rotulo_hidro(round(Int, f.c_ep[2])))
        registra_p!("PBootDegrau" * k, f.p_deg); registra_p!("PBootEpisodio" * k, f.p_ep)
        if k == "Rec"
            ao = sort(aa); yo = yy[sortperm(aa)]; nn = length(yo)
            Xr = hcat(ones(nn), ao); βr = Xr \ yo; σr = sqrt(sum(abs2, yo .- Xr * βr) / nn)
            Tn(μ, σ) = 1 / cdf(Normal(μ, σ), E.x1213)
            registra!("TFormaRetaUlt", Tn(βr[1] + βr[2] * ult, σr); sig = 2)
            registra!("TFormaRetaFim", Tn(βr[1] + βr[2] * maximum(ao), σr); sig = 2)
            g = ao .>= f.c_deg
            σd = sqrt((sum(abs2, yo[g] .- mean(yo[g])) + sum(abs2, yo[.!g] .- mean(yo[.!g]))) / nn)
            registra!("TFormaDegrau", Tn(mean(yo[g]), σd); sig = 2); registra!("MediaFormaDegrau", mean(yo[g]); dig = 0)
            e = f.c_ep[1] .<= ao .<= f.c_ep[2]
            σe = sqrt((sum(abs2, yo[e] .- mean(yo[e])) + sum(abs2, yo[.!e] .- mean(yo[.!e]))) / nn)
            registra!("TFormaEpDentro", Tn(mean(yo[e]), σe); sig = 2); registra!("MediaFormaEpDentro", mean(yo[e]); dig = 0)
            registra!("TFormaEpFora", Tn(mean(yo[.!e]), σe); sig = 2); registra!("MediaFormaEpFora", mean(yo[.!e]); dig = 0)
        end
        push!(ftabs, (amostra = rot, n = length(yy), AIC_reta = f.lin, AIC_degrau = f.deg,
                      corte_degrau = rotulo_hidro(round(Int, f.c_deg)), p_degrau = f.p_deg, AIC_episodio = f.ep,
                      episodio = rotulo_hidro(round(Int, f.c_ep[1])) * "–" * rotulo_hidro(round(Int, f.c_ep[2])), p_episodio = f.p_ep))
    end
    mostra(DataFrame(ftabs); titulo = "Forma da tendência (AIC; p por bootstrap paramétrico sob a reta)")
    salva_ptab(DataFrame(ftabs), "forma_tendencia_chuva")

    # --- probabilidade de ≥ 1 evento em n anos, por cenário
    cen = NamedTuple[]
    g = fit_mle(Gamma, y); p_est = cdf(g, E.x1213)
    for nh in HORIZONTES
        congel = [1 - (1 - p25[i])^nh for i in 1:nd]
        continua = [prob_n_variavel(j -> pdraw(i, ult + j), nh) for i in 1:nd]
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
    return (tend = tend, tab = tab, oni = oni, nt = nt, nf = nf, po = po, cen = tcen, tT = tT, ma = ma, rhat = rh,
            pbm = pbm, pns_ult = ma.p_em(E.x1213, tc(ult)))
end

# ============================================================================= P1: validação fora da amostra
"""
Ajusta os modelos só com os anos válidos até `corte` e compara o número esperado de anos com
chuva ≤ limiar nos anos válidos seguintes com o observado. Dois limiares: E1 (2012/13), que tem
viés de seleção por ser um ano do próprio período de teste, e o percentil 10 do período de
ajuste, que não depende do período de teste. A tendência dos modelos não estacionários é
extrapolada para os anos de teste.
"""
function bloco_validacao(L, E; corte = 2005)
    titulo("Validação fora da amostra (ajuste até $(rotulo_hidro(corte)))")
    v = sort(E.v, :ano_hidro); y = Float64.(v.prec); a = v.ano_hidro
    tc(x) = (x - mean(a)) / 10
    tr = a .<= corte; te = .!tr
    ytr, ttr = y[tr], tc.(a[tr]); tte = tc.(a[te]); yte = y[te]
    registra!("AnoCorteValidacao", rotulo_hidro(corte)); registra!("NTreino", count(tr)); registra!("NTeste", count(te))
    registra!("PrimeiroAnoTeste", rotulo_hidro(minimum(a[te]))); registra!("UltimoAnoTeste", rotulo_hidro(maximum(a[te])))
    # modelos ajustados ao período de treino
    est = Dict(string(nameof(D)) => ajusta_mle(D, ytr) for D in CANDIDATAS)
    aic = Dict(k => 2length(params(d)) - 2loglikelihood(d, ytr) for (k, d) in est)
    amin = minimum(values(aic)); w = Dict(k => exp(-(x - amin) / 2) for (k, x) in aic); sw = sum(values(w))
    ma = ns_media_akaike(ytr, ttr)
    gfix = only(f for f in ma.fits if f.D === Gamma && !f.forma_var)
    nfix = only(f for f in ma.fits if f.D === Normal && !f.forma_var)
    modelos = [("EstGamma", "Gamma estacionária"), ("EstMedia", "Média de Akaike estacionária"),
               ("NSGamma", "Gamma não estacionária (CV fixo)"), ("NSNormal", "Normal não estacionária (DP fixo)"),
               ("NSMedia", "Média de Akaike não estacionária (8)")]
    pfun = Dict("EstGamma" => (x0, t) -> cdf(est["Gamma"], x0),
                "EstMedia" => (x0, t) -> sum(w[k] * cdf(est[k], x0) for k in keys(est)) / sw,
                "NSGamma" => (x0, t) -> gfix.F(x0, t),
                "NSNormal" => (x0, t) -> nfix.F(x0, t),
                "NSMedia" => (x0, t) -> ma.p_em(x0, t))
    p10 = quantile(ytr, 0.10)
    registra!("LimiarPdez", p10; dig = 0)
    linhas = NamedTuple[]
    x1112 = only(E.v.prec[E.v.rotulo .== "2011/12"])
    for (rotL, x0) in (("Eum", E.x1213), ("Pdez", p10), ("Eumb", x1112))
        obs = count(<=(x0), yte)
        registra!("Obs" * rotL, obs)
        for (k, nome) in modelos
            ps = [pfun[k](x0, t) for t in tte]
            esp = sum(ps); pk = prob_ao_menos_k(ps, obs)
            registra!("Esp" * rotL * k, esp; dig = 2)
            registra_p!("PObs" * rotL * k, pk)
            push!(linhas, (limiar = rotL, x0 = x0, modelo = nome, esperado = esp, observado = obs, p_ao_menos_obs = pk))
        end
    end
    # p preditivo: o p acima trata os parâmetros como conhecidos. Aqui eles variam por bootstrap
    # dos anos de treino (pares ano–chuva), e o p é a média de P(N ≥ obs | θ*) entre as réplicas.
    # O limiar P10 fica fixo no valor do treino original.
    rngv = MersenneTwister(SEMENTE); ntr = length(ytr)
    acum = Dict((rotL, k) => Float64[] for rotL in ("Eum", "Pdez", "Eumb") for (k, _) in modelos)
    for _ in 1:NBOOT_VALID
        i = rand(rngv, 1:ntr, ntr); yb = ytr[i]; tb = ttr[i]
        res = try
            (Dict(string(nameof(D)) => ajusta_mle(D, yb) for D in CANDIDATAS), ns_media_akaike(yb, tb))
        catch
            nothing
        end
        isnothing(res) && continue
        estb, mab = res
        aicb = Dict(k => 2length(params(d)) - 2loglikelihood(d, yb) for (k, d) in estb)
        mb = minimum(values(aicb)); wb = Dict(k => exp(-(x - mb) / 2) for (k, x) in aicb); swb = sum(values(wb))
        gfb = only(f for f in mab.fits if f.D === Gamma && !f.forma_var)
        nfb = only(f for f in mab.fits if f.D === Normal && !f.forma_var)
        pb = Dict("EstGamma" => (x0, t) -> cdf(estb["Gamma"], x0),
                  "EstMedia" => (x0, t) -> sum(wb[k] * cdf(estb[k], x0) for k in keys(estb)) / swb,
                  "NSGamma" => (x0, t) -> gfb.F(x0, t), "NSNormal" => (x0, t) -> nfb.F(x0, t),
                  "NSMedia" => (x0, t) -> mab.p_em(x0, t))
        for (rotL, x0) in (("Eum", E.x1213), ("Pdez", p10), ("Eumb", x1112)), (k, _) in modelos
            push!(acum[(rotL, k)], prob_ao_menos_k([pb[k](x0, t) for t in tte], count(<=(x0), yte)))
        end
    end
    registra!("NBootValid", length(acum[("Pdez", "EstGamma")]))
    tab = DataFrame(linhas)
    tab.p_preditivo = [mean(acum[(r.limiar, k)]) for (r, (k, _)) in zip(eachrow(tab), Iterators.cycle(modelos))]
    for (rotL, k) in keys(acum)
        registra_p!("PPredObs" * rotL * k, mean(acum[(rotL, k)]))
    end
    mostra(tab; titulo = "Esperado × observado nos anos de teste (P = probabilidade de ≥ observado)")
    salva_ptab(tab, "validacao_fora_amostra")
    return tab
end

# ============================================================================= P1: anos chuvosos recentes
"""
Os anos de 2021/22 a 2023/24 ficaram fora do ajuste pelo critério de completude. Eles são um teste
fora da amostra *depois* do fim da série: sob a média de Akaike não estacionária ajustada até o
último ano válido, calcula-se P(Y ≥ total observado) em cada ano. O total observado é um limite
inferior (dias sem registro contam como zero), então essas probabilidades são limites superiores.
"""
function bloco_anos_recentes(L, E, N)
    titulo("Anos chuvosos recentes, fora do ajuste, contra o modelo não estacionário")
    v = sort(E.v, :ano_hidro); tc(a) = (a - mean(v.ano_hidro)) / 10
    ma = N.ma
    Fm(x, t) = sum(ma.w[j] * ma.fits[j].F(x, t) for j in eachindex(ma.fits))
    function mediana_modelo(t)                        # bisseção de F(x) = 0,5
        lo, hi = 0.0, 3000.0
        for _ in 1:60
            m = (lo + hi) / 2; Fm(m, t) < 0.5 ? (lo = m) : (hi = m)
        end
        (lo + hi) / 2
    end
    # todos os anos depois do último válido com total observado (critério fixado antes do resultado),
    # e não só os chuvosos: 2024/25 entra, embora tenha FCA alta, como 2022/23
    RECENTES = (2021, 2022, 2023, 2024)
    sufixo = Dict(2021 => "VinteUm", 2022 => "VinteDois", 2023 => "VinteTres", 2024 => "VinteQuatro")
    linhas = NamedTuple[]
    for a in RECENTES
        r = only(eachrow(L.dh[L.dh.ano_hidro .== a, :]))
        x = Float64(r.prec_obs); p = 1 - Fm(x, tc(a)); med = mediana_modelo(tc(a))
        o = get(N.oni, a, missing)
        registra!("PRec" * sufixo[a], 100p; dig = 1); registra!("MedianaModelo" * sufixo[a], med; dig = 0)
        ismissing(o) || registra!("ONIRec" * sufixo[a], o; dig = 1)
        registra!("ChuvaObsRec" * sufixo[a], x; dig = 0); registra!("FCARec" * sufixo[a], 100r.fca; dig = 1)
        push!(linhas, (ano = rotulo_hidro(a), chuva_obs = x, dias_aus = r.dias_aus, fca = r.fca,
                       mediana_modelo = med, p_ao_menos_obs = p, oni_nov_abr = o))
    end
    tab = DataFrame(linhas)
    pj = prod(tab.p_ao_menos_obs)
    registra!("PRecConjuntaInv", 1 / pj; sig = 1)
    # versão preditiva: o p acima trata os parâmetros como conhecidos, justamente no ponto de
    # maior alavancagem da reta. Aqui os parâmetros variam pelo bootstrap dos anos (pares
    # ano–chuva), e P(Y ≥ obs) é a média entre réplicas. A conjunta é a média do produto (os três
    # anos compartilham os parâmetros de cada réplica), e não o produto das médias.
    yv = Float64.(v.prec); tv = tc.(v.ano_hidro); nv = length(yv)
    rngr = MersenneTwister(SEMENTE); pr = Vector{Vector{Float64}}()
    xs_rec = [(Float64(only(L.dh.prec_obs[L.dh.ano_hidro .== a])), tc(a)) for a in RECENTES]
    for _ in 1:NBOOT_VALID
        i = rand(rngr, 1:nv, nv)
        mb = try ns_media_akaike(yv[i], tv[i]) catch; continue end
        push!(pr, [1 - mb.p_em(x, tt) for (x, tt) in xs_rec])
    end
    P = permutedims(reduce(hcat, pr))
    tab.p_preditivo = vec(mean(P, dims = 1))
    for (j, a) in enumerate(RECENTES)
        registra!("PRecPred" * sufixo[a], 100tab.p_preditivo[j]; dig = 1)
    end
    pjp = mean(prod(P, dims = 2))
    registra!("PRecConjuntaPredInv", 1 / pjp; sig = 1); registra!("NBootRec", size(P, 1))
    # o produto de três probabilidades de cauda não é um p-valor (qualquer produto específico é
    # pequeno); a combinação de Fisher, −2 Σ log pᵢ ~ χ²_{2k} sob H₀, é a medida calibrada de surpresa
    fisher(ps) = ccdf(Chisq(2length(ps)), -2sum(log, ps))
    registra!("PRecFisherInv", 1 / fisher(tab.p_ao_menos_obs); sig = 2)
    registra!("PRecFisherPredInv", 1 / fisher(tab.p_preditivo); sig = 2)
    # bilateral: a direção (anos chuvosos) foi notada depois de ver os dados, então a versão que não
    # escolhe a cauda é a de referência: pᵢ ← 2 min(pᵢ, 1 − pᵢ)
    bil(ps) = 2 .* min.(ps, 1 .- ps)
    registra!("PRecFisherBilInv", 1 / fisher(bil(tab.p_ao_menos_obs)); sig = 2)
    registra!("PRecFisherBilPredInv", 1 / fisher(bil(tab.p_preditivo)); sig = 2)
    registra!("NAnosRec", length(RECENTES))
    # média do modelo no último ano válido, com IC do bootstrap dos anos (normal com DP fixo, a
    # candidata de menor AIC; serve para comparar com a normal oficial de 1991–2020)
    nfix = only(f for f in ma.fits if f.D === Normal && !f.forma_var)
    ult = maximum(v.ano_hidro); y = Float64.(v.prec); t = tc.(v.ano_hidro)
    rng = MersenneTwister(SEMENTE); mb = Float64[]
    for _ in 1:NBOOT_VALID
        i = rand(rng, 1:length(y), length(y))
        f = try ns_candidata(Normal, y[i], t[i]) catch; continue end
        push!(mb, exp(f.θ[2] + f.θ[3] * tc(ult)))
    end
    registra!("NBootMediaNormal", length(mb))
    registra_ic!("MediaNormalUlt", exp(nfix.θ[2] + nfix.θ[3] * tc(ult)), quantile(mb, 0.025), quantile(mb, 0.975); dig = 0)
    mostra(tab; titulo = "P(Y ≥ observado) sob a média de Akaike ajustada até $(rotulo_hidro(ult))")
    salva_ptab(tab, "anos_recentes_modelo")
    return tab
end

# ============================================================================= M4: pluviômetros independentes (ANA)
const PREC_MAX_DIA = 400.0          # mm/dia: acima disso, o valor diário é tratado como erro (como no BDMEP)
"""
Verificação independente da tendência da chuva com pluviômetros da rede da ANA (outros
operadores, outros instrumentos). Cada posto passa pelo mesmo critério de completude da 82983
(FCA ≤ 5 %, com a climatologia diária do próprio posto em 1991–2020). Entram os postos com
≥ `min_comum` anos válidos em comum com a 82983. Nos anos comuns, compara-se a tendência dos
dois e testa-se a razão Q = 82983/posto: tendência (Sen, MK), quebra (SNHT) e diferença de
média antes e depois de 2003/04, quando a estação automática foi instalada no sítio da 82983.
"""
function bloco_ana(L, E; min_comum = 30, ano_a307 = 2003)
    titulo("Pluviômetros independentes da rede da ANA")
    ref = Dict(a => Float64(p) for (a, p) in zip(E.v.ano_hidro, E.v.prec))
    linhas = NamedTuple[]
    for cod in sort(collect(keys(POSTOS_ANA)))
        arq = joinpath(RAIZ, "data", "raw", "ana", "ana_$(cod).csv"); isfile(arq) || continue
        info = POSTOS_ANA[cod]
        d = CSV.read(arq, DataFrame; missingstring = "")
        d.prec = allowmissing(Float64.(coalesce.(d.prec, NaN))); d.prec[isnan.(d.prec)] .= missing
        d.prec[coalesce.(d.prec .> PREC_MAX_DIA, false)] .= missing      # mesmo limite de plausibilidade do BDMEP
        d.ano = year.(d.data); d.mes = month.(d.data); d.ano_hidro = ano_hidro.(d.data)
        dh, _ = hidro_longo(d; ref = REF)
        v = dh[dh.valido, :]
        com = sort([a for a in v.ano_hidro if haskey(ref, a)])
        dist = dist_km(info.lat, info.lon, LAT_82983, LON_82983)
        push!(linhas, (codigo = cod, posto = info.nome, operador = info.resp, dist_km = dist, n_validos = nrow(v),
                       n_comum = length(com), usado = length(com) >= min_comum))
        length(com) >= min_comum || continue
        yp = [only(v.prec[v.ano_hidro .== a]) for a in com]; yr = [ref[a] for a in com]
        tp = theil_sen(Float64.(com), yp); tr = theil_sen(Float64.(com), yr)
        q = yr ./ yp; tq = theil_sen(Float64.(com), q)
        mq = (p = analise_tendencia(com, q).p_mk_tfpw,)                # teste principal na razão Q
        pposto = analise_tendencia(com, yp).p_mk_tfpw                    # teste principal no próprio posto
        # anos hidrológicos inteiramente consistidos (nível 2) pela ANA
        cons = Set(g.ano_hidro[1] for g in groupby(d, :ano_hidro) if all(==(2), coalesce.(g.nivel_consistencia, 1)))
        comc = [a for a in com if a in cons]
        sc = length(comc) >= 15 ? (theil_sen(Float64.(comc), [only(v.prec[v.ano_hidro .== a]) for a in comc]).declive,
                                    theil_sen(Float64.(comc), [ref[a] for a in comc]).declive,
                                    analise_tendencia(comc, [ref[a] / only(v.prec[v.ano_hidro .== a]) for a in comc]).p_mk_tfpw,
                                    maximum(comc)) : (NaN, NaN, NaN, missing)
        # anos antes da quebra (SNHT) da razão com o posto
        sn = snht_teste(Float64.(com), q)
        c11 = [a for a in com if a <= Int(sn.ultimo_antes)]
        s11 = length(c11) >= 8 ? (10theil_sen(Float64.(c11), [only(v.prec[v.ano_hidro .== a]) for a in c11]).declive,
                                  10theil_sen(Float64.(c11), [ref[a] for a in c11]).declive) : (NaN, NaN)
        antes = q[com .< ano_a307]; depois = q[com .>= ano_a307]
        pw = length(depois) >= 3 ? pvalue(UnequalVarianceTTest(antes, depois)) : NaN
        k = info.chave
        registra!("Dist" * k, dist; dig = 0); registra!("NComum" * k, length(com))
        registra_ic!("SenAna" * k, 10tp.declive, 10tp.ic_inf, 10tp.ic_sup; dig = 0)
        registra!("SenRef" * k, 10tr.declive; dig = 0)
        registra!("SenQAna" * k, 1000tq.declive; dig = 1); registra_p!("PQAna" * k, mq.p)
        registra_p!("PSnhtAna" * k, sn.p); registra!("QuebraAna" * k, rotulo_hidro(Int(sn.ultimo_antes) + 1))
        registra!("QAntesAna" * k, mean(antes); dig = 2); registra!("QDepoisAna" * k, mean(depois); dig = 2)
        registra_p!("PAntesDepoisAna" * k, pw); registra_p!("PTesteAna" * k, pposto)
        registra!("SenAnaAteOnze" * k, s11[1]; dig = 0); registra!("SenRefAteOnze" * k, s11[2]; dig = 0)
        registra!("NComumAteOnze" * k, length(c11)); registra!("AnoAntesQuebraAna" * k, rotulo_hidro(Int(sn.ultimo_antes)))
        linhas[end] = merge(linhas[end], (p_posto = pposto, n_cons = length(comc), ult_cons = sc[4],
                                          razao_sen = tp.declive / tr.declive, razao_sen_cons = sc[1] / sc[2],
                                          p_Q_cons = sc[3], sen_posto = 10tp.declive, sen_posto_inf = 10tp.ic_inf, sen_posto_sup = 10tp.ic_sup,
                                          sen_82983 = 10tr.declive, sen_Q_pct = 1000tq.declive, p_Q_mk = mq.p,
                                          p_snht = sn.p, quebra_apos = Int(sn.ultimo_antes), Q_antes = mean(antes),
                                          Q_depois = mean(depois), p_antes_depois = pw))
    end
    tab = DataFrame([merge((p_posto = missing, n_cons = missing, ult_cons = missing, razao_sen = missing,
                            razao_sen_cons = missing, p_Q_cons = missing, sen_posto = missing, sen_posto_inf = missing, sen_posto_sup = missing, sen_82983 = missing,
                            sen_Q_pct = missing, p_Q_mk = missing, p_snht = missing, quebra_apos = missing,
                            Q_antes = missing, Q_depois = missing, p_antes_depois = missing), l) for l in linhas])
    sort!(tab, :dist_km)
    mostra(tab; titulo = "Postos da ANA: tendência nos anos comuns com a 82983 e razão Q = 82983/posto")
    salva_ptab(tab, "ana_postos")
    u = tab[tab.usado, :]
    registra!("NPostosAnaCandidatos", nrow(tab)); registra!("NPostosAna", nrow(u))
    registra!("MinComumAna", min_comum)
    registra!("NAnaNeg", count(<(0), u.sen_posto)); registra!("NAnaNegSig", count(<(0), u.sen_posto_sup))
    registra!("SenAnaMin", minimum(u.sen_posto); dig = 0); registra!("SenAnaMax", maximum(u.sen_posto); dig = 0)
    sg = u[u.sen_posto_sup .< 0, :]
    registra!("SenAnaSigMin", minimum(sg.sen_posto); dig = 0); registra!("SenAnaSigMax", maximum(sg.sen_posto); dig = 0)
    registra!("NAnaQSig", count(<(0.05), u.p_Q_mk)); registra!("NAnaAntesDepoisSig", count(<(0.05), u.p_antes_depois))
    registra!("NAnaNegSigTeste", count(r -> r.sen_posto < 0 && r.p_posto < 0.05, eachrow(u)))
    # quanto da tendência da 82983 os postos compartilham: razão dos declives nos anos comuns
    registra!("RazaoSenAnaMediana", 100median(u.razao_sen); dig = 0)
    uc = u[.!isnan.(u.razao_sen_cons), :]
    registra!("RazaoSenAnaConsMediana", 100median(uc.razao_sen_cons); dig = 0)
    registra!("NAnaCons", nrow(uc)); registra!("NAnaConsMenor", count(<(1), uc.razao_sen_cons))
    registra!("NAnaConsQSig", count(<(0.05), uc.p_Q_cons))
    registra!("UltimoAnoConsAna", rotulo_hidro(maximum(skipmissing(uc.ult_cons))))
    return tab
end

"""
Período de retorno de um ano tão seco quanto E1 nos postos da ANA usados na
verificação, com o mesmo esquema da 82983 (médias por pesos de Akaike, estacionária
e não estacionária, IC por bootstrap dos anos), cada posto com todos os seus anos
válidos. O limiar é relativo: E1 × (média do posto / média da 82983) nos anos
comuns, porque o 2012/13 do posto mais próximo (Juazeiro) não é válido. O 2012/13
do próprio posto, quando válido, entra como sensibilidade. O valor não
estacionário é avaliado no mesmo ano da 82983 (último ano válido dela).
"""
function bloco_ana_retorno(L, E, AN; B = 1000)
    titulo("Período de retorno de E1 nos postos da ANA")
    ref = Dict(a => Float64(p) for (a, p) in zip(E.v.ano_hidro, E.v.prec))
    ult = maximum(E.v.ano_hidro)
    function media_akaike(x, x0)
        fits = [ajusta_mle(D, x) for D in CANDIDATAS]
        aic = [2length(params(f)) - 2loglikelihood(f, x) for f in fits]
        w = exp.(-(aic .- minimum(aic)) ./ 2); w ./= sum(w)
        return sum(w .* [cdf(f, x0) for f in fits])
    end
    linhas = NamedTuple[]
    for r in eachrow(AN[AN.usado, :])
        info = POSTOS_ANA[r.codigo]; k = info.chave
        d = CSV.read(joinpath(RAIZ, "data", "raw", "ana", "ana_$(r.codigo).csv"), DataFrame; missingstring = "")
        d.prec = allowmissing(Float64.(coalesce.(d.prec, NaN))); d.prec[isnan.(d.prec)] .= missing
        d.prec[coalesce.(d.prec .> PREC_MAX_DIA, false)] .= missing      # mesmo limite de plausibilidade do BDMEP
        d.ano = year.(d.data); d.mes = month.(d.data); d.ano_hidro = ano_hidro.(d.data)
        dh, _ = hidro_longo(d; ref = REF)
        # mesmo início da 82983: antes de 1963 há segmentos inomogêneos (p.ex. Lagoa Grande
        # registra 1 250–1 550 mm/ano em 1955–61, cerca do triplo do restante da série)
        vtodos = sort(dh[dh.valido .& (dh.ano_hidro .>= minimum(E.v.ano_hidro)), :], :ano_hidro)
        # os modelos usam só os anos até o último válido da 82983, para que os T sejam comparáveis
        # aos dela; as tendências "até o fim" usam todos os anos (vtodos)
        v = vtodos[vtodos.ano_hidro .<= ult, :]
        y = Float64.(v.prec); n = length(y); m = mean(v.ano_hidro)
        t = (v.ano_hidro .- m) ./ 10; tc(a) = (a - m) / 10
        com = [a for a in v.ano_hidro if haskey(ref, a)]
        razao = mean(y[v.ano_hidro .∈ Ref(Set(com))]) / mean(ref[a] for a in com)
        x0 = E.x1213 * razao
        # estacionário
        pe = media_akaike(y, x0)
        rng = MersenneTwister(SEMENTE); pbe = Float64[]
        for _ in 1:B
            yb = y[rand(rng, 1:n, n)]
            p = try media_akaike(yb, x0) catch; continue end
            push!(pbe, p)
        end
        # não estacionário, no último ano válido da 82983 e no do próprio posto
        ma = ns_media_akaike(y, t)
        pn = ma.p_em(x0, tc(ult)); pn_posto = ma.p_em(x0, tc(maximum(v.ano_hidro)))
        pbn = ns_media_akaike_boot(y, t, x0, [tc(ult)]; B)
        tend = analise_tendencia(vtodos.ano_hidro, vtodos.prec; nome = info.nome)
        tult = analise_tendencia(v.ano_hidro, v.prec; nome = info.nome * " até " * rotulo_hidro(ult))
        # sensibilidade: limiar = 2012/13 do próprio posto, quando válido
        own = v.prec[v.ano_hidro .== 2012]
        x0_own = isempty(own) ? missing : Float64(only(own))
        Tn_own = ismissing(x0_own) ? missing : 1 / ma.p_em(x0_own, tc(ult))
        Te_own = ismissing(x0_own) ? missing : 1 / media_akaike(y, x0_own)
        l = (codigo = r.codigo, posto = info.nome, dist_km = r.dist_km, n = n,
             periodo = string(v.rotulo[1], "–", vtodos.rotulo[end]), razao_media = razao, limiar_mm = x0,
             sen_decada = tend.sen_decada, sen_inf = tend.sen_ic_inf, sen_sup = tend.sen_ic_sup, p_tfpw = tend.p_mk_tfpw,
             sen_ate_ult = tult.sen_decada, p_tfpw_ate_ult = tult.p_mk_tfpw,
             T_est = 1 / pe, T_est_inf = 1 / quantile(pbe, 0.975), T_est_sup = 1 / quantile(pbe, 0.025),
             T_ns = 1 / pn, T_ns_inf = 1 / quantile(pbn[:, 1], 0.975), T_ns_sup = 1 / quantile(pbn[:, 1], 0.025),
             T_ns_ult_posto = 1 / pn_posto, nb_est = length(pbe), nb_ns = size(pbn, 1), limiar_proprio = x0_own, T_est_proprio = Te_own, T_ns_proprio = Tn_own)
        push!(linhas, l)
        registra!("NAnaRet" * k, n); registra!("PeriodoAnaRet" * k, l.periodo); registra!("LimiarAna" * k, x0; dig = 0)
        registra_ic!("SenAnaTodos" * k, l.sen_decada, l.sen_inf, l.sen_sup; dig = 0); registra_p!("PTFPWAnaTodos" * k, l.p_tfpw)
        registra!("SenAnaAteUlt" * k, l.sen_ate_ult; dig = 0); registra_p!("PTFPWAnaAteUlt" * k, l.p_tfpw_ate_ult)
        registra_ic!("TEstAna" * k, l.T_est, l.T_est_inf, l.T_est_sup; sig = 2)
        registra_ic!("TNSAna" * k, l.T_ns, l.T_ns_inf, l.T_ns_sup; sig = 2)
        registra!("TNSAnaUltPosto" * k, l.T_ns_ult_posto; sig = 2)
        ismissing(Tn_own) || (registra!("TNSAnaProprio" * k, Tn_own; sig = 2); registra!("TEstAnaProprio" * k, Te_own; sig = 2))
    end
    tab = sort(DataFrame(linhas), :dist_km)
    mostra(tab; titulo = "Período de retorno de E1 (limiar relativo) nos postos da ANA")
    salva_ptab(tab, "ana_retorno")
    registra!("NAnaRetorno", nrow(tab))
    registra!("TEstAnaMin", minimum(tab.T_est); sig = 2); registra!("TEstAnaMax", maximum(tab.T_est); sig = 2)
    registra!("TNSAnaMin", minimum(tab.T_ns); sig = 2); registra!("TNSAnaMax", maximum(tab.T_ns); sig = 2)
    registra!("TNSAnaMediana", median(tab.T_ns); sig = 2)
    registra!("SenAnaTodosMediana", median(tab.sen_decada); dig = 0); registra!("SenAnaAteUltMediana", median(tab.sen_ate_ult); dig = 0)
    registra!("NAnaTodosNegSig", count(<(0), tab.sen_sup)); registra!("NAnaTodosPSig", count(<(0.05), tab.p_tfpw))
    registra!("UltimoAnoAna", rotulo_hidro(maximum(parse(Int, split(p, "–")[2][1:4]) for p in tab.periodo)))
    registra!("NAnaRetornoCai", count(tab.T_ns .< tab.T_est))
    registra!("NAnaProprio", count(!ismissing, tab.T_ns_proprio))
    registra!("NBootAnaRetorno", minimum(min.(tab.nb_est, tab.nb_ns)))      # réplicas válidas (mínimo entre postos)
    return tab
end

"""
Análise regional de frequência (índice de cheia com tendência): a 82983 e os postos
da ANA usados na verificação, desde 1963/64, num modelo com distribuição padronizada
e tendência comuns e um índice μᵢ por estação (μ = 1 na 82983). Médias por pesos de
Akaike das candidatas, estacionárias e não estacionárias; IC por bootstrap por ano.
A homogeneidade da tendência é avaliada contra um modelo com inclinação por estação.
"""
function bloco_regional(L, E, AN; B = 400)
    titulo("Análise regional de frequência (82983 + postos da ANA)")
    v82 = sort(E.v, :ano_hidro)
    y = Float64.(v82.prec); ano = collect(v82.ano_hidro); s = fill(1, nrow(v82))
    nomes = ["82983"]; chaves = ["Est"]
    u = sort(AN[AN.usado, :], :dist_km)
    for (k, r) in enumerate(eachrow(u))
        d = CSV.read(joinpath(RAIZ, "data", "raw", "ana", "ana_$(r.codigo).csv"), DataFrame; missingstring = "")
        d.prec = allowmissing(Float64.(coalesce.(d.prec, NaN))); d.prec[isnan.(d.prec)] .= missing
        d.prec[coalesce.(d.prec .> PREC_MAX_DIA, false)] .= missing      # mesmo limite de plausibilidade do BDMEP
        d.ano = year.(d.data); d.mes = month.(d.data); d.ano_hidro = ano_hidro.(d.data)
        dh, _ = hidro_longo(d; ref = REF)
        v = dh[dh.valido .& (minimum(v82.ano_hidro) .<= dh.ano_hidro .<= maximum(v82.ano_hidro)), :]   # mesmos anos da 82983
        append!(y, Float64.(v.prec)); append!(ano, v.ano_hidro); append!(s, fill(k + 1, nrow(v)))
        push!(nomes, POSTOS_ANA[r.codigo].nome); push!(chaves, POSTOS_ANA[r.codigo].chave)
    end
    S = length(nomes); m = mean(v82.ano_hidro); x = (ano .- m) ./ 10; tc(a) = (a - m) / 10
    ult = maximum(v82.ano_hidro); fim = maximum(ano)
    mn = regional_media_akaike(y, x, s; nsites = S)
    me = regional_media_akaike(y, x, s; nsites = S, estacionario = true)
    sort!(mn.tab, :AIC); mostra(mn.tab; titulo = "Candidatas regionais não estacionárias")
    salva_ptab(mn.tab, "regional_candidatas")
    bs = regional_media_akaike_boot(y, x, s, ano, E.x1213, tc(ult); nsites = S, B)
    pe = me.p_em(E.x1213, 0.0); pn = mn.p_em(E.x1213, tc(ult))
    registra!("NRegObs", length(y)); registra!("NRegEstacoes", S); registra!("NBootReg", length(bs.ns))
    registra_ic!("TRegEst", 1 / pe, 1 / quantile(bs.est, 0.975), 1 / quantile(bs.est, 0.025); sig = 2)
    registra_ic!("TRegNS", 1 / pn, 1 / quantile(bs.ns, 0.975), 1 / quantile(bs.ns, 0.025); sig = 2)
    registra!("TRegNSFim", 1 / mn.p_em(E.x1213, tc(fim)); sig = 2); registra!("UltimoAnoReg", rotulo_hidro(fim))
    # tendência regional (média dos efeitos ponderada pelos pesos) e homogeneidade da tendência
    registra!("RegEfeito", sum(mn.w .* [100(exp(f.b) - 1) for f in mn.fits]); dig = 1)
    imelhor = argmax(mn.w); fm = mn.fits[imelhor]
    registra!("RegMelhor", lowercase(string(nameof(fm.D)))); registra!("RegMelhorForma", fm.forma_var ? "variável" : "fixa")
    registra!("RegMelhorPeso", mn.w[imelhor]; dig = 2)
    fc = regional_candidata(fm.D, y, x, s; nsites = S, forma_var = fm.forma_var)
    fs = regional_candidata(fm.D, y, x, s; nsites = S, forma_var = fm.forma_var, b_por_site = true)
    registra!("DAICRegSite", (2fs.npar - 2fs.ll) - (2fc.npar - 2fc.ll); dig = 1)
    ef = 100 .* (exp.(fs.b) .- 1)
    for (k, c) in enumerate(chaves)
        registra!("RegEfeitoSite" * c, ef[k]; dig = 1); registra!("RegIndice" * c, fm.μ[k]; dig = 2)
        registra!("NReg" * c, count(==(k), s))
    end
    registra!("RegEfeitoSiteMin", minimum(ef); dig = 1); registra!("RegEfeitoSiteMax", maximum(ef); dig = 1)
    # largura relativa dos intervalos (razão sup/inf de T), regional × 82983 sozinha
    tab = DataFrame(estacao = nomes, n = [count(==(k), s) for k in 1:S], indice_mu = fm.μ, efeito_pct_proprio = ef)
    mostra(tab; titulo = "Estações na análise regional")
    salva_ptab(tab, "regional_estacoes")
    @printf("T regional: estacionário %.1f [%.1f; %.1f]; não estacionário (%s) %.1f [%.1f; %.1f]\n",
            1 / pe, 1 / quantile(bs.est, 0.975), 1 / quantile(bs.est, 0.025), rotulo_hidro(ult),
            1 / pn, 1 / quantile(bs.ns, 0.975), 1 / quantile(bs.ns, 0.025))
    return (tab = tab, mn = mn, me = me, boot = bs)
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
    dm.spei12_pp, _ = spei(dm.balanco, dm.mes, dm.ano; k = 12, ref = REF, pwm = :pp)     # PWM de Vicente-Serrano
    registra!("NRefJunSPEI", count(i -> dm.mes[i] == 6 && REF[1] <= dm.ano[i] <= REF[2] && !ismissing(dm.spei12[i]), 1:nrow(dm)))

    # Stagge et al. (2015): qualidade do ajuste pela normalidade do índice resultante
    sw = DataFrame(map([(:spi12, "SPI-12 (Gamma)"), (:spei12, "SPEI-12 (log-logística)"), (:spei12_gev, "SPEI-12 (GEV)")]) do (v, nm)
        ref = [z for (z, a) in zip(dm[!, v], dm.ano) if !ismissing(z) && REF[1] <= a <= REF[2]]
        jun = [z for (z, a, m) in zip(dm[!, v], dm.ano, dm.mes) if !ismissing(z) && REF[1] <= a <= REF[2] && m == 6]
        t = ShapiroWilkTest(Float64.(jun))                 # 1 valor por ano (junho), sem sobreposição
        (indice = nm, n_meses = count(!ismissing, dm[!, v]), media_ref = mean(ref), dp_ref = std(ref),
         W = t.W, p_shapiro = pvalue(t))
    end)
    mostra(sw; titulo = "Normalidade dos índices no período de referência")
    salva_ptab(sw, "indices_normalidade")

    # índice no fim do ano hidrológico (junho) = SPI-12/SPEI-12 do ano hidrológico
    jun = dm[dm.mes .== 6, [:ano, :spi12, :spei12, :spei12_gev, :spei12_pp]]
    jun.ano_hidro = jun.ano .- 1
    for (rot, a) in (("DozeTreze", 2012), ("DezesseisDezessete", 2016), ("OnzeDoze", 2011))
        r = jun[jun.ano_hidro .== a, :]
        nrow(r) == 1 || continue
        registra!("SPI" * rot, r.spi12[1]; dig = 2); registra!("SPEI" * rot, r.spei12[1]; dig = 2)
        registra!("SPEIpp" * rot, r.spei12_pp[1]; dig = 2)
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
        e3 = e[e.duracao .>= 3, :]
        T3 = 1 / (nrow(e3) / anos_cobertos * ccdf(fit_mle(Exponential, e3.severidade), r.severidade))
        registra!("TSeveridadeRefMinTres", T3; sig = 2); registra!("NeventosSPIMinTres", nrow(e3))
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
    # mesma análise com a normal não estacionária (DP constante), a candidata de menor AIC:
    # verifica que a conclusão não depende da Gamma usada para remover a tendência
    nn = ns_candidata(Normal, Float64.(v.prec), tc.(v.ano_hidro))
    est4 = [ismissing(p) ? missing : (nn.F(p, tc(a)) <= 1 / 3 ? "S" : "N") for (p, a) in zip(L.dh.prec, L.dh.ano_hidro)]
    m4 = markov_posterior(transicoes(L.dh.ano_hidro, est4))
    registra_ic!("PSSdetNormal", m4.pSS[2], m4.pSS[1], m4.pSS[3]; dig = 2)
    registra_ic!("PNSdetNormal", m4.pNS[2], m4.pNS[1], m4.pNS[3]; dig = 2)
    registra!("ProbPersistenciaDetNormal", m4.prob_dif_pos; dig = 2)
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
    tend.q_teste = bh_ajuste(tend.p_mk_tfpw)
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
    # R̄ é fixo (normal de 1991–2020): num ano mais seco, o início fica mais tarde por construção.
    # Controlando pelo total do ano, quanto sobra da tendência das datas?
    sd_ = dropmissing(s[s.estacao_definida, [:decada, :prec, :inicio, :duracao]])
    for (vv, nm) in ((:inicio, "Inicio"), (:duracao, "Duracao"))
        mm_ = lm(Term(vv) ~ Term(:decada) + Term(:prec), sd_)
        registra!("Sen" * nm * "ControleTotal", coef(mm_)[2]; dig = 1)
        registra_p!("P" * nm * "ControleTotal", coeftable(mm_).cols[4][2])
    end
    registra!("FracNovAbr", 100mean(s.frac_nov_abr); dig = 0)
    for (i, nm) in enumerate(("SI", "Centroide", "Conc", "Inicio", "Fim", "Duracao", "FracNovAbr"))
        r = tend[i, :]; dg = nm in ("SI", "Conc", "FracNovAbr") ? 3 : 1
        registra_ic!("Sen" * nm, r.sen_decada, r.sen_ic_inf, r.sen_ic_sup; dig = dg)
        registra_p!("PMK" * nm, r.p_mk); registra_p!("PTFPW" * nm, r.p_mk_tfpw); registra_p!("PHR" * nm, r.p_hamed_rao)
        registra!("MDT" * nm, r.mdt_decada; dig = dg)
        registra!("FatorHR" * nm, r.fator_hr; dig = 2); registra_p!("QTeste" * nm, r.q_teste)
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
    MESES = ("janeiro", "fevereiro", "março", "abril", "maio", "junho", "julho", "agosto", "setembro", "outubro", "novembro", "dezembro")
    registra!("MesMaisChuvoso", MESES[argmax(clim.prec_9120)])
    registra!("ChuvaMarNormal", clim.prec_9120[3]; dig = 0)

    # normais oficiais do INMET (scripts/extrai_normais_inmet.py) × climatologia calculada
    no = CSV.read(joinpath(RAIZ, "data", "raw", "inmet_normais", "normais_82983.csv"), DataFrame)
    of(p, v) = (g = no[(no.periodo .== p) .& (no.variavel .== v), :]; (mensal = g.valor[sortperm(g.mes)][2:13], anual = only(g.valor[g.mes .== 0])))
    p91, x91, n91 = of("1991-2020", "prec"), of("1991-2020", "tmax"), of("1991-2020", "tmin")
    registra!("NormalOficialChuva", p91.anual; dig = 0)
    registra!("NormalOficialChuvaSeisUm", of("1961-1990", "prec").anual; dig = 0)
    registra!("NormalOficialChuvaOitoUm", of("1981-2010", "prec").anual; dig = 0)
    dp, ip = findmax(abs.(clim.prec_9120 .- p91.mensal))
    registra!("DifNormalChuvaMensalMax", dp; dig = 1); registra!("MesDifNormalChuva", MESES[ip])
    registra!("DifNormalTempMax", maximum(abs.(vcat(clim.tmax_9120 .- x91.mensal, clim.tmin_9120 .- n91.mensal))); dig = 2)
    registra!("ChuvaJanNormal", clim.prec_9120[1]; dig = 0)
    registra!("ChuvaJanNormalOficial", p91.mensal[1]; dig = 1); registra!("ChuvaMarNormalOficial", p91.mensal[3]; dig = 1)
    registra!("ChuvaMarSeisUmOficial", of("1961-1990", "prec").mensal[3]; dig = 0)
    registra!("TmaxNormalOficialSeisUm", of("1961-1990", "tmax").anual; dig = 1)
    t7190 = dm[1971 .<= dm.ano .<= 1990, :]
    registra!("TmaxSetentaUmNoventa", mean(mean(skipmissing(g.tmax)) for g in groupby(t7190, :mes)); dig = 1)
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
    # sensibilidade ao controle de qualidade: Tmín anual da série sem CQ (com o segmento de 2006–08)
    dmb = mensal_longo(L.d_bruto; ref = REF)
    clb = Dict(m => mean(skipmissing(dmb[(dmb.mes .== m) .& (REF[1] .<= dmb.ano .<= REF[2]), :tmin])) for m in 1:12)
    tminb = map(anu.ano) do a
        g = dmb[dmb.ano .== a, :]
        an = [x - clb[m] for (x, m) in zip(g.tmin, g.mes) if !ismissing(x)]
        length(an) >= MIN_MESES_ANO ? mean(values(clb)) + mean(an) : missing
    end
    tb = tendencia_completa(anu.ano, tminb; nome = "Tmín sem CQ")
    registra_ic!("SenTminSemQC", tb.sen_decada, tb.sen_ic_inf, tb.sen_ic_sup; dig = 2)
    # extremos de temperatura (percentis por dia do calendário, base 1991–2020)
    p90 = percentil_calendario(d.data, d.tmax, 0.9); p10 = percentil_calendario(d.data, d.tmin, 0.1)
    ext = combine(groupby(d, :ano; sort = true)) do g
        a = g.ano[1]; nd = daysinyear(a)
        # completude dos índices de extremos (como no climdex): no máximo 15 dias ausentes no ano
        ausx = nd - count(!ismissing, g.tmax); ausn = nd - count(!ismissing, g.tmin)
        okx = !ismissing(anu.tmax[anu.ano .== a][1]) && ausx <= 15
        okn = !ismissing(anu.tmin[anu.ano .== a][1]) && ausn <= 15
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
    # Tmáx e Tmín nos MESMOS anos (os da DTR): na Tab. 7 cada série usa seus anos válidos
    cm = dropmissing(anu[:, [:ano, :tmax, :tmin]])
    for (vv, k) in ((:tmax, "Tmax"), (:tmin, "Tmin"))
        tt = analise_tendencia(cm.ano, cm[!, vv]; nome = k * " anos comuns")
        registra_ic!("Sen" * k * "Comum", tt.sen_decada, tt.sen_ic_inf, tt.sen_ic_sup; dig = 2)
    end
    registra!("NanosTempComum", nrow(cm))
    # sensibilidade à quebra A307 × 82983 (≈ 2009): tendência só até 2008
    for (vv, k) in ((:tmin, "Tmin"), (:tmax, "Tmax"))
        sub = anu[anu.ano .<= 2008, :]
        tt = tendencia_completa(sub.ano, sub[!, vv]; nome = k * " até 2008")
        registra_ic!("Sen" * k * "AteDoisMilOito", tt.sen_decada, tt.sen_ic_inf, tt.sen_ic_sup; dig = 2)
        registra_p!("PTFPW" * k * "AteDoisMilOito", tt.p_mk_tfpw)
    end

    # demanda evaporativa, fora da família BH: ETo de Hargreaves (a do SPEI) e evaporação Piche,
    # médias anuais via anomalias mensais (mesma regra das demais médias)
    linhas = readlines(ARQ_LONGO); i = findfirst(startswith("Data Medicao"), linhas)
    bp = CSV.read(IOBuffer(join(linhas[i:end], '\n')), DataFrame; delim = ';', header = 1, types = String,
                  missingstring = ["", "null"], silencewarnings = true)
    col = findfirst(n -> startswith(uppercase(n), "EVAPORACAO DO PICHE"), names(bp))
    ev = DataFrame(data = Date.(bp[!, 1]), evap = parse_num.(bp[!, col]))
    ev = ev[year.(ev.data) .>= INICIO_TEMP, :]; ev.ano = year.(ev.data); ev.mes = month.(ev.data)
    evm = combine(groupby(ev, [:ano, :mes]; sort = true)) do g
        nd = daysinmonth(Date(g.ano[1], g.mes[1])); f = ismissing.(g.evap)
        ok = count(f) + nd - nrow(g) <= MAX_DIAS_FALTA && max_consecutivo(f) <= MAX_DIAS_CONSEC
        (evap = ok ? mean(skipmissing(g.evap)) : missing,)
    end
    dme = leftjoin(select(dm, :ano, :mes, :tmax, :tmin), evm, on = [:ano, :mes])
    dme.eto = [eto_mensal(a, m, x, n, LAT_82983) for (a, m, x, n) in zip(dme.ano, dme.mes, dme.tmax, dme.tmin)]
    for (vv, k) in ((:eto, "ETo"), (:evap, "Piche"))
        cl = Dict(m => mean(skipmissing(dme[(dme.mes .== m) .& (REF[1] .<= dme.ano .<= REF[2]), vv])) for m in 1:12)
        anos = sort(unique(dme.ano))
        y = map(anos) do a
            g = dme[dme.ano .== a, :]
            an = [x - cl[m] for (x, m) in zip(g[!, vv], g.mes) if !ismissing(x)]
            length(an) >= MIN_MESES_ANO ? mean(values(cl)) + mean(an) : missing
        end
        tt = tendencia_completa(anos, y; nome = k)
        media = mean(skipmissing(y))
        registra_ic!("Sen" * k * "Pct", 100tt.sen_decada / media, 100tt.sen_ic_inf / media, 100tt.sen_ic_sup / media; dig = k == "ETo" ? 2 : 1)
        registra_p!("PTFPW" * k, tt.p_mk_tfpw); registra!("N" * k, tt.n)
        @printf("%s: n = %d, Sen = %.2f%% por década, p TFPW = %.3f\n", k, tt.n, 100tt.sen_decada / media, tt.p_mk_tfpw)
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
            if rot == "Tmin"     # segmento suspeito da Tmín (out/2006–dez/2008), excluído pelo CQ
                noseg(dt) = SEG_TMIN_SUSPEITO[1] <= dt <= SEG_TMIN_SUSPEITO[2]
                ds = [x - y for (dt, x, y) in zip(j.data, j.tmin, j.tn) if !ismissing(x) && !ismissing(y) && noseg(dt)]
                df_ = [x - y for (dt, x, y) in zip(j.data, j.tmin, j.tn) if !ismissing(x) && !ismissing(y) && !noseg(dt) && year(dt) >= INICIO_TEMP]
                registra!("DifSegTmin" * nm, mean(ds); dig = 1); registra!("DifForaTmin" * nm, mean(df_); dig = 1)
            end
            @printf("%s × %s: diferença média 1963–65 = %.2f °C; ≥ %d = %.2f °C\n", rot, nm, mean(ini), INICIO_TEMP, mean(pos))
        end
    end
    a307 = CSV.read(joinpath(DIR_PROC, "diario.csv"), DataFrame)
    a307.tmin = [p ? missing : x for (x, p) in zip(a307.tmin, a307.tmin_pre)]
    ja = innerjoin(select(a307, :data, :tmin => :ta), select(bruto, :data, :tmin), on = :data)
    noseg2(dt) = SEG_TMIN_SUSPEITO[1] <= dt <= SEG_TMIN_SUSPEITO[2]
    dsa = [a - b for (dt, a, b) in zip(ja.data, ja.ta, ja.tmin) if !ismissing(a) && !ismissing(b) && noseg2(dt)]
    dfa = [a - b for (dt, a, b) in zip(ja.data, ja.ta, ja.tmin) if !ismissing(a) && !ismissing(b) && !noseg2(dt)]
    registra!("DifSegTminAtres", mean(dsa); dig = 1); registra!("DifForaTminAtres", mean(dfa); dig = 1)
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
    # Tmín: as quebras frente às duas vizinhas caem na mesma direção e perto da mesma data (por volta
    # de 1990), o que aponta para a 82983. Mesma sensibilidade da Tmáx, com o salto frente à vizinha
    # em que a quebra é mais significativa.
    sn_ = dropmissing(anu[:, [:ano, :tmin]])
    hn = sort(h[h.variavel .== "Tmín", :], :q)
    if nrow(hn) >= 1
        a0n = hn.quebra_snht[1] + 1; saln = hn.salto[1]
        registra!("VizTminAjuste", hn.vizinha[1] == "Cabrobó" ? "Cabrobó" : hn.vizinha[1])
        tn_ = tendencia_completa(sn_.ano, sn_.tmin .- saln .* (sn_.ano .>= a0n); nome = "Tmín ajustada")
        registra_ic!("SenTminAjustada", tn_.sen_decada, tn_.sen_ic_inf, tn_.sen_ic_sup; dig = 2)
        registra_p!("PTFPWTminAjustada", tn_.p_mk_tfpw)
        ddn = DataFrame(y = sn_.tmin, t = (sn_.ano .- mean(sn_.ano)) ./ 10, d = Float64.(sn_.ano .>= a0n))
        mn_ = lm(@formula(y ~ t + d), ddn); icn = confint(mn_)
        registra_ic!("BTminDegrau", coef(mn_)[2], icn[2, 1], icn[2, 2]; dig = 2)
        registra!("DegrauTmin", coef(mn_)[3]; dig = 2)
        # mesma regra da Tmáx (salto frente a Cabrobó), para a correção não depender da escolha
        hcn = h[(h.variavel .== "Tmín") .& (h.vizinha .== "Cabrobó"), :]
        if nrow(hcn) == 1
            a0c = hcn.quebra_snht[1] + 1
            tc_ = tendencia_completa(sn_.ano, sn_.tmin .- hcn.salto[1] .* (sn_.ano .>= a0c); nome = "Tmín ajustada (Cabrobó)")
            registra!("SenTminAjustadaCabrobo", tc_.sen_decada; dig = 2)
        end
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
    qs = Dict(v => quantilica_tempo(L.d.data, L.d[!, v]) for v in (:tmax, :tmin))
    for (v, k) in ((:tmax, "Tmax"), (:tmin, "Tmin"))
        dq = metadata(qs[v], "dif")                       # quantil 0,9 menos quantil 0,1
        registra_ic!("QRDif" * k, dq.est, dq.inf, dq.sup; dig = 2); registra_p!("PQRDif" * k, dq.p)
    end
    q = vcat([insertcols!(copy(qs[v]), 1, :variavel => String(v)) for v in (:tmax, :tmin)]...)
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
    # sensibilidade: chuva sem tendência pela normal não estacionária (DP constante)
    nn = ns_candidata(Normal, Float64.(v.prec), tc.(v.ano_hidro))
    df.z_normal = [quantile(Normal(), clamp(nn.F(p, tc(a)), 1e-6, 1 - 1e-6)) for (p, a) in zip(v.prec, v.ano_hidro)]
    for (ind, k) in ((:oni, "Oni"), (:grad, "Grad"))
        s = dropmissing(df[:, [ind, :z_normal]]); r = cor(s[!, ind], s.z_normal)
        lo, hi = ic_pearson(r, nrow(s))
        registra_ic!("R" * k * "DetNormal", r, lo, hi; dig = 2)
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
    for (crit, k) in (("bruto", ""), ("sem tendência", "Det"))
        a = only(eachrow(tcd[(tcd.fase .== "El Niño") .& (tcd.criterio .== crit), :]))
        b = only(eachrow(tcd[(tcd.fase .== "neutro") .& (tcd.criterio .== crit), :]))
        ft_ = FisherExactTest(a.secos, a.n - a.secos, b.secos, b.n - b.secos)
        registra_p!("PFisherNinoNeutro" * k, pvalue(ft_))                       # bilateral central
        registra_p!("PFisherNinoNeutroMinlike" * k, pvalue(ft_; method = :minlike))  # soma das tabelas menos prováveis
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
    registra!("RhatMaxTodos", max(N.rhat, maximum(B.res.rhat_logk), maximum(B.res.rhat_logmu)); dig = 3)
    VF = bloco_validacao(L, E)
    RA = bloco_anos_recentes(L, E, N)
    AN = bloco_ana(L, E)
    AR = bloco_ana_retorno(L, E, AN)
    RG = bloco_regional(L, E, AN)
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
    return (L = L, H = H, V = V, E = E, R = R, B = B, N = N, VF = VF, RA = RA, AN = AN, AR = AR, RG = RG, I = I, M = M, S = S, Z = Z, X = X, HT = HT, QR = QR, ST = ST, T = T)
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
