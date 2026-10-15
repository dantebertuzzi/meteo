# =============================================================================
# 04_temporal.jl — tendência e pontos de mudança em séries anuais curtas
#
# Por que estes métodos:
# • Mann-Kendall (MK): teste de postos para tendência monotônica. Não exige
#   normalidade nem linearidade — adequado à chuva anual (assimétrica) e robusto
#   a um ou dois anos extremos. Supõe observações independentes.
# • Theil-Sen: estimador do declive = mediana dos declives entre pares. É o
#   complemento natural do MK (mesma base de postos) e tem ponto de ruptura ≈29 %.
# • TFPW (Yue et al., 2002): se houver autocorrelação lag-1, o MK tem tamanho
#   inflado; o pré-branqueamento livre de tendência remove a AR(1) preservando o declive.
# • Pettitt: detecta um único ponto de mudança na posição (teste de postos).
#   Serve para cruzar possíveis quebras com trocas de instrumento.
# • Médias móveis centradas de 5 anos: suavização descritiva (não inferencial).
# =============================================================================

"""
    mann_kendall(t, y)

S = Σ_{i<j} sgn(y_j − y_i)
Var(S) = [n(n−1)(2n+5) − Σ_k t_k(t_k−1)(2t_k+5)] / 18   (correção de empates)
Z = (S − sgn(S))/√Var(S)   (correção de continuidade)
H₀: ausência de tendência monotônica (y i.i.d.); H₁: tendência (bilateral).
τ de Kendall = S / [n(n−1)/2].
"""
function mann_kendall(t, y)
    o = sortperm(t); y = y[o]; n = length(y)
    S = 0
    for i in 1:n-1, j in i+1:n
        S += sign(y[j] - y[i])
    end
    empates = values(countmap(y))
    varS = (n * (n - 1) * (2n + 5) - sum(k * (k - 1) * (2k + 5) for k in empates)) / 18
    Z = S == 0 ? 0.0 : (S - sign(S)) / sqrt(varS)
    p = 2 * ccdf(Normal(), abs(Z))
    return (n = n, S = S, varS = varS, Z = Z, p = p, tau = S / (n * (n - 1) / 2))
end

"""
    theil_sen(t, y; nivel)

β̂ = mediana{(y_j − y_i)/(t_j − t_i)}, i < j. IC não paramétrico (Gilbert, 1987):
C = z_{1−α/2}√Var(S); posições M₁ = (N − C)/2 e M₂ = (N + C)/2 + 1 nos N declives ordenados.
Intercepto: mediana(y − β̂t).
"""
function theil_sen(t, y; nivel = 1 - ALFA)
    n = length(y)
    d = sort!([(y[j] - y[i]) / (t[j] - t[i]) for i in 1:n-1 for j in i+1:n if t[j] != t[i]])
    N = length(d); b = median(d)
    C = quantile(Normal(), 1 - (1 - nivel) / 2) * sqrt(mann_kendall(t, y).varS)
    m1 = clamp(round(Int, (N - C) / 2), 1, N); m2 = clamp(round(Int, (N + C) / 2) + 1, 1, N)
    return (declive = b, ic_inf = d[m1], ic_sup = d[m2], intercepto = median(y .- b .* t))
end

"""
MK com pré-branqueamento livre de tendência (TFPW):
1) remove a tendência de Sen; 2) estima r₁ (só pares consecutivos no tempo);
3) yʷₜ = xₜ − r₁xₜ₋₁ + β̂t; 4) aplica MK em yʷ.
"""
function mk_tfpw(t, y)
    ts = theil_sen(t, y)
    x = y .- ts.declive .* t
    r1 = acf_irregular(t, x, 1)[1]
    pos = Dict(round(Int, ti) => xi for (ti, xi) in zip(t, x))
    tt = Float64[]; yw = Float64[]
    for (ti, xi) in zip(t, x)
        k = round(Int, ti) - 1
        haskey(pos, k) || continue
        push!(tt, ti); push!(yw, xi - r1 * pos[k] + ts.declive * ti)
    end
    return (r1 = r1, mk = mann_kendall(tt, yw))
end

"""
Teste de Pettitt: U_k = Σ_{i≤k}Σ_{j>k} sgn(y_j − y_i); K = max|U_k|;
p ≈ 2·exp(−6K²/(n³ + n²)) (aproximação válida para p ≲ 0,5).
H₀: sem mudança de posição; H₁: uma mudança em algum ponto desconhecido.
"""
function pettitt(t, y)
    o = sortperm(t); t = t[o]; y = y[o]; n = length(y)
    U = [sum(sign(y[j] - y[i]) for i in 1:k for j in k+1:n) for k in 1:n-1]
    K, k = findmax(abs.(U))
    p = min(1.0, 2exp(-6K^2 / (n^3 + n^2)))
    return (K = K, p = p, ultimo_antes = t[k], media_antes = mean(y[1:k]), media_depois = mean(y[k+1:end]))
end

"Média móvel centrada de janela `w` anos, exigindo `minobs` valores na janela."
function media_movel(t, y; w = 5, minobs = 3)
    h = w ÷ 2
    map(t) do ti
        v = [yi for (tj, yi) in zip(t, y) if abs(tj - ti) <= h && !ismissing(yi)]
        length(v) >= minobs ? mean(v) : missing
    end
end

"Bateria de tendência para uma série anual (valores ausentes descartados)."
function analise_tendencia(t, y; nome = "", escala = 10)
    ok = .!ismissing.(y); t = Float64.(t[ok]); y = Float64.(y[ok])
    mk = mann_kendall(t, y); ts = theil_sen(t, y); tf = mk_tfpw(t, y); pt = pettitt(t, y)
    (variavel = nome, n = mk.n, periodo = "$(Int(minimum(t)))–$(Int(maximum(t)))",
     tau = mk.tau, Z_mk = mk.Z, p_mk = mk.p, r1 = tf.r1, p_mk_tfpw = tf.mk.p,
     sen_decada = escala * ts.declive, sen_ic_inf = escala * ts.ic_inf, sen_ic_sup = escala * ts.ic_sup,
     pettitt_K = pt.K, p_pettitt = pt.p, mudanca_apos = Int(pt.ultimo_antes),
     media_antes = pt.media_antes, media_depois = pt.media_depois)
end
