# =============================================================================
# 08_modelagem.jl — modelos de tendência com diagnóstico e comparação de períodos
#
# Modelo base:  Y_t = β₀ + β₁(t − t̄) + ε_t   (MQO; t centrado ⇒ β₀ = média no ano central)
# Pressupostos e como são verificados:
#   normalidade dos resíduos  → Shapiro–Wilk  (H₀: resíduos normais)
#   homocedasticidade         → Breusch–Pagan (H₀: Var(ε) constante; LM = n·R²_aux ~ χ²₁)
#   independência             → r₁ dos resíduos, Durbin–Watson e Ljung–Box
#                               (H₀: sem autocorrelação até a defasagem h)
# Alternativa se houver autocorrelação: MQG com erros AR(1),
#   Corr(ε_s, ε_t) = ρ^{|s−t|}, estimado por máxima verossimilhança perfilada em ρ.
#   Essa forma lida exatamente com anos ausentes (a correlação decai com a distância
#   real no tempo), algo que Cochrane–Orcutt/Prais–Winsten não fazem.
# Para chuva (positiva, variância crescente com a média) compara-se também um
#   GLM Gamma com ligação log: E[Y_t] = exp(β₀ + β₁t) ⇒ variação de 100(e^{10β₁} − 1) %/década.
# =============================================================================

function breusch_pagan(m)
    e2 = residuals(m) .^ 2; X = modelmatrix(m)
    aux = X \ e2; ŷ = X * aux
    R2 = 1 - sum(abs2, e2 .- ŷ) / sum(abs2, e2 .- mean(e2))
    LM = length(e2) * R2
    return LM, ccdf(Chisq(size(X, 2) - 1), LM)
end

durbin_watson(e) = sum(abs2, diff(e)) / sum(abs2, e)

"MQG com erros AR(1) em tempo contínuo (anos com lacunas), ρ por MV perfilada."
function gls_ar1(t, y)
    n = length(y); X = [ones(n) t]
    D = abs.(t .- t')
    melhor = (ll = -Inf,)
    for ρ in -0.9:0.005:0.95
        Σ = ρ .^ D
        C = cholesky(Symmetric(Σ); check = false); issuccess(C) || continue
        Xs = C.L \ X; ys = C.L \ y
        β = Xs \ ys; r = ys - Xs * β
        σ2 = sum(abs2, r) / n
        ll = -n / 2 * (log(2π * σ2) + 1) - logdet(C) / 2
        ll > melhor.ll && (melhor = (ll = ll, ρ = ρ, β = β, σ2 = σ2, V = inv(Xs' * Xs)))
    end
    ep = sqrt.(diag(melhor.V) .* melhor.σ2 .* n ./ (n - 2))
    return (ρ = melhor.ρ, β = melhor.β, ep = ep, ll = melhor.ll, aic = 2 * 4 - 2melhor.ll)
end

"""
    modelo_tendencia(t, y; nome, escala=10)

Ajusta MQO e (sempre) o MQG-AR(1) para comparação. Coeficientes por década.
"""
function modelo_tendencia(t, y; nome = "", escala = 10)
    ok = .!ismissing.(y)
    tc = mean(Float64.(t[ok]))
    d = DataFrame(t = Float64.(t[ok]) .- tc, y = Float64.(y[ok]), ano = Float64.(t[ok]))
    m = lm(@formula(y ~ t), d); n = nrow(d)
    e = residuals(m)
    tq = quantile(TDist(n - 2), 0.975)
    b, se = coef(m), stderror(m)
    sw = ShapiroWilkTest(e); bp = breusch_pagan(m)
    h = min(5, n ÷ 4); lb = LjungBoxTest(e, h)
    r1 = acf_irregular(d.ano, e, 1)[1]
    ar = gls_ar1(d.t, d.y)
    ll_ols = loglikelihood(m)
    LR = 2(ar.ll - ll_ols)
    resumo = (variavel = nome, n = n, beta0 = b[1], ep_beta0 = se[1],
        beta1_decada = escala * b[2], ep_decada = escala * se[2],
        ic_inf = escala * (b[2] - tq * se[2]), ic_sup = escala * (b[2] + tq * se[2]),
        p_beta1 = coeftable(m).cols[4][2], R2 = r2(m), R2aj = adjr2(m),
        AIC = aic(m), BIC = bic(m),
        p_shapiro = pvalue(sw), BP = bp[1], p_bp = bp[2],
        DW = durbin_watson(e), r1_resid = r1, p_ljungbox = pvalue(lb), lags_lb = h,
        rho_ar1 = ar.ρ, beta1_ar1_decada = escala * ar.β[2], ep_ar1_decada = escala * ar.ep[2],
        p_beta1_ar1 = 2ccdf(TDist(n - 2), abs(ar.β[2] / ar.ep[2])),
        AIC_ar1 = ar.aic, p_LR_rho = ccdf(Chisq(1), max(LR, 0)))
    return (resumo = resumo, m = m, dados = d, tc = tc)
end

"Faixa de confiança da reta MQO para a figura."
function faixa_confianca(mt, tgrid)
    nd = DataFrame(t = tgrid .- mt.tc)
    p = predict(mt.m, nd; interval = :confidence, level = 0.95)
    return Float64.(p.prediction), Float64.(p.lower), Float64.(p.upper)
end

"GLM Gamma (ligação log) para chuva anual: variação percentual por década."
function tendencia_gamma(t, y)
    ok = .!ismissing.(y)
    d = DataFrame(t = Float64.(t[ok]) .- mean(t[ok]), y = Float64.(y[ok]))
    m = glm(@formula(y ~ t), d, Gamma(), LogLink())
    β, ep = coef(m)[2], stderror(m)[2]; n = nrow(d); tq = quantile(TDist(n - 2), 0.975)
    (pct_decada = 100(exp(10β) - 1), ic_inf = 100(exp(10(β - tq * ep)) - 1),
     ic_sup = 100(exp(10(β + tq * ep)) - 1), p = 2ccdf(TDist(n - 2), abs(β / ep)), AIC = aic(m), m = m)
end

"""
Comparação de dois grupos (ex.: período antigo × recente).
• Welch: diferença de médias com IC t, sem supor variâncias iguais.
  H₀: μ_B = μ_A.
• Mann–Whitney: H₀: P(B > A) = 1/2. Efeito: deslocamento de Hodges–Lehmann
  (mediana das diferenças Bᵢ − Aⱼ) e δ de Cliff = P(B > A) − P(B < A).
• g de Hedges: (x̄_B − x̄_A)/s_pooled × J, J = 1 − 3/(4(n_A+n_B) − 9).
• IC bootstrap (percentil, reamostragem estratificada) para diferença de médias e HL.
Ambos supõem observações independentes — com autocorrelação positiva, os p-valores
são otimistas (ver r₁ nos modelos de tendência).
"""
function compara_grupos(a, b; nome = "", rotulo_a = "A", rotulo_b = "B", rng = MersenneTwister(SEMENTE))
    a = Float64.(collect(skipmissing(a))); b = Float64.(collect(skipmissing(b)))
    na, nb = length(a), length(b)
    w = UnequalVarianceTTest(b, a); ci = confint(w)
    sp = sqrt(((na - 1) * var(a) + (nb - 1) * var(b)) / (na + nb - 2))
    J = 1 - 3 / (4(na + nb) - 9); g = J * (mean(b) - mean(a)) / sp
    se_g = sqrt((na + nb) / (na * nb) + g^2 / (2(na + nb)))
    hl(x, y) = median(vec([yi - xi for xi in x, yi in y]))
    cliff = (sum(yi > xi for xi in a, yi in b) - sum(yi < xi for xi in a, yi in b)) / (na * nb)
    bd = [(xa = a[rand(rng, 1:na, na)]; xb = b[rand(rng, 1:nb, nb)]; (mean(xb) - mean(xa), hl(xa, xb))) for _ in 1:NBOOT]
    (variavel = nome, grupo_A = rotulo_a, n_A = na, media_A = mean(a), grupo_B = rotulo_b, n_B = nb, media_B = mean(b),
     dif_medias = mean(b) - mean(a), ic_welch_inf = ci[1], ic_welch_sup = ci[2], p_welch = pvalue(w),
     ic_boot_inf = quantile(first.(bd), 0.025), ic_boot_sup = quantile(first.(bd), 0.975),
     hedges_g = g, g_ic_inf = g - 1.96se_g, g_ic_sup = g + 1.96se_g,
     hodges_lehmann = hl(a, b), hl_ic_inf = quantile(last.(bd), 0.025), hl_ic_sup = quantile(last.(bd), 0.975),
     cliff_delta = cliff, p_mann_whitney = pvalue(MannWhitneyUTest(b, a)))
end
