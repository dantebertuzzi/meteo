# =============================================================================
# 05_distribuicoes.jl — ajuste e comparação de distribuições (máxima verossimilhança)
#
# Candidatas para a chuva anual (positiva, possivelmente assimétrica):
#   Normal(μ,σ), Lognormal(μ,σ), Gamma(k,θ), Weibull(k,λ).
# A Exponencial NÃO é considerada: tem moda em zero e CV = 1 fixo, incompatível
# com totais anuais (que ficam longe de zero e têm CV ≈ 0,3–0,5).
#
# Critérios: ℓ̂ (log-verossimilhança), AIC = 2k − 2ℓ̂, BIC = k·ln n − 2ℓ̂,
# pesos de Akaike wᵢ = exp(−Δᵢ/2)/Σexp(−Δⱼ/2) e o estatístico de
# Anderson–Darling A² = −n − (1/n)Σ(2i−1)[ln F(x₍ᵢ₎) + ln(1 − F(x₍ₙ₊₁₋ᵢ₎))],
# sensível às caudas. Como os parâmetros são estimados dos mesmos dados, o p-valor
# tabelado de A² é inválido; usa-se bootstrap paramétrico (reajustando a cada réplica).
# =============================================================================

const CANDIDATAS = [Normal, LogNormal, Gamma, Weibull]

"""
Ajuste por máxima verossimilhança. Para a log-normal usa-se a forma fechada
(μ̂, σ̂) = (média, desvio-padrão com divisor n) de log x: o `fit_mle(LogNormal, x)`
do Distributions.jl usa o divisor n − 1 (verificado contra o SciPy), o que não é o EMV
e distorce levemente o AIC.
"""
ajusta_mle(D, x) = D === LogNormal ? (l = log.(x); LogNormal(mean(l), std(l; corrected = false))) : fit_mle(D, x)

function anderson_darling(d, x)
    x = sort(x); n = length(x)
    F = clamp.(cdf.(d, x), 1e-12, 1 - 1e-12)
    -n - sum((2i - 1) * (log(F[i]) + log(1 - F[n+1-i])) for i in 1:n) / n
end

function p_ad_bootstrap(D, x; B = 1999, rng = MersenneTwister(SEMENTE))
    dhat = ajusta_mle(D, x); A = anderson_darling(dhat, x); n = length(x)
    maiores = 0
    for _ in 1:B
        xb = rand(rng, dhat, n)
        db = try ajusta_mle(D, xb) catch; continue end
        maiores += anderson_darling(db, xb) >= A
    end
    return A, (maiores + 1) / (B + 1)
end

"Tabela comparativa de ajuste para o vetor `x`."
function compara_distribuicoes(x; candidatas = CANDIDATAS)
    x = collect(Float64, skipmissing(x)); n = length(x)
    linhas = map(candidatas) do D
        d = ajusta_mle(D, x); k = length(params(d)); ll = loglikelihood(d, x)
        A, pA = p_ad_bootstrap(D, x)
        (distribuicao = string(nameof(D)), parametros = join(map(p -> @sprintf("%.4g", p), params(d)), "; "),
         k = k, loglik = ll, AIC = 2k - 2ll, BIC = k * log(n) - 2ll, AD = A, p_AD_boot = pA, modelo = d)
    end
    t = DataFrame(linhas)
    t.dAIC = t.AIC .- minimum(t.AIC)
    t.peso_akaike = exp.(-t.dAIC ./ 2) ./ sum(exp.(-t.dAIC ./ 2))
    sort!(t, :AIC)
end

"Testes de normalidade (Shapiro–Wilk) para várias séries — p-valores ajustados por BH."
function normalidade(df, vars)
    t = DataFrame(map(vars) do v
        x = collect(Float64, skipmissing(df[!, v]))
        sw = ShapiroWilkTest(x)
        (variavel = String(v), n = length(x), W = sw.W, p = pvalue(sw))
    end)
    t.p_bh = bh_ajuste(t.p)
    t
end
