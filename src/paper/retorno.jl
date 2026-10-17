# =============================================================================
# retorno.jl — período de retorno de anos secos (cauda inferior da chuva anual)
#
# Para um evento de referência x₀ (total anual), a probabilidade de não
# excedência é p = F(x₀) = P(X ≤ x₀) e o período de retorno médio é T = 1/p
# (anos). Supondo anos independentes e estacionariedade, o número de anos até o
# próximo evento é geométrico, e
#     P(≥ 1 evento em n anos) = 1 − (1 − p)ⁿ.
# A incerteza de p é propagada para T e para essa probabilidade de três formas:
#   1) bootstrap paramétrico (reamostra da Gamma ajustada e reajusta);
#   2) verossimilhança perfilada de p (IC de Wilks, χ²₁);
#   3) posterior bayesiana (Turing, NUTS) com prioris fracamente informativas.
# =============================================================================

# ----------------------------------------------------------------------------- empírico
"Posição de plotagem de Weibull: P(X ≤ x₍ᵢ₎) ≈ i/(n + 1) — não viesada para a probabilidade."
pos_weibull(n) = (1:n) ./ (n + 1)

# ----------------------------------------------------------------------------- bootstrap
"""
    bootstrap_p(D, x, x0; B)

Bootstrap paramétrico de p = F(x₀) sob a família `D`. Devolve as B réplicas de p.
"""
function bootstrap_p(D, x, x0; B = NBOOT, rng = MersenneTwister(SEMENTE))
    d = fit_mle(D, x); n = length(x)
    ps = Float64[]
    for _ in 1:B
        xb = rand(rng, d, n)
        db = try fit_mle(D, xb) catch; continue end
        push!(ps, cdf(db, x0))
    end
    return ps
end

# ----------------------------------------------------------------------------- perfilada
"""
    perfil_p_gamma(x, x0; nivel)

Reparametriza a Gamma por (k, p): para k fixo, o θ que faz F(x₀) = p é
θ = x₀ / Q_{Gamma(k,1)}(p). O perfil ℓₚ(p) = max_k ℓ(k, θ(k, p)) e o IC é
{p : 2[ℓ̂ − ℓₚ(p)] ≤ χ²₁(nível)}. Os limites são achados por bisseção em logit(p).
"""
function perfil_p_gamma(x, x0; nivel = 0.95)
    ll(k, p) = (θ = x0 / quantile(Gamma(k, 1.0), p); loglikelihood(Gamma(k, θ), x))
    function lp(p)
        r = optimize(lk -> -ll(exp(lk), p), log(0.05), log(200.0))
        return -Optim.minimum(r)
    end
    dhat = fit_mle(Gamma, x); phat = cdf(dhat, x0); lmax = loglikelihood(dhat, x)
    corte = quantile(Chisq(1), nivel) / 2
    g(p) = lmax - lp(p) - corte
    logit(p) = log(p / (1 - p)); expit(u) = 1 / (1 + exp(-u))
    function raiz(a, b)                          # bisseção em logit(p)
        fa = g(expit(a))
        for _ in 1:80
            m = (a + b) / 2; fm = g(expit(m))
            sign(fm) == sign(fa) ? (a = m; fa = fm) : (b = m)
        end
        expit((a + b) / 2)
    end
    u = logit(phat)
    inf = g(expit(u - 15)) > 0 ? raiz(u - 15, u) : 0.0
    sup = g(expit(u + 10)) > 0 ? raiz(u, u + 10) : 1.0
    grade = expit.(range(u - 4, u + 3, length = 120))
    return (p = phat, inf = inf, sup = sup, perfil = DataFrame(p = grade, desvio = [2(lmax - lp(p)) for p in grade]))
end

# ----------------------------------------------------------------------------- resumo
"Resumo de p, T e P(≥1 em n anos) a partir de réplicas (bootstrap ou posterior)."
function resumo_replicas(ps; horizontes = (5, 10, 20, 30), nivel = 0.95)
    a = (1 - nivel) / 2
    q(v) = quantile(v, [a, 0.5, 1 - a])
    out = (p = q(ps), T = q(1 ./ ps),
           prob_n = Dict(n => q(1 .- (1 .- ps) .^ n) for n in horizontes))
    return out
end

prob_ao_menos_um(p, n) = 1 - (1 - p)^n

# ----------------------------------------------------------------------------- bayes
"""
Modelo Gamma parametrizado por média μ e forma k (CV = 1/√k):
    log k ~ Normal(log 4, 1)     ⇒ CV entre ~0,2 e ~1,3 com 95 % de prob.
    log μ ~ Normal(log 500, 1)   ⇒ média entre ~70 e ~3 500 mm.
São fracamente informativas: cobrem todo o semiárido e excluem só valores
fisicamente absurdos. A sensibilidade é testada com desvios-padrão 2× maiores.
"""
@model function modelo_gamma(y; sd_k = 1.0, sd_mu = 1.0)
    logk ~ Normal(log(4.0), sd_k)
    logmu ~ Normal(log(500.0), sd_mu)
    k = exp(logk); μ = exp(logmu)
    for i in eachindex(y)
        y[i] ~ Gamma(k, μ / k)
    end
end

"Amostra a posterior e devolve as réplicas de p = F(x₀) (uma por draw)."
function posterior_p(y, x0; n = 4000, cadeias = 4, sd_k = 1.0, sd_mu = 1.0, rng = MersenneTwister(SEMENTE))
    ch = sample(rng, modelo_gamma(y; sd_k, sd_mu), NUTS(0.9), MCMCThreads(), n ÷ cadeias, cadeias; progress = false)
    lk = vec(Array(ch[:logk])); lm = vec(Array(ch[:logmu]))
    ps = [cdf(Gamma(exp(a), exp(b) / exp(a)), x0) for (a, b) in zip(lk, lm)]
    return ps, ch
end

"R̂ e ESS básicos (Gelman–Rubin com divisão das cadeias), sem depender da API do objeto de cadeias."
function diagnostico_cadeias(ch, nome)
    M = Array(ch[nome])                       # iterações × cadeias
    M = ndims(M) == 1 ? reshape(M, :, 1) : M
    n, m = size(M); h = n ÷ 2
    sub = hcat([M[1:h, j] for j in 1:m]..., [M[h+1:2h, j] for j in 1:m]...)
    médias = vec(mean(sub, dims = 1)); W = mean(var(sub, dims = 1)); B = h * var(médias)
    rhat = sqrt(((h - 1) / h * W + B / h) / W)
    return rhat
end

# ----------------------------------------------------------------------------- não estacionário
"""
    gamma_nao_estacionaria(y, x; x0)

Gamma com forma k constante e média log-linear num covariável:
    log μᵢ = a + b·xᵢ   (x = tempo centrado em décadas, ou ONI).
Ajuste por MV (Optim, BFGS). Compara-se com o modelo estacionário por razão de
verossimilhança (χ²₁) e AIC.
"""
function gamma_nao_estacionaria(y, x)
    nll(θ) = (k = exp(θ[1]); -sum(logpdf(Gamma(k, exp(θ[2] + θ[3] * xi) / k), yi) for (yi, xi) in zip(y, x)))
    d0 = fit_mle(Gamma, y)
    θ0 = [log(shape(d0)), log(mean(y)), 0.0]
    r = optimize(nll, θ0, BFGS())
    θ = Optim.minimizer(r); ll1 = -Optim.minimum(r); ll0 = loglikelihood(d0, y)
    H = hessiana(nll, θ; h = 1e-5)
    ep = sqrt.(max.(diag(inv(H)), 0.0))
    LR = 2(ll1 - ll0)
    return (k = exp(θ[1]), a = θ[2], b = θ[3], ep_b = ep[3], ll = ll1, ll0 = ll0,
            AIC = 6 - 2ll1, AIC0 = 4 - 2ll0, LR = LR, p_LR = ccdf(Chisq(1), LR),
            p_em = (x0, xi) -> cdf(Gamma(exp(θ[1]), exp(θ[2] + θ[3] * xi) / exp(θ[1])), x0))
end

# ----------------------------------------------------------------------------- não estacionário: extensões
"Gamma com log μᵢ = a + b_fixo·xᵢ (b imposto, p.ex. a taxa regional); estima a e k por MV."
function gamma_inclinacao_fixa(y, x, b)
    nll(θ) = (k = exp(θ[1]); -sum(logpdf(Gamma(k, exp(θ[2] + b * xi) / k), yi) for (yi, xi) in zip(y, x)))
    r = optimize(nll, [log(shape(fit_mle(Gamma, y))), log(mean(y))], NelderMead())
    θ = Optim.minimizer(r)
    return (k = exp(θ[1]), a = θ[2], b = b, ll = -Optim.minimum(r),
            p_em = (x0, xi) -> cdf(Gamma(exp(θ[1]), exp(θ[2] + b * xi) / exp(θ[1])), x0))
end

# ----------------------------------------------------------------------------- não estacionário: média entre candidatas
"""
    ns_candidata(D, y, x; b, forma_var)

Ajusta por MV a família `D` com parâmetro de posição log-linear no covariável x:
    Normal e Gamma: log(média) = a + b·x;   Weibull: log(escala) = a + b·x;
    LogNormal: média de log y = a + b·x.
Com `forma_var = false`, a forma (ou o desvio-padrão) é constante: a normal mantém o
DP e as demais mantêm o coeficiente de variação, de modo que todas têm a mesma
tendência multiplicativa, exp(b) − 1 por unidade de x. Com `forma_var = true`, o log
da forma também varia: log s = s₀ + c·x, o que libera a relação entre a
variabilidade e a média. Com `b` dado, a inclinação da posição é fixada (p.ex. na
taxa regional).
"""
function ns_candidata(D, y, x; b = nothing, forma_var = false)
    livre_b = isnothing(b)
    function dist(θ, xi)
        η = θ[2] + (livre_b ? θ[3] : b) * xi
        s = exp(θ[1] + (forma_var ? θ[end] * xi : 0.0))
        D === Normal    ? Normal(exp(η), s) :
        D === LogNormal ? LogNormal(η, s) :
        D === Gamma     ? Gamma(s, exp(η) / s) :
        D === Weibull   ? Weibull(s, exp(η)) : error("família não suportada: $D")
    end
    nll(θ) = (v = -sum(logpdf(dist(θ, xi), yi) for (yi, xi) in zip(y, x)); isfinite(v) ? v : 1e12)
    if forma_var
        θ0 = [ns_candidata(D, y, x; b).θ; 0.0]          # parte do ajuste com forma fixa
    else
        d0 = ajusta_mle(D, y)
        s0 = D === Normal ? std(y) : D === LogNormal ? std(log.(y)) : shape(d0)
        a0 = D === LogNormal ? mean(log.(y)) : D === Weibull ? log(scale(d0)) : log(mean(y))
        θ0 = livre_b ? [log(s0), a0, 0.0] : [log(s0), a0]
    end
    r = optimize(nll, θ0, NelderMead(), Optim.Options(iterations = 20_000, g_tol = 1e-12))
    r = optimize(nll, Optim.minimizer(r), BFGS())         # refina a partir do simplex
    θ = Optim.minimizer(r)
    return (D = D, forma_var = forma_var, θ = θ, ll = -Optim.minimum(r), npar = length(θ),
            b = livre_b ? θ[3] : b, c = forma_var ? θ[end] : 0.0,
            F = (x0, xi) -> cdf(dist(θ, xi), x0))
end

"""
    ns_media_akaike(y, x; b, formas)

Ajusta as candidatas não estacionárias (cada família com forma fixa e, se pedido,
com forma variável) e devolve a tabela (AIC e pesos de Akaike) e `p_em(x0, xi)`, a
média de F(x0) no ponto xi ponderada pelos pesos — o mesmo esquema do caso
estacionário, para que os dois sejam comparáveis.
"""
function ns_media_akaike(y, x; b = nothing, candidatas = CANDIDATAS, formas = (false, true))
    fits = [ns_candidata(D, y, x; b, forma_var = fv) for fv in formas for D in candidatas]
    aic = [2f.npar - 2f.ll for f in fits]
    w = exp.(-(aic .- minimum(aic)) ./ 2); w ./= sum(w)
    p_em(x0, xi) = sum(w[j] * fits[j].F(x0, xi) for j in eachindex(fits))
    tab = DataFrame(distribuicao = [string(nameof(f.D)) for f in fits],
                    forma = [f.forma_var ? "variável" : "fixa" for f in fits],
                    efeito_pct = [100(exp(f.b) - 1) for f in fits],
                    efeito_forma_pct = [100(exp(f.c) - 1) for f in fits], loglik = [f.ll for f in fits],
                    AIC = aic, dAIC = aic .- minimum(aic), peso_akaike = w)
    return (fits = fits, tab = tab, w = w, p_em = p_em)
end

"""
Bootstrap não paramétrico dos anos (pares (x, y)) da média por pesos de Akaike:
as candidatas e os pesos são reajustados a cada réplica. Devolve uma matriz
réplicas × pontos, com p = F(x₀) em cada xi de `xs`.
"""
function ns_media_akaike_boot(y, x, x0, xs; b = nothing, formas = (false, true), B = 2000,
                              rng = MersenneTwister(SEMENTE))
    n = length(y); R = Vector{Vector{Float64}}()
    for _ in 1:B
        i = rand(rng, 1:n, n)
        m = try ns_media_akaike(y[i], x[i]; b, formas) catch; continue end
        push!(R, [m.p_em(x0, xi) for xi in xs])
    end
    return permutedims(reduce(hcat, R))
end

"P(N ≥ k) para N = Σ Bernoulli(pᵢ) independentes (Poisson-binomial), por programação dinâmica."
function prob_ao_menos_k(ps, k)
    f = [1.0]                                         # f[j+1] = P(N = j)
    for p in ps
        f = [f; 0.0] .* (1 - p) .+ [0.0; f] .* p
    end
    return sum(f[k+1:end])
end

"""
P(≥ 1 ano com chuva ≤ x₀ nos próximos n anos) quando p varia ano a ano:
    1 − Πⱼ (1 − pⱼ).
`pfun(j)` devolve p no ano j = 1…n.
"""
prob_n_variavel(pfun, n) = 1 - prod(1 - pfun(j) for j in 1:n)

"""
Gamma não estacionária bayesiana: log μᵢ = a + b·xᵢ, x em décadas centradas.
    log k ~ Normal(log 4, 1);  a ~ Normal(log 500, 1);  b ~ Normal(0, 0,2)
A priori de b admite ±40 %/década com 95 % de probabilidade — muito mais larga
que qualquer tendência plausível de chuva anual, mas exclui valores absurdos.
"""
@model function modelo_gamma_ns(y, x)
    logk ~ Normal(log(4.0), 1.0)
    a ~ Normal(log(500.0), 1.0)
    b ~ Normal(0.0, 0.2)
    k = exp(logk)
    θ = [exp(a + b * xi) / k for xi in x]
    if !(isfinite(k) && k > 0 && all(t -> isfinite(t) && t > 0, θ))
        Turing.@addlogprob! -Inf          # rejeita propostas numericamente inválidas no aquecimento
        return
    end
    for i in eachindex(y)
        y[i] ~ Gamma(k, θ[i])
    end
end

function posterior_ns(y, x; n = 4000, cadeias = 4, rng = MersenneTwister(SEMENTE))
    ch = sample(rng, modelo_gamma_ns(y, x), NUTS(0.9), MCMCThreads(), n ÷ cadeias, cadeias; progress = false)
    return (logk = vec(Array(ch[:logk])), a = vec(Array(ch[:a])), b = vec(Array(ch[:b])), ch = ch)
end

# ----------------------------------------------------------------------------- diagnósticos bayesianos
"""
Tamanho efetivo de amostra (ESS) de uma matriz iterações × cadeias:
ESS = N / (1 + 2 Σₖ ρ̂ₖ), com ρ̂ₖ a autocorrelação média entre cadeias, somada
enquanto os pares (ρ̂₂ₖ₋₁ + ρ̂₂ₖ) forem positivos (sequência inicial positiva de Geyer).
"""
function ess(M::AbstractMatrix)
    n, m = size(M)
    x = M .- mean(M, dims = 1); v = vec(mean(x .^ 2, dims = 1))
    ρ(k) = mean(sum(x[1:n-k, j] .* x[1+k:n, j]) / (n * v[j]) for j in 1:m)
    s = 0.0; k = 1
    while k + 1 < n
        par = ρ(k) + ρ(k + 1)
        par <= 0 && break
        s += par; k += 2
    end
    return n * m / (1 + 2s)
end
ess(ch, nome) = (A = Array(ch[nome]); ess(ndims(A) == 1 ? reshape(A, :, 1) : A))

"""
Checagem preditiva posterior do mínimo: para cada draw θ⁽ˢ⁾, simula uma amostra
do mesmo tamanho (e, no modelo não estacionário, nos mesmos anos) e guarda o
mínimo. p_ppc = P(min_sim ≤ min_obs). Valores próximos de 0 ou 1 indicam que o
modelo não reproduz a cauda seca observada.
"""
function ppc_minimo(dist_por_draw, nobs, min_obs; rng = MersenneTwister(SEMENTE))
    mins = [minimum(rand(rng, d) for d in dist_por_draw(s)) for s in 1:nobs.ndraws]
    return mean(mins .<= min_obs), mins
end

# ----------------------------------------------------------------------------- análise regional de frequência
"""
    regional_candidata(D, y, x, s; nsites, forma_var = false, b_por_site = false)

Modelo de índice de cheia com tendência: Yᵢₜ = μᵢ·Zₜ, com Z da família `D` com
posição log-linear no tempo (como em `ns_candidata`) e forma fixa ou log-linear no
tempo, comum a todas as estações. `s` é o índice da estação de cada observação e
μ₁ = 1 (a estação 1 é a de referência). Com `b_por_site`, cada estação tem sua
própria inclinação (usado para testar a homogeneidade da tendência).
Parâmetros: [log s₀, a, b (ou b₁…b_S), (c), log μ₂ … log μ_S].
"""
function regional_candidata(D, y, x, s; nsites, forma_var = false, b_por_site = false, θ0 = nothing)
    nb = b_por_site ? nsites : 1
    ib = 3:(2 + nb); ic = forma_var ? 3 + nb : 0
    iμ = (3 + nb + (forma_var ? 1 : 0)):(3 + nb + (forma_var ? 1 : 0) + nsites - 2)
    function distz(θ, xi, si)
        η = θ[2] + θ[b_por_site ? ib[si] : ib[1]] * xi
        sh = exp(θ[1] + (forma_var ? θ[ic] * xi : 0.0))
        D === Normal    ? Normal(exp(η), sh) :
        D === LogNormal ? LogNormal(η, sh) :
        D === Gamma     ? Gamma(sh, exp(η) / sh) :
        D === Weibull   ? Weibull(sh, exp(η)) : error("família não suportada: $D")
    end
    logμ(θ, si) = si == 1 ? 0.0 : θ[iμ[si - 1]]
    function nll(θ)
        v = 0.0
        for (yi, xi, si) in zip(y, x, s)
            lm = logμ(θ, si)
            v -= logpdf(distz(θ, xi, si), yi / exp(lm)) - lm
        end
        isfinite(v) ? v : 1e12
    end
    if isnothing(θ0)
        # partida: índices = médias por estação relativas à de referência; Z = y/μ
        m1 = mean(y[s .== 1]); μ0 = [mean(y[s .== k]) / m1 for k in 1:nsites]
        z = y ./ μ0[s]
        d0 = ajusta_mle(D, z)
        s0 = D === Normal ? std(z) : D === LogNormal ? std(log.(z)) : shape(d0)
        a0 = D === LogNormal ? mean(log.(z)) : D === Weibull ? log(scale(d0)) : log(mean(z))
        θ0 = [log(s0); a0; zeros(nb); forma_var ? 0.0 : Float64[]; log.(μ0[2:end])]
    end
    r = optimize(nll, θ0, NelderMead(), Optim.Options(iterations = 40_000, g_tol = 1e-12))
    r = optimize(nll, Optim.minimizer(r), BFGS())
    θ = Optim.minimizer(r)
    return (D = D, forma_var = forma_var, θ = θ, ll = -Optim.minimum(r), npar = length(θ),
            b = b_por_site ? θ[ib] : θ[ib[1]], c = forma_var ? θ[ic] : 0.0,
            μ = [1.0; exp.(θ[iμ])],
            F = (z0, xi) -> cdf(distz(θ, xi, 1), z0))     # na escala da estação de referência
end

"""
Média por pesos de Akaike das candidatas regionais. Com `estacionario = true`, a
posição e a forma não variam no tempo (x = 0 para todos). Os pesos usam a
verossimilhança combinada, que trata as estações como independentes; a incerteza
vem do bootstrap por ano (`regional_media_akaike_boot`), que preserva a correlação
entre estações.
"""
function regional_media_akaike(y, x, s; nsites, estacionario = false, candidatas = CANDIDATAS,
                               formas = (false, true))
    xx = estacionario ? zero(x) : x
    fs = estacionario ? (false,) : formas
    fits = [regional_candidata(D, y, xx, s; nsites, forma_var = fv) for fv in fs for D in candidatas]
    npar = [f.npar - (estacionario ? 1 : 0) for f in fits]           # sem tendência, b não conta
    aic = [2k - 2f.ll for (k, f) in zip(npar, fits)]
    w = exp.(-(aic .- minimum(aic)) ./ 2); w ./= sum(w)
    p_em(z0, xi) = sum(w[j] * fits[j].F(z0, estacionario ? 0.0 : xi) for j in eachindex(fits))
    tab = DataFrame(distribuicao = [string(nameof(f.D)) for f in fits],
                    forma = [f.forma_var ? "variável" : "fixa" for f in fits],
                    efeito_pct = [estacionario ? 0.0 : 100(exp(f.b) - 1) for f in fits],
                    loglik = [f.ll for f in fits], AIC = aic, dAIC = aic .- minimum(aic), peso_akaike = w)
    return (fits = fits, tab = tab, w = w, p_em = p_em)
end

"""
Bootstrap por ano: sorteia anos hidrológicos com reposição e leva todas as
observações das estações naquele ano, preservando a correlação espacial. Devolve
p = F(z₀) no ponto `xa` (não estacionário) e o p estacionário, por réplica.
"""
function regional_media_akaike_boot(y, x, s, ano, z0, xa; nsites, B = 300, rng = MersenneTwister(SEMENTE))
    anos = unique(ano); idx = Dict(a => findall(==(a), ano) for a in anos)
    sorteios = [rand(rng, anos, length(anos)) for _ in 1:B]
    out = Vector{Union{Nothing, NTuple{2, Float64}}}(nothing, B)
    Threads.@threads for r in 1:B
        i = reduce(vcat, [idx[a] for a in sorteios[r]])
        length(unique(s[i])) == nsites || continue                    # todas as estações precisam aparecer
        out[r] = try
            mn = regional_media_akaike(y[i], x[i], s[i]; nsites)
            me = regional_media_akaike(y[i], x[i], s[i]; nsites, estacionario = true)
            (mn.p_em(z0, xa), me.p_em(z0, 0.0))
        catch
            nothing
        end
    end
    ok = filter(!isnothing, out)
    return (ns = [o[1] for o in ok], est = [o[2] for o in ok])
end
