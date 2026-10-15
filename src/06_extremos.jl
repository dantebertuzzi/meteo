# =============================================================================
# 06_extremos.jl — anos anômalos, chuva diária extrema (POT/GPD) e índices de contagem
#
# • Anos secos/chuvosos: SPI anual (McKee et al., 1993). Ajusta-se uma Gamma
#   por MV aos totais do ano hidrológico e transforma-se zᵢ = Φ⁻¹(F_Gamma(xᵢ)).
#   Isso corrige a assimetria da chuva: um desvio de +300 mm e um de −300 mm não
#   são igualmente raros. Categorias: |z| ≥ 1 moderado, ≥ 1,5 severo, ≥ 2 extremo.
# • Anos quentes/frios: anomalia padronizada z = (x − x̄)/s. A temperatura anual
#   é aproximadamente simétrica (ver Shapiro–Wilk), então dispensa transformação.
# • Chuva diária extrema: picos acima de limiar (POT) com Pareto Generalizada.
#   Com ~16 anos válidos, uma GEV de máximos anuais teria só ~16 pontos; o POT usa
#   todos os eventos acima do limiar (dezenas), reduzindo muito a variância.
#   F(y) = 1 − (1 + ξy/σ)^{−1/ξ},  y = x − u > 0.
#   Nível de retorno de T anos: x_T = u + (σ/ξ)[(λT)^ξ − 1], λ = eventos/ano.
# =============================================================================

function spi_anual(x)
    ok = .!ismissing.(x); g = fit_mle(Gamma, Float64.(x[ok]))
    z = Vector{Union{Missing, Float64}}(missing, length(x))
    z[ok] = quantile.(Normal(), cdf.(g, Float64.(x[ok])))
    return z, g
end

function classe_spi(z)
    ismissing(z) && return "sem dados suficientes"
    z <= -2 ? "extremamente seco" : z <= -1.5 ? "severamente seco" : z <= -1 ? "moderadamente seco" :
    z >= 2 ? "extremamente chuvoso" : z >= 1.5 ? "severamente chuvoso" : z >= 1 ? "moderadamente chuvoso" : "normal"
end

anomalia_padronizada(x) = (μ = mean(skipmissing(x)); s = std(skipmissing(x)); [ismissing(v) ? missing : (v - μ) / s for v in x])

"Cercas de Tukey: fora de [Q₁ − 1,5·IQR, Q₃ + 1,5·IQR] é valor atípico (descritivo, não removido)."
function cercas_tukey(x)
    q1, q3 = quantile(collect(skipmissing(x)), [0.25, 0.75]); h = 1.5(q3 - q1)
    return q1 - h, q3 + h
end

# ----------------------------------------------------------------------------- otimização
"Nelder–Mead mínimo (evita dependência de Optim.jl para um problema de 2 parâmetros)."
function nelder_mead(f, x0; passo = 0.2, tol = 1e-10, maxit = 10_000)
    n = length(x0)
    S = [copy(float.(x0)) for _ in 1:n+1]
    for i in 1:n
        S[i+1][i] += passo
    end
    F = f.(S)
    for _ in 1:maxit
        o = sortperm(F); S = S[o]; F = F[o]
        abs(F[end] - F[1]) < tol && break
        c = sum(S[1:n]) / n
        xr = c .+ (c .- S[end]); fr = f(xr)
        if fr < F[1]
            xe = c .+ 2 .* (c .- S[end]); fe = f(xe)
            fe < fr ? (S[end] = xe; F[end] = fe) : (S[end] = xr; F[end] = fr)
        elseif fr < F[n]
            S[end] = xr; F[end] = fr
        else
            xc = fr < F[end] ? c .+ 0.5 .* (xr .- c) : c .+ 0.5 .* (S[end] .- c); fc = f(xc)
            if fc < min(fr, F[end])
                S[end] = xc; F[end] = fc
            else
                for i in 2:n+1
                    S[i] = S[1] .+ 0.5 .* (S[i] .- S[1]); F[i] = f(S[i])
                end
            end
        end
    end
    i = argmin(F)
    return S[i], F[i]
end

"Hessiana por diferenças finitas centradas."
function hessiana(f, x; h = 1e-4)
    n = length(x); H = zeros(n, n)
    for i in 1:n, j in 1:n
        ei = zeros(n); ej = zeros(n); ei[i] = h; ej[j] = h
        H[i, j] = (f(x .+ ei .+ ej) - f(x .+ ei .- ej) - f(x .- ei .+ ej) + f(x .- ei .- ej)) / (4h^2)
    end
    return H
end

# ----------------------------------------------------------------------------- GPD
"−log-verossimilhança da GPD para excessos y > 0; θ = (σ, ξ)."
function gpd_nll(θ, y)
    σ, ξ = θ
    σ <= 0 && return Inf
    n = length(y)
    abs(ξ) < 1e-6 && return n * log(σ) + sum(y) / σ
    z = 1 .+ ξ .* y ./ σ
    any(<=(0), z) && return Inf
    return n * log(σ) + (1 + 1 / ξ) * sum(log, z)
end

"Ajuste da GPD por MV; EP pelos inversos da informação observada."
function ajusta_gpd(y)
    f(θ) = gpd_nll(θ, y)
    θ, nll = nelder_mead(f, [mean(y), 0.1]; passo = 0.1 * mean(y))
    H = hessiana(f, θ; h = 1e-4 * max(1.0, θ[1]))
    V = try inv(H) catch; fill(NaN, 2, 2) end
    ep = sqrt.(max.(diag(V), 0.0))
    return (σ = θ[1], ξ = θ[2], ep_σ = ep[1], ep_ξ = ep[2], nll = nll, n = length(y))
end

"""
Desagrupamento por corridas: dias acima de `u` separados por ≤ `r` dias formam um
cluster; mantém-se o máximo. Isso aproxima a independência exigida pelo POT.
"""
function desagrupa(datas, x, u; r = 1)
    idx = findall(>(u), x); isempty(idx) && return Float64[], Date[]
    picos = Float64[]; dpico = Date[]
    atual = [idx[1]]
    for k in idx[2:end]
        if Dates.value(datas[k] - datas[atual[end]]) <= r
            push!(atual, k)
        else
            m = atual[argmax(x[atual])]; push!(picos, x[m]); push!(dpico, datas[m]); atual = [k]
        end
    end
    m = atual[argmax(x[atual])]; push!(picos, x[m]); push!(dpico, datas[m])
    return picos, dpico
end

"""
Exposição efetiva em 'anos climatológicos': Σ_dias válidos c̄(mês)/P̄_anual.
Como as falhas se concentram na estação chuvosa, dividir por (dias válidos/365)
subestimaria a taxa λ; ponderar pela chuva climatológica do mês corrige isso.
"""
anos_efetivos(dd, clim_mes, clim_anual) =
    sum(clim_mes[m] / daysinmonth(Date(2001, m)) for (m, p) in zip(dd.mes, dd.prec) if !ismissing(p)) / clim_anual

"Diagnósticos de escolha de limiar: vida residual média e estabilidade de ξ̂ e σ* = σ − ξu."
function diagnostico_limiar(datas, x, limiares; r = 1)
    DataFrame(map(limiares) do u
        pk, _ = desagrupa(datas, x, u; r = r); y = pk .- u
        g = length(y) >= 10 ? ajusta_gpd(y) : (σ = NaN, ξ = NaN, ep_σ = NaN, ep_ξ = NaN)
        (u = u, n_exc = length(y), mrl = mean(y), mrl_ep = std(y) / sqrt(length(y)),
         ξ = g.ξ, ξ_ep = g.ep_ξ, σ_mod = g.σ - g.ξ * u)
    end)
end

nivel_retorno(u, σ, ξ, λ, T) = abs(ξ) < 1e-6 ? u + σ * log(λ * T) : u + σ / ξ * ((λ * T)^ξ - 1)

"Níveis de retorno com IC por bootstrap (reamostra excessos; nº de eventos ~ Poisson)."
function niveis_retorno(y, u, λ, anos, Ts; B = NBOOT, rng = MersenneTwister(SEMENTE))
    g = ajusta_gpd(y); n = length(y)
    pontual = [nivel_retorno(u, g.σ, g.ξ, λ, T) for T in Ts]
    boot = fill(NaN, B, length(Ts))
    for b in 1:B
        nb = max(rand(rng, Poisson(n)), 5)
        yb = y[rand(rng, 1:n, nb)]
        gb = try ajusta_gpd(yb) catch; continue end
        boot[b, :] = [nivel_retorno(u, gb.σ, gb.ξ, nb / anos, T) for T in Ts]
    end
    ic = [quantile(filter(isfinite, boot[:, j]), [0.025, 0.975]) for j in eachindex(Ts)]
    DataFrame(T_anos = Ts, nivel_mm = pontual, ic_inf = first.(ic), ic_sup = last.(ic),
              extrapolacao = Ts .> 2anos)
end

# ----------------------------------------------------------------------------- contagens
"""
Tendência em contagens: GLM Poisson com offset log(exposição) (ou Binomial com
n tentativas). Sobredispersão estimada por φ̂ = Σr²_Pearson/(n − p); se φ̂ > 1 os
EP são multiplicados por √φ̂ (quase-verossimilhança). Efeito: razão de taxas por década.
"""
function tendencia_contagem(t, k, n; familia = :poisson, nome = "")
    ok = .!ismissing.(k) .& .!ismissing.(n)
    df = DataFrame(t = Float64.(t[ok]) .- mean(t[ok]), k = Float64.(k[ok]), n = Float64.(n[ok]))
    if familia == :poisson
        m = glm(@formula(k ~ t), df, Poisson(), LogLink(); offset = log.(df.n))
        μ = predict(m, df; offset = log.(df.n)); rp = (df.k .- μ) ./ sqrt.(μ)
    else
        df.prop = df.k ./ df.n
        m = glm(@formula(prop ~ t), df, Binomial(), LogitLink(); wts = df.n)
        π̂ = predict(m); rp = (df.k .- df.n .* π̂) ./ sqrt.(df.n .* π̂ .* (1 .- π̂))
    end
    φ = sum(abs2, rp) / (nrow(df) - 2)
    β = coef(m)[2]; ep = stderror(m)[2] * sqrt(max(φ, 1.0))
    z = β / ep; tq = quantile(TDist(nrow(df) - 2), 0.975)
    (indice = nome, familia = String(familia), n_anos = nrow(df), dispersao = φ,
     razao_decada = exp(10β), ic_inf = exp(10(β - tq * ep)), ic_sup = exp(10(β + tq * ep)),
     p = 2ccdf(TDist(nrow(df) - 2), abs(z)), modelo = m)
end
