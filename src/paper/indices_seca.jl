# =============================================================================
# indices_seca.jl — SPI, ETo de Hargreaves–Samani, SPEI e teoria de corridas
#
# SPI-k (McKee et al., 1993; Guttman, 1999):
#   Xₜ = Σ_{j=0}^{k−1} Pₜ₋ⱼ (soma móvel de k meses; exige os k meses válidos).
#   Para cada mês calendário m ajusta-se G = Gamma(α, β) por MV aos Xₜ > 0 do
#   período de referência; com q = P(X = 0):
#       H(x) = q + (1 − q)·G(x),   SPIₜ = Φ⁻¹(H(Xₜ)).
#   O ajuste por mês calendário remove a sazonalidade: dezembro é comparado só
#   com dezembros.
#
# ETo de Hargreaves & Samani (1985), na forma da FAO-56 (eq. 52):
#   ETo = 0,0023 · 0,408·Ra · (T̄ + 17,8) · √(Tmáx − Tmín)        [mm/dia]
#   Ra = radiação no topo da atmosfera (FAO-56, eq. 21), função só da latitude
#   e do dia do ano — por isso dispensa a radiação medida, que foi excluída.
#
# SPEI-k (Vicente-Serrano et al., 2010): mesma construção do SPI sobre o
# balanço Dₜ = Pₜ − EToₜ, que pode ser negativo; usa-se a log-logística de 3
# parâmetros ajustada por momentos ponderados por probabilidade (PWM):
#   F(x) = [1 + (α/(x − γ))^β]⁻¹.
# Stagge et al. (2015) recomendam a GEV; as duas são comparadas.
#
# Teoria de corridas (Yevjevich, 1967): um evento é uma sequência de meses com
# índice < z₀. Duração D = nº de meses; severidade S = −Σ zₜ; intensidade I = S/D.
# =============================================================================

# ----------------------------------------------------------------------------- somas móveis
"Soma móvel de `k` meses; `missing` se algum dos `k` meses for ausente."
function soma_movel(x::AbstractVector, k::Int)
    n = length(x); out = Vector{Union{Missing, Float64}}(missing, n)
    for t in k:n
        janela = @view x[t-k+1:t]
        any(ismissing, janela) || (out[t] = sum(janela))
    end
    return out
end

# ----------------------------------------------------------------------------- SPI
"""
    spi(prec, meses, anos; k, ref = (1991, 2020))

`prec`: série mensal contínua (sem saltos de calendário; ausentes = missing).
Ajuste por mês calendário usando só os anos de `ref` (período normal da OMM).
Devolve (z, ajustes) — `ajustes[m] = (q, Gamma)`.
"""
function spi(prec, meses, anos; k::Int, ref = (1991, 2020))
    X = soma_movel(prec, k)
    z = Vector{Union{Missing, Float64}}(missing, length(X))
    ajustes = Dict{Int, Any}()
    for m in 1:12
        idx = findall(i -> meses[i] == m && !ismissing(X[i]), eachindex(X))
        cal = [X[i] for i in idx if ref[1] <= anos[i] <= ref[2]]
        length(cal) < 10 && continue
        pos = filter(>(0), cal); q = 1 - length(pos) / length(cal)
        length(pos) < 5 && continue
        g = fit_mle(Gamma, pos)
        ajustes[m] = (q = q, gamma = g, n = length(cal))
        for i in idx
            H = q + (1 - q) * cdf(g, X[i])
            z[i] = quantile(Normal(), clamp(H, 1e-6, 1 - 1e-6))
        end
    end
    return z, ajustes
end

# ----------------------------------------------------------------------------- Hargreaves
"Radiação extraterrestre diária Ra (MJ m⁻² d⁻¹), FAO-56 eqs. 21–25. `lat` em graus."
function radiacao_extraterrestre(lat, J)
    φ = deg2rad(lat)
    dr = 1 + 0.033cos(2π * J / 365)
    δ = 0.409sin(2π * J / 365 - 1.39)
    ωs = acos(clamp(-tan(φ) * tan(δ), -1.0, 1.0))
    return 24 * 60 / π * 0.0820 * dr * (ωs * sin(φ) * sin(δ) + cos(φ) * cos(δ) * sin(ωs))
end

"ETo de Hargreaves–Samani (mm/dia), FAO-56 eq. 52."
eto_hargreaves(tmax, tmin, Ra) = 0.0023 * 0.408Ra * ((tmax + tmin) / 2 + 17.8) * sqrt(max(tmax - tmin, 0.0))

"""
ETo mensal (mm/mês) a partir das médias mensais de Tmáx e Tmín e da média de Ra
nos dias do mês. É a forma mensal usual da equação (Hargreaves & Samani, 1985).
"""
function eto_mensal(ano, mes, tmax, tmin, lat)
    (ismissing(tmax) || ismissing(tmin)) && return missing
    dias = Date(ano, mes, 1):Day(1):lastdayofmonth(Date(ano, mes, 1))
    Ra = mean(radiacao_extraterrestre(lat, dayofyear(d)) for d in dias)
    return length(dias) * eto_hargreaves(tmax, tmin, Ra)
end

# ----------------------------------------------------------------------------- log-logística
"PWM não viesados αₛ = E[X(1 − F)ˢ] (Hosking, 1986), s = 0, 1, 2."
function pwm_alfa(x)
    x = sort(x); n = length(x)
    [sum(binomial(n - i, s) / binomial(n - 1, s) * x[i] for i in 1:n) / n for s in 0:2]
end

"""
Log-logística de 3 parâmetros. `espelhada = true` quando a amostra tem
assimetria negativa: nesse caso o PWM dá β < 0 (distribuição limitada por
cima), e o ajuste é feito em −x, com F(x) = 1 − F₋(−x).
"""
struct LogLogistica3
    α::Float64; β::Float64; γ::Float64; espelhada::Bool
end
LogLogistica3(α, β, γ) = LogLogistica3(α, β, γ, false)

function _pwm_ll(x)
    w0, w1, w2 = pwm_alfa(x)
    β = (2w1 - w0) / (6w1 - w0 - 6w2)
    Γ = gamma(1 + 1 / β) * gamma(1 - 1 / β)
    α = (w0 - 2w1) * β / Γ
    γ = w0 - α * Γ
    return α, β, γ
end

"Ajuste por PWM conforme Vicente-Serrano et al. (2010), eqs. 4–6."
function ajusta_loglogistica(x)
    α, β, γ = _pwm_ll(x)
    β > 1 && return LogLogistica3(α, β, γ, false)
    α, β, γ = _pwm_ll(-x)
    β > 1 || error("log-logística sem solução por PWM (β = $β)")
    return LogLogistica3(α, β, γ, true)
end

_cdf_ll(α, β, γ, x) = x <= γ ? 0.0 : 1 / (1 + (α / (x - γ))^β)
cdf_ll(d::LogLogistica3, x) = d.espelhada ? 1 - _cdf_ll(d.α, d.β, d.γ, -x) : _cdf_ll(d.α, d.β, d.γ, x)
function quantil_ll(d::LogLogistica3, p)
    d.espelhada && return -(d.γ + d.α * ((1 - p) / p)^(1 / d.β))
    d.γ + d.α * (p / (1 - p))^(1 / d.β)
end

"GEV por L-momentos (Hosking, 1990). Devolve GeneralizedExtremeValue de Distributions (ξ = −k)."
function ajusta_gev_lmom(x)
    a0, a1, a2 = pwm_alfa(x)
    b0, b1, b2 = a0, (a0 - a1), (a0 - 2a1 + a2)     # βᵣ = PWM em F (Hosking): β₁ = α₀ − α₁ …
    λ1 = b0; λ2 = 2b1 - b0; λ3 = 6b2 - 6b1 + b0
    t3 = λ3 / λ2
    c = 2 / (3 + t3) - log(2) / log(3)
    k = 7.8590c + 2.9554c^2
    σ = λ2 * k / ((1 - 2.0^(-k)) * gamma(1 + k))
    μ = λ1 - σ * (1 - gamma(1 + k)) / k
    return GeneralizedExtremeValue(μ, σ, -k)
end

"""
    spei(balanco, meses, anos; k, ref, dist = :loglogistica)

Como `spi`, mas sobre o balanço hídrico (pode ser negativo) e com distribuição
de 3 parâmetros (`:loglogistica` ou `:gev`).
"""
function spei(balanco, meses, anos; k::Int, ref = (1991, 2020), dist = :loglogistica)
    X = soma_movel(balanco, k)
    z = Vector{Union{Missing, Float64}}(missing, length(X))
    ajustes = Dict{Int, Any}()
    for m in 1:12
        idx = findall(i -> meses[i] == m && !ismissing(X[i]), eachindex(X))
        cal = Float64[X[i] for i in idx if ref[1] <= anos[i] <= ref[2]]
        length(cal) < 10 && continue
        if dist == :loglogistica
            d = ajusta_loglogistica(cal); F = x -> cdf_ll(d, x)
        else
            d = ajusta_gev_lmom(cal); F = x -> cdf(d, x)
        end
        ajustes[m] = d
        for i in idx
            z[i] = quantile(Normal(), clamp(F(X[i]), 1e-6, 1 - 1e-6))
        end
    end
    return z, ajustes
end

# ----------------------------------------------------------------------------- corridas
"""
    corridas(datas, z; z0 = -1.0)

Eventos = meses consecutivos com z < z0. Um mês ausente encerra o evento e o
marca como `censurado` (a duração real pode ser maior).
"""
function corridas(datas, z; z0 = -1.0)
    ev = NamedTuple[]
    i = 1; n = length(z)
    while i <= n
        if !ismissing(z[i]) && z[i] < z0
            j = i
            while j < n && !ismissing(z[j+1]) && z[j+1] < z0
                j += 1
            end
            zz = Float64.(z[i:j]); S = -sum(zz)
            cens = (i > 1 && ismissing(z[i-1])) || (j < n && ismissing(z[j+1])) || j == n
            push!(ev, (inicio = datas[i], fim = datas[j], duracao = j - i + 1, severidade = S,
                       intensidade = S / (j - i + 1), pico = minimum(zz), censurado = cens))
            i = j + 1
        else
            i += 1
        end
    end
    return DataFrame(ev)
end
