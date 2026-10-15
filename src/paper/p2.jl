# =============================================================================
# p2.jl — mudança no regime: sazonalidade, índices ETCCDI, tendências com
# autocorrelação, poder, regressão quantílica e decomposição STL.
# =============================================================================
using QuantileRegressions: qreg
using Loess: loess, predict as lpredict

# ============================================================================= sazonalidade
"""
Índice de sazonalidade de Walsh & Lawler (1981) de um ano:
    SI = (1/R)·Σₘ |xₘ − R/12|,   R = total anual.
0 = chuva uniforme; ≈ 1,83 = toda a chuva num só mês.
"""
indice_sazonalidade(x12) = (R = sum(x12); sum(abs.(x12 .- R / 12)) / R)

"""
Estatística circular da chuva diária num ano hidrológico de N dias:
θ_d = 2π(d − ½)/N; C = Σ P_d cos θ_d; S = Σ P_d sin θ_d.
Centroide (dia) = N·atan(S, C)/(2π) (mod N); concentração r = √(C² + S²)/Σ P_d ∈ [0, 1]
(0 = chuva espalhada uniformemente; 1 = toda num só dia). O ano hidrológico começa
em julho, então o centroide (≈ fev–mar) fica longe da costura e pode ser tratado
como variável linear nos testes de tendência.
"""
function circular_chuva(p::AbstractVector)
    N = length(p); θ = 2π .* ((1:N) .- 0.5) ./ N
    C = sum(p .* cos.(θ)); S = sum(p .* sin.(θ)); T = sum(p)
    ang = atan(S, C); ang < 0 && (ang += 2π)
    return (centroide = N * ang / (2π), concentracao = sqrt(C^2 + S^2) / T)
end

"""
Início e fim da estação chuvosa (Liebmann et al., 2007): anomalia acumulada
    A(d) = Σ_{n ≤ d} (Rₙ − R̄),   R̄ = chuva diária média climatológica.
Início = dia seguinte ao mínimo de A; fim = dia do máximo de A após o início.
Dias ausentes entram como R̄ (anomalia nula), o que não desloca o mínimo nem o máximo.
`p` é a série diária do ano hidrológico (jul–jun).
"""
function liebmann(p::AbstractVector, Rbar)
    a = cumsum([ismissing(x) ? 0.0 : x - Rbar for x in p])
    i0 = argmin(a)
    i1 = i0 - 1 + argmax(a[i0:end])
    return (inicio = i0 + 1, fim = i1, duracao = i1 - i0)
end

# ============================================================================= ETCCDI
"Maior sequência de dias consecutivos que satisfazem `f`; `missing` interrompe a corrida."
function maior_corrida(x, f)
    best = run = 0
    for v in x
        run = (!ismissing(v) && f(v)) ? run + 1 : 0
        best = max(best, run)
    end
    return best
end

"""
Índices de chuva por ano (Zhang et al., 2011), sobre os dias do ano:
R1mm (dias ≥ 1 mm), Rx1day, SDII (total dos dias úmidos / nº de dias úmidos),
R95p (total em dias acima do P95 dos dias úmidos do período-base), CDD e CWD.
"""
function indices_chuva(p, p95)
    pv = collect(skipmissing(p)); u = pv[pv .>= 1]
    (R1mm = length(u), Rx1day = isempty(pv) ? missing : maximum(pv),
     SDII = isempty(u) ? missing : sum(u) / length(u), R95p = sum(pv[pv .> p95]; init = 0.0),
     CDD = maior_corrida(p, <(1.0)), CWD = maior_corrida(p, >=(1.0)))
end

"""
Percentil por dia do calendário com janela de ±2 dias, no período-base
(Zhang et al., 2005, para TX90p/TN10p). Devolve Dict(dia_do_ano => limiar).
"""
function percentil_calendario(datas, x, q; base = (1991, 2020), janela = 2)
    porDia = Dict{Int, Vector{Float64}}(d => Float64[] for d in 1:366)
    for (dt, v) in zip(datas, x)
        (ismissing(v) || !(base[1] <= year(dt) <= base[2])) && continue
        d = dayofyear(dt)
        for k in -janela:janela
            push!(porDia[mod1(d + k, 366)], v)
        end
    end
    Dict(d => (length(v) >= 30 ? quantile(v, q) : NaN) for (d, v) in porDia)
end

# ============================================================================= tendência
"""
MK com variância corrigida de Hamed & Rao (1998): Var*(S) = Var(S)·n/n*,
    n/n* = 1 + 2/[n(n−1)(n−2)] · Σ_{k=1}^{n−1} (n−k)(n−k−1)(n−k−2)·ρₖ,
com ρₖ a autocorrelação dos postos da série sem tendência (Sen), usando só os
lags significativos a 5 %. Para séries com lacunas usa-se a série compactada.
"""
function mk_hamed_rao(t, y)
    o = sortperm(t); t = Float64.(t[o]); y = Float64.(y[o]); n = length(y)
    mk = mann_kendall(t, y); ts = theil_sen(t, y)
    r = tiedrank(y .- ts.declive .* t)
    rc = r .- mean(r); den = sum(abs2, rc)
    soma = 0.0
    for k in 1:n-1
        ρ = sum(rc[1:n-k] .* rc[1+k:n]) / den
        abs(ρ) > quantile(Normal(), 0.975) / sqrt(n) || continue
        soma += (n - k) * (n - k - 1) * (n - k - 2) * ρ
    end
    fator = max(1 + 2soma / (n * (n - 1) * (n - 2)), 1e-3)
    v = mk.varS * fator
    Z = mk.S == 0 ? 0.0 : (mk.S - sign(mk.S)) / sqrt(v)
    return (fator = fator, Z = Z, p = 2ccdf(Normal(), abs(Z)))
end

"""
Tendência mínima detectável (poder 80 %, MK bilateral a 5 %) por simulação:
yₜ = β·t + eₜ, eₜ ~ AR(1) com ρ e σ estimados dos resíduos de Sen da série real.
Busca β numa grade e devolve o menor com poder ≥ 0,8 (em unidades/década).
"""
function tendencia_minima_detectavel(t, y; poder = 0.8, B = 400, rng = MersenneTwister(SEMENTE))
    o = sortperm(t); t = Float64.(t[o]); y = Float64.(y[o]); n = length(y)
    ts = theil_sen(t, y); e = y .- ts.declive .* t .- ts.intercepto
    ρ = clamp(cor(e[1:end-1], e[2:end]), -0.5, 0.9); σ = std(e)
    simula(β) = begin
        rej = 0
        for _ in 1:B
            u = zeros(n); u[1] = randn(rng) * σ
            for i in 2:n
                u[i] = ρ * u[i-1] + sqrt(1 - ρ^2) * σ * randn(rng)
            end
            rej += mann_kendall(t, β .* (t .- t[1]) .+ u).p < 0.05
        end
        rej / B
    end
    escala = σ / (t[end] - t[1])
    for f in 0.25:0.25:12
        simula(f * escala) >= poder && return (mdt_decada = 10f * escala, rho = ρ, sigma = σ)
    end
    return (mdt_decada = NaN, rho = ρ, sigma = σ)
end

"Bateria completa para uma série anual: MK, TFPW, Hamed–Rao, Sen (IC), Pettitt e MDT."
function tendencia_completa(t, y; nome = "", escala = 10)
    ok = .!ismissing.(y); t = Float64.(t[ok]); y = Float64.(y[ok])
    a = analise_tendencia(t, y; nome = nome, escala = escala)
    hr = mk_hamed_rao(t, y); md = tendencia_minima_detectavel(t, y)
    return merge(a, (p_hamed_rao = hr.p, fator_hr = hr.fator, mdt_decada = md.mdt_decada))
end

# ============================================================================= regressão quantílica
"""
Regressão quantílica da temperatura diária no tempo, controlando a sazonalidade
por 2 harmônicos: Q_τ(y | t, d) = β₀ + β₁·t + Σₖ [aₖ sin(2πkd/365,25) + bₖ cos(…)].
β₁ (°C/década) é a mudança do quantil τ. Como os dias são autocorrelacionados,
o IC vem de bootstrap em blocos de um ano (reamostram-se anos inteiros).
"""
function quantilica_tempo(datas, y; taus = (0.1, 0.5, 0.9), B = 300, rng = MersenneTwister(SEMENTE))
    ok = .!ismissing.(y)
    df = DataFrame(y = Float64.(y[ok]), data = datas[ok])
    df.t = (year.(df.data) .+ (dayofyear.(df.data) .- 0.5) ./ 365.25 .- 2000) ./ 10
    w = 2π .* dayofyear.(df.data) ./ 365.25
    df.s1 = sin.(w); df.c1 = cos.(w); df.s2 = sin.(2w); df.c2 = cos.(2w)
    df.ano = year.(df.data)
    anos = unique(df.ano); blocos = Dict(a => findall(==(a), df.ano) for a in anos)
    f = @formula(y ~ t + s1 + c1 + s2 + c2)
    map(taus) do τ
        b = coef(qreg(f, df, τ))[2]
        bs = Float64[]
        for _ in 1:B
            idx = reduce(vcat, [blocos[a] for a in rand(rng, anos, length(anos))])
            push!(bs, coef(qreg(f, df[idx, :], τ))[2])
        end
        (tau = τ, beta_decada = b, ic_inf = quantile(bs, 0.025), ic_sup = quantile(bs, 0.975), n_dias = nrow(df), n_anos = length(anos))
    end |> DataFrame
end

# ============================================================================= STL
"""
STL simplificada (Cleveland et al., 1990), laço interno sem pesos de robustez,
para série mensal completa x (período 12):
 1) destendencia; 2) suaviza cada subsérie de mês (loess, span ns);
 3) passa-baixa (médias móveis 12, 12, 3 + loess) e subtrai → sazonal;
 4) dessazonaliza e suaviza com loess (span nt) → tendência. Repete `iter` vezes.
Lacunas mensais são preenchidas com a climatologia do mês antes da decomposição
(o texto informa quantos meses foram preenchidos).
"""
function stl_simples(x::Vector{Float64}; ns = 0.6, nt = 0.25, iter = 3)
    n = length(x); tt = collect(1.0:n)
    T = zeros(n); Sz = zeros(n)
    mm(v, k) = [mean(v[i:i+k-1]) for i in 1:length(v)-k+1]
    for _ in 1:iter
        d = x .- T
        C = zeros(n + 24)                             # subséries estendidas 1 período de cada lado
        for m in 1:12
            idx = m:12:n; xs = collect(Float64, 1:length(idx))
            fit = loess(xs, d[idx]; span = ns, degree = 1)
            ext = vcat(0.0, xs, length(idx) + 1.0)
            pred = lpredict(fit, clamp.(ext, 1.0, Float64(length(idx))))
            pos = vcat(m, (12 .+ collect(idx)), 12 + n + m)
            for (p, v) in zip(pos, pred)
                p <= n + 24 && (C[p] = v)
            end
        end
        L = mm(mm(mm(C, 12), 12), 3)                  # comprimento n
        L = lpredict(loess(tt, L[1:n]; span = 0.15, degree = 1), tt)
        Sz = C[13:12+n] .- L
        T = lpredict(loess(tt, x .- Sz; span = nt, degree = 1), tt)
    end
    return (tendencia = T, sazonal = Sz, resto = x .- T .- Sz)
end
