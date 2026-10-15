# =============================================================================
# 03_descritiva.jl — estatística descritiva clássica e robusta + utilitários
# =============================================================================

"""
    descreve(x; nome, razao=true)

Medidas de posição, dispersão e forma.
- variância amostral s² = Σ(xᵢ − x̄)²/(n − 1)
- CV = s/x̄ — só faz sentido em escala de razão (zero absoluto). Para °C (escala
  intervalar) o CV depende da origem arbitrária da escala e é omitido (`razao=false`).
- assimetria g₁ = m₃/m₂^{3/2}; EP aproximado √[6n(n−1)/((n−2)(n+1)(n+3))]
- curtose em excesso g₂ = m₄/m₂² − 3 (0 para a Normal)
- robustas: MAD normalizado (1,4826·mediana|xᵢ − mediana|, consistente com σ
  sob normalidade), média aparada a 10 % e assimetria de Bowley
  (Q₃ + Q₁ − 2·Q₂)/(Q₃ − Q₁), que não depende das caudas.
"""
function descreve(x; nome = "", razao = true)
    x = collect(Float64, skipmissing(x)); n = length(x)
    q1, q2, q3 = quantile(x, [0.25, 0.5, 0.75])
    s = std(x)
    ep_g1 = n > 3 ? sqrt(6n * (n - 1) / ((n - 2) * (n + 1) * (n + 3))) : NaN
    (variavel = nome, n = n, media = mean(x), mediana = q2, minimo = minimum(x),
     maximo = maximum(x), amplitude = maximum(x) - minimum(x), variancia = var(x),
     dp = s, cv_pct = razao ? 100s / mean(x) : NaN, q1 = q1, q3 = q3, iqr = q3 - q1,
     assimetria = skewness(x), ep_assimetria = ep_g1, curtose_excesso = kurtosis(x),
     mad = mad(x; normalize = true), media_aparada10 = mean(trim(x; prop = 0.1)),
     assim_bowley = (q3 + q1 - 2q2) / (q3 - q1))
end

tabela_descritiva(df, vars; razao = Dict()) =
    DataFrame([descreve(df[!, v]; nome = String(v), razao = get(razao, v, true)) for v in vars])

"Ajuste de Benjamini–Hochberg (controle da taxa de falsas descobertas)."
function bh_ajuste(p::AbstractVector)
    n = length(p); o = sortperm(p); q = similar(p, Float64)
    acc = 1.0
    for k in n:-1:1
        acc = min(acc, p[o[k]] * n / k)
        q[o[k]] = acc
    end
    return q
end

"Intervalo percentil bootstrap de `estat(x)`."
function ic_bootstrap(estat, x; B = NBOOT, nivel = 0.95, rng = MersenneTwister(SEMENTE))
    x = collect(skipmissing(x)); n = length(x)
    bs = [estat(x[rand(rng, 1:n, n)]) for _ in 1:B]
    a = (1 - nivel) / 2
    return quantile(bs, [a, 1 - a])
end

"IC de Pearson por transformação z de Fisher: z = atanh(r), EP = 1/√(n−3)."
function ic_pearson(r, n; nivel = 0.95)
    z = atanh(r); e = quantile(Normal(), 1 - (1 - nivel) / 2) / sqrt(n - 3)
    return tanh(z - e), tanh(z + e)
end

"Autocorrelação amostral até `maxlag`, usando apenas pares realmente defasados de k no tempo."
function acf_irregular(t, y, maxlag)
    μ = mean(y); v = sum(abs2, y .- μ) / length(y)
    pos = Dict(round(Int, ti) => yi for (ti, yi) in zip(t, y))
    map(1:maxlag) do k
        prs = [(pos[ti], pos[ti + k]) for ti in keys(pos) if haskey(pos, ti + k)]
        isempty(prs) ? NaN : sum((a - μ) * (b - μ) for (a, b) in prs) / length(prs) / v
    end
end

salva_tabela(df, nome) = CSV.write(joinpath(DIR_TAB, nome * ".csv"), df)

"Imprime DataFrame com números arredondados (apenas exibição)."
function mostra(df; digitos = 3, titulo = "")
    isempty(titulo) || println("\n", "="^length(titulo), "\n", titulo, "\n", "="^length(titulo))
    arred(x) = x isa AbstractFloat ? round(x; sigdigits = digitos + 2) : x
    show(stdout, MIME("text/plain"), mapcols(c -> arred.(c), df); allrows = true, allcols = true,
         summary = false, eltypes = false)
    println()
end
