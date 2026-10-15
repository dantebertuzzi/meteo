# =============================================================================
# homogeneidade.jl — homogeneidade relativa entre A307 e 82983
#
# Sem uma estação de referência independente, testa-se a série de diferenças
# (temperatura) ou razões (chuva) entre as duas estações do mesmo sítio no período
# comum. O clima regional se cancela; o que sobra é diferença instrumental.
# Uma quebra indica que uma das séries mudou (troca de sensor, observador,
# entupimento do pluviômetro) — não diz qual.
#
# SNHT (Alexandersson, 1986): zᵢ = (qᵢ − q̄)/s;
#   T(k) = k·z̄₁² + (n − k)·z̄₂²,   T₀ = max_k T(k), k = 1…n−1.
# Os valores críticos dependem de n; o p-valor é obtido por Monte Carlo sob H₀
# (série i.i.d. normal do mesmo tamanho), o que evita tabelas aproximadas.
# Pettitt (já implementado em 04_temporal.jl) é o complemento não paramétrico.
# =============================================================================

function snht(y)
    n = length(y); z = (y .- mean(y)) ./ std(y)
    T = [k * mean(z[1:k])^2 + (n - k) * mean(z[k+1:n])^2 for k in 1:n-1]
    T0, k = findmax(T)
    return (T0 = T0, k = k, T = T)
end

function snht_teste(t, y; B = 20_000, rng = MersenneTwister(SEMENTE))
    o = sortperm(t); t = t[o]; y = y[o]; n = length(y)
    r = snht(y)
    nulos = [snht(randn(rng, n)).T0 for _ in 1:B]
    p = (count(>=(r.T0), nulos) + 1) / (B + 1)
    return (T0 = r.T0, p = p, T95 = quantile(nulos, 0.95), ultimo_antes = t[r.k],
            media_antes = mean(y[1:r.k]), media_depois = mean(y[r.k+1:end]), n = n)
end

# ----------------------------------------------------------------------------- referência regional
"""
    serie_referencia(cand, vizinhas; min_comum = 15, rmin = 0.5)

Série de razões Qᵢ = Yᵢ / Rᵢ (Alexandersson & Moberg, 1997) para chuva anual:
    Rᵢ = Σⱼ ρⱼ² · Xⱼᵢ·(Ȳ/X̄ⱼ) / Σⱼ ρⱼ²,
com ρⱼ a correlação entre as *primeiras diferenças* da candidata e da vizinha j
(as diferenças removem tendências comuns e realçam a covariação ano a ano).
Só entram vizinhas com ≥ `min_comum` anos em comum e ρⱼ ≥ `rmin`; em cada ano,
Rᵢ usa as vizinhas válidas naquele ano (pesos renormalizados).
`cand` e cada vizinha: Dict(ano => total).
"""
function serie_referencia(cand::Dict, vizinhas::Dict; min_comum = 15, rmin = 0.5)
    pesos = Dict{String, Float64}(); info = NamedTuple[]
    for (nome, viz) in vizinhas
        comuns = sort([a for a in keys(cand) if haskey(viz, a)])
        dif = [(cand[a+1] - cand[a], viz[a+1] - viz[a]) for a in comuns if haskey(cand, a + 1) && haskey(viz, a + 1)]
        ρ = length(dif) >= 8 ? cor(first.(dif), last.(dif)) : NaN
        usa = length(comuns) >= min_comum && isfinite(ρ) && ρ >= rmin
        push!(info, (estacao = nome, anos_comuns = length(comuns), pares_dif = length(dif), rho_dif = ρ, usada = usa))
        usa && (pesos[nome] = ρ^2)
    end
    anos = sort(collect(keys(cand)))
    escala = Dict(n => mean(cand[a] for a in anos if haskey(vizinhas[n], a)) /
                      mean(vizinhas[n][a] for a in anos if haskey(vizinhas[n], a)) for n in keys(pesos))
    linhas = NamedTuple[]
    for a in anos
        ns = [n for n in keys(pesos) if haskey(vizinhas[n], a)]
        isempty(ns) && continue
        R = sum(pesos[n] * vizinhas[n][a] * escala[n] for n in ns) / sum(pesos[n] for n in ns)
        push!(linhas, (ano = a, candidata = cand[a], referencia = R, Q = cand[a] / R, n_viz = length(ns)))
    end
    return DataFrame(linhas), DataFrame(info)
end
