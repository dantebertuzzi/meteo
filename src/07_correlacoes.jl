# =============================================================================
# 07_correlacoes.jl — associação entre variáveis
#
# • Pearson (r): associação *linear*; sensível a valores extremos.
# • Spearman (ρ): Pearson sobre postos; capta relações monotônicas e é robusto —
#   preferível para chuva (assimétrica) e para n pequeno.
#   Ambos testados por t = r√[(n−2)/(1−r²)] ~ t_{n−2} sob H₀: correlação nula
#   (para Spearman a aproximação t é adequada a partir de n ≈ 10).
# • Muitas correlações simultâneas → p-valores ajustados por Benjamini–Hochberg.
# • Escala mensal: usam-se ANOMALIAS (remove o ciclo sazonal, que geraria
#   correlação espúria — ex.: UR e chuva sobem juntas no verão por causa da
#   estação, não por relação interanual). Anomalias mensais são autocorrelacionadas;
#   o n efetivo de Bretherton et al. (1999) corrige o teste:
#   n_ef = n(1 − r₁ˣr₁ʸ)/(1 + r₁ˣr₁ʸ).
# =============================================================================

p_correlacao(r, n) = n > 2 && abs(r) < 1 ? 2ccdf(TDist(n - 2), abs(r) * sqrt((n - 2) / (1 - r^2))) : NaN

function par_completo(x, y)
    ok = .!ismissing.(x) .& .!ismissing.(y)
    return Float64.(x[ok]), Float64.(y[ok])
end

"Tabela longa de correlações par a par (Pearson e Spearman) com IC e BH."
function tabela_correlacoes(df, vars)
    linhas = NamedTuple[]
    for i in 1:length(vars)-1, j in i+1:length(vars)
        x, y = par_completo(df[!, vars[i]], df[!, vars[j]]); n = length(x)
        n < 8 && continue
        r = cor(x, y); ρ = corspearman(x, y); li, ls = ic_pearson(r, n)
        push!(linhas, (x = String(vars[i]), y = String(vars[j]), n = n,
                       pearson = r, ic_inf = li, ic_sup = ls, p_pearson = p_correlacao(r, n),
                       spearman = ρ, p_spearman = p_correlacao(ρ, n)))
    end
    t = DataFrame(linhas)
    t.q_pearson = bh_ajuste(t.p_pearson); t.q_spearman = bh_ajuste(t.p_spearman)
    return t
end

"Matriz de correlação (Spearman por padrão) para a figura."
function matriz_correlacao(df, vars; f = corspearman)
    k = length(vars); R = fill(NaN, k, k); N = zeros(Int, k, k)
    for i in 1:k, j in 1:k
        x, y = par_completo(df[!, vars[i]], df[!, vars[j]])
        N[i, j] = length(x); R[i, j] = length(x) >= 8 ? f(x, y) : NaN
    end
    return R, N
end

"Correlação entre anomalias mensais com n efetivo de Bretherton."
function correlacao_mensal(dm, vx, vy; metodo = :pearson)
    ok = .!ismissing.(dm[!, vx]) .& .!ismissing.(dm[!, vy])
    idx = dm.ano[ok] .* 12 .+ dm.mes[ok]
    x = Float64.(dm[ok, vx]); y = Float64.(dm[ok, vy]); n = length(x)
    r = metodo == :pearson ? cor(x, y) : corspearman(x, y)
    r1x = acf_irregular(idx, x, 1)[1]; r1y = acf_irregular(idx, y, 1)[1]
    nef = clamp(n * (1 - r1x * r1y) / (1 + r1x * r1y), 3, n)
    li, ls = ic_pearson(r, nef)
    (x = String(vx), y = String(vy), metodo = String(metodo), n = n, r1x = r1x, r1y = r1y,
     n_efetivo = nef, r = r, ic_inf = li, ic_sup = ls, p_nef = p_correlacao(r, round(Int, nef)))
end

"""
Erros-padrão HAC de Newey–West para MQO com resíduos autocorrelacionados:
V = (X'X)⁻¹ Ŝ (X'X)⁻¹, Ŝ = Σₜ eₜ²xₜxₜ' + Σₗ wₗ Σₜ eₜeₜ₋ₗ(xₜxₜ₋ₗ' + xₜ₋ₗxₜ'),
pesos de Bartlett wₗ = 1 − l/(L+1), L = ⌊4(n/100)^{2/9}⌋.
"""
function newey_west(X, e; L = floor(Int, 4 * (size(X, 1) / 100)^(2 / 9)))
    n = size(X, 1); XtXi = inv(X' * X)
    S = X' * Diagonal(e .^ 2) * X
    for l in 1:L
        w = 1 - l / (L + 1)
        Γ = sum(e[t] * e[t-l] * (X[t, :] * X[t-l, :]') for t in l+1:n)
        S += w * (Γ + Γ')
    end
    return sqrt.(diag(XtXi * S * XtXi)), L
end

"""
Regressão MQO y ~ x₁ + … com diagnósticos e EP HAC (quando `hac=true`).
Devolve tabela de coeficientes, R², VIF (multicolinearidade: VIFⱼ = 1/(1 − R²ⱼ)).
"""
function regressao(df, y::Symbol, xs::Vector{Symbol}; hac = false, nome = "")
    d = dropmissing(df[:, [y; xs]])
    for c in [y; xs]; d[!, c] = Float64.(d[!, c]); end
    f = term(y) ~ sum(term.(xs))
    m = lm(f, d)
    X = modelmatrix(m); e = residuals(m); n = size(X, 1); p = size(X, 2)
    ep = hac ? newey_west(X, e)[1] : stderror(m)
    tq = quantile(TDist(n - p), 0.975)
    vif = map(xs) do xj
        length(xs) == 1 && return 1.0
        outros = setdiff(xs, [xj])
        1 / (1 - r2(lm(term(xj) ~ sum(term.(outros)), d)))
    end
    coefs = DataFrame(modelo = nome, termo = coefnames(m), estimativa = coef(m), ep = ep,
                      ic_inf = coef(m) .- tq .* ep, ic_sup = coef(m) .+ tq .* ep,
                      p = 2 .* ccdf.(TDist(n - p), abs.(coef(m) ./ ep)),
                      vif = [NaN; vif], ep_tipo = hac ? "HAC Newey-West" : "MQO clássico")
    return (m = m, dados = d, coefs = coefs, r2 = r2(m), r2aj = adjr2(m), n = n,
            aic = aic(m), bic = bic(m))
end
