# =============================================================================
# markov.jl — persistência de anos secos: cadeia de Markov de 1ª ordem
#
# Estados anuais sₜ ∈ {seco, não seco} (seco: SPI do ano hidrológico ≤ limiar).
# Só contam transições entre anos *consecutivos* válidos (lacunas quebram a
# cadeia). Com contagens nᵢⱼ e priori de Jeffreys Dirichlet(½, ½) por linha, a
# posterior de cada linha da matriz é Dirichlet(½ + nᵢ₁, ½ + nᵢ₂) (conjugada).
#
# Quantidades derivadas (por draw da posterior):
#   a = P(N→S), b = P(S→N); distribuição estacionária π_S = a/(a + b);
#   tempo médio de recorrência ao estado seco = 1/π_S (anos);
#   duração média de uma sequência seca = 1/b;
#   espera média, partindo de um ano não seco, até o próximo seco = 1/a.
# Independência (sem persistência) ⇔ P(S|S) = P(S|N).
# =============================================================================

const ESTADOS = ("S", "N")

"Contagens de transição entre anos consecutivos. `estado[t]` ∈ {\"S\",\"N\",missing}."
function transicoes(anos, estado)
    pos = Dict(a => e for (a, e) in zip(anos, estado) if !ismissing(e))
    C = zeros(Int, 2, 2)
    for (a, e) in pos
        haskey(pos, a + 1) || continue
        C[e == "S" ? 1 : 2, pos[a+1] == "S" ? 1 : 2] += 1
    end
    return C
end

function markov_posterior(C; α0 = 0.5, B = 20_000, rng = MersenneTwister(SEMENTE))
    dS = Dirichlet(α0 .+ C[1, :]); dN = Dirichlet(α0 .+ C[2, :])
    pSS = [rand(rng, dS)[1] for _ in 1:B]
    pNS = [rand(rng, dN)[1] for _ in 1:B]
    b = 1 .- pSS; a = pNS
    πS = a ./ (a .+ b)
    q(v) = quantile(v, [0.025, 0.5, 0.975])
    return (pSS = q(pSS), pNS = q(pNS), dif = q(pSS .- pNS), prob_dif_pos = mean(pSS .> pNS),
            piS = q(πS), recorrencia = q(1 ./ πS), duracao_seca = q(1 ./ b), espera = q(1 ./ a),
            draws = (pSS = pSS, pNS = pNS))
end
