# Testes unitários dos métodos implementados à mão, contra valores publicados
# ou propriedades conhecidas. Uso: julia --project=. test/runtests.jl
using Test
include(joinpath(@__DIR__, "..", "src", "paper", "carrega.jl"))

@testset "Radiação extraterrestre e Hargreaves (FAO-56)" begin
    # FAO-56, Exemplo 8: 20° S, 3 de setembro (J = 246) → Ra = 32,2 MJ m⁻² d⁻¹
    @test radiacao_extraterrestre(-20.0, 246) ≈ 32.2 atol = 0.05
    # Exemplo 8 também dá Ra(20° N, J = 246) — simetria: hemisfério oposto no mesmo dia difere
    @test radiacao_extraterrestre(20.0, 246) > radiacao_extraterrestre(-20.0, 246)
    # Hargreaves: verificação algébrica direta da eq. 52
    @test eto_hargreaves(30.0, 20.0, 32.2) ≈ 0.0023 * 0.408 * 32.2 * 42.8 * sqrt(10) atol = 1e-12
end

@testset "SPI" begin
    rng = MersenneTwister(1)
    anos = repeat(1961:2020, inner = 12); meses = repeat(1:12, 60)
    p = rand(rng, Gamma(2.0, 50.0), length(anos))
    z, aj = spi(p, meses, anos; k = 1, ref = (1961, 2020))
    zz = collect(skipmissing(z))
    @test abs(mean(zz)) < 0.05          # SPI ≈ N(0, 1) quando o modelo é correto
    @test abs(std(zz) - 1) < 0.05
    @test soma_movel([1.0, 2.0, missing, 4.0, 5.0], 2) |> x -> isequal(x, [missing, 3.0, missing, missing, 9.0])
    # massa em zero: com 30 % de zeros, H(0) = q ⇒ SPI(0) = Φ⁻¹(0,3)
    p0 = copy(p); p0[rand(rng, length(p0)) .< 0.3] .= 0.0
    z0, aj0 = spi(p0, meses, anos; k = 1, ref = (1961, 2020))
    i = findfirst(==(0.0), p0)
    @test z0[i] ≈ quantile(Normal(), aj0[meses[i]].q) atol = 1e-8
end

@testset "Log-logística e GEV por PWM" begin
    rng = MersenneTwister(2)
    d = LogLogistica3(50.0, 4.0, -100.0)
    x = [quantil_ll(d, u) for u in rand(rng, 20_000)]
    dh = ajusta_loglogistica(x)
    @test dh.α ≈ 50 rtol = 0.05
    @test dh.β ≈ 4 rtol = 0.05
    @test dh.γ ≈ -100 atol = 3
    @test cdf_ll(d, quantil_ll(d, 0.3)) ≈ 0.3
    # assimetria negativa → ajuste espelhado; F deve ser ≈ uniforme nos dados
    xn = -x
    dn = ajusta_loglogistica(xn)
    @test dn.espelhada
    @test cdf_ll(dn, quantil_ll(dn, 0.3)) ≈ 0.3
    @test mean(cdf_ll.(Ref(dn), xn) .< 0.25) ≈ 0.25 atol = 0.02
    g = GeneralizedExtremeValue(10.0, 3.0, 0.1)
    gh = ajusta_gev_lmom(rand(rng, g, 20_000))
    @test params(gh)[1] ≈ 10 atol = 0.15
    @test params(gh)[2] ≈ 3 atol = 0.15
    @test params(gh)[3] ≈ 0.1 atol = 0.03
end

@testset "Corridas" begin
    dt = Date(2000, 1, 1):Month(1):Date(2000, 10, 1)
    z = [0.0, -1.5, -2.0, 0.0, -1.2, missing, -1.1, -1.3, 0.5, -0.5]
    e = corridas(collect(dt), z; z0 = -1.0)
    @test nrow(e) == 3
    @test e.duracao == [2, 1, 2]
    @test e.severidade[1] ≈ 3.5
    @test e.censurado == [false, true, true]
end

@testset "Período de retorno: perfil e bootstrap" begin
    rng = MersenneTwister(3)
    d = Gamma(4.0, 110.0); x = rand(rng, d, 60); x0 = 150.0
    pr = perfil_p_gamma(x, x0)
    @test pr.inf < pr.p < pr.sup
    @test pr.p ≈ cdf(fit_mle(Gamma, x), x0)
    # cobertura aproximada: o IC de 95 % deve conter o p verdadeiro na maioria das réplicas
    cobre = 0
    for r in 1:60
        xr = rand(rng, d, 40); q = perfil_p_gamma(xr, x0)
        cobre += q.inf <= cdf(d, x0) <= q.sup
    end
    @test cobre / 60 > 0.85
    ps = bootstrap_p(Gamma, x, x0; B = 500)
    @test 0.4 < mean(ps .< pr.p) < 0.7
    @test prob_ao_menos_um(0.1, 10) ≈ 1 - 0.9^10
end

@testset "SNHT" begin
    rng = MersenneTwister(4)
    y = vcat(randn(rng, 15), randn(rng, 15) .+ 2.0)
    r = snht_teste(collect(1.0:30.0), y; B = 2000)
    @test r.p < 0.01
    @test r.ultimo_antes == 15.0
    r0 = snht_teste(collect(1.0:30.0), randn(rng, 30); B = 2000)
    @test r0.p > 0.01
end

@testset "Markov" begin
    anos = 1:8; est = ["S", "S", "N", "N", "S", missing, "S", "N"]
    C = transicoes(anos, est)
    # 1→2 S→S, 2→3 S→N, 3→4 N→N, 4→5 N→S, 7→8 S→N; 5→6 e 6→7 quebram no ausente
    @test C == [1 2; 1 1]
    m = markov_posterior([30 10; 10 50])
    @test m.pSS[2] ≈ 30.5 / 41 atol = 0.01      # média da Beta(30,5; 10,5) ≈ mediana
end

@testset "Referência regional (Alexandersson & Moberg)" begin
    rng = MersenneTwister(5)
    anos = 1961:2020
    clima = Dict(a => 500 * exp(0.3randn(rng)) for a in anos)          # sinal regional comum
    viz = Dict("v$j" => Dict(a => clima[a] * exp(0.1randn(rng)) * (0.8 + 0.1j) for a in anos) for j in 1:4)
    cand = Dict(a => clima[a] * exp(0.1randn(rng)) * (a > 1995 ? 0.75 : 1.0) for a in anos)   # quebra só na candidata
    Q, info = serie_referencia(cand, viz)
    @test all(info.usada)
    r = snht_teste(Float64.(Q.ano), Q.Q; B = 2000)
    @test r.p < 0.01
    @test 1993 <= r.ultimo_antes <= 1997
    # sem quebra na candidata: o teste não deve rejeitar com frequência
    cand2 = Dict(a => clima[a] * exp(0.1randn(rng)) for a in anos)
    Q2, _ = serie_referencia(cand2, viz)
    @test snht_teste(Float64.(Q2.ano), Q2.Q; B = 2000).p > 0.01
end

@testset "P2: sazonalidade, circular, Liebmann, ETCCDI" begin
    @test indice_sazonalidade(fill(10.0, 12)) ≈ 0
    @test indice_sazonalidade([120.0; zeros(11)]) ≈ 11 / 6        # máximo teórico = 1,83
    p = zeros(365); p[200] = 50.0
    c = circular_chuva(p); @test c.concentracao ≈ 1; @test c.centroide ≈ 199.5 atol = 1e-6
    @test circular_chuva(ones(365)).concentracao < 1e-10
    q = vcat(zeros(100), fill(5.0, 120), zeros(145)); L = liebmann(q, sum(q) / 365)
    @test L.inicio == 101          # mínimo de A no dia 100 ⇒ início no 1º dia chuvoso
    @test L.fim == 220
    @test maior_corrida([0, 0, 2, 0, missing, 0, 0, 0], <(1)) == 3
    ic = indices_chuva([0.0, 5.0, 10.0, 0.5, 30.0], 20.0)
    @test ic.R1mm == 3 && ic.Rx1day == 30 && ic.SDII ≈ 15 && ic.R95p == 30 && ic.CWD == 2 && ic.CDD == 1
end

@testset "Tendência: Hamed–Rao e MDT" begin
    rng = MersenneTwister(6)
    t = collect(1.0:50); y = randn(rng, 50)
    hr = mk_hamed_rao(t, y); @test 0 < hr.p <= 1
    # AR(1) de memória curta (ρ = 0,6, n = 300): em média o fator de variância é > 1
    fat = map(1:20) do _
        u = zeros(300); for i in 2:300; u[i] = 0.6u[i-1] + randn(rng); end
        mk_hamed_rao(collect(1.0:300), u).fator
    end
    @test median(fat) > 1.5
    # ruído branco: em geral nenhum lag é significativo ⇒ fator = 1
    @test median([mk_hamed_rao(collect(1.0:300), randn(rng, 300)).fator for _ in 1:20]) ≈ 1 atol = 0.3
    md = tendencia_minima_detectavel(t, y; B = 200)
    @test 0 < md.mdt_decada < 2
end

@testset "STL simplificada" begin
    n = 240; tt = 1:n
    x = 0.01 .* tt .+ 2 .* sin.(2π .* tt ./ 12) .+ 0.1 .* randn(MersenneTwister(7), n)
    s = stl_simples(Float64.(x))
    @test cor(s.sazonal, 2 .* sin.(2π .* tt ./ 12)) > 0.98
    @test cor(s.tendencia, 0.01 .* tt) > 0.95
end

@testset "ESS" begin
    rng = MersenneTwister(8)
    iid = randn(rng, 1000, 4)
    @test 3000 < ess(iid) < 5000                     # amostras independentes: ESS ≈ N
    ar = zeros(1000, 4)
    for j in 1:4, i in 2:1000
        ar[i, j] = 0.9ar[i-1, j] + randn(rng)
    end
    @test ess(ar) < 600                              # AR(1) com ρ = 0,9: ESS ≈ N(1−ρ)/(1+ρ) ≈ 210
end

@testset "Leitura BDMEP diária: chuva > 150 mm é mantida" begin
    d = le_bdmep_diario(ARQ_LONGO)
    i = findfirst(==(Date(1978, 1, 1)), d.data)
    @test d.prec[i] ≈ 151.3            # dia real acima do limite horário de 150 mm
    @test maximum(skipmissing(d.prec)) <= PREC_DIARIA_MAX
end

@testset "EMV da log-normal" begin
    x = [480.0, 300.0, 610.0, 150.0, 420.0, 700.0, 260.0]
    d = ajusta_mle(LogNormal, x)
    @test params(d)[1] ≈ mean(log.(x))
    @test params(d)[2] ≈ sqrt(mean((log.(x) .- mean(log.(x))) .^ 2))     # divisor n
    # o EMV maximiza a verossimilhança: perturbar σ não pode aumentá-la
    @test loglikelihood(d, x) >= loglikelihood(LogNormal(params(d)[1], 1.02params(d)[2]), x)
    @test loglikelihood(d, x) >= loglikelihood(LogNormal(params(d)[1], 0.98params(d)[2]), x)
end

@testset "Não estacionário: candidatas e média de Akaike" begin
    rng = MersenneTwister(3)
    x = collect(range(-3, 3, length = 60)); y = [rand(rng, Gamma(6.0, 500exp(-0.12xi) / 6.0)) for xi in x]
    # a Gamma genérica reproduz o ajuste dedicado já usado no artigo
    g = gamma_nao_estacionaria(y, x); c = ns_candidata(Gamma, y, x)
    @test c.ll ≈ g.ll rtol = 1e-6
    @test c.F(200.0, 1.0) ≈ g.p_em(200.0, 1.0) rtol = 1e-4
    # inclinação fixada em zero reproduz o EMV estacionário de cada família
    for D in CANDIDATAS
        @test ns_candidata(D, y, x; b = 0.0).ll ≈ loglikelihood(ajusta_mle(D, y), y) rtol = 1e-6
    end
    # forma variável: aninha a fixa (verossimilhança não menor) e, com CV que cresce no
    # tempo nos dados simulados, recupera o sinal do efeito na forma
    fv = ns_candidata(Gamma, y, x; forma_var = true)
    @test fv.ll >= c.ll - 1e-6
    z = [rand(rng, Gamma(8exp(-0.3xi), 500 / (8exp(-0.3xi)))) for xi in x for _ in 1:5]
    xz = repeat(x, inner = 5)
    @test ns_candidata(Gamma, z, xz; forma_var = true).c < 0
    # a média é uma combinação convexa das candidatas (4 famílias × forma fixa e variável)
    m = ns_media_akaike(y, x)
    @test length(m.fits) == 2length(CANDIDATAS)
    @test sum(m.w) ≈ 1
    ps = [f.F(200.0, 2.0) for f in m.fits]
    @test minimum(ps) - 1e-12 <= m.p_em(200.0, 2.0) <= maximum(ps) + 1e-12
end

@testset "Poisson-binomial" begin
    # p iguais: reduz à binomial
    @test prob_ao_menos_k(fill(0.2, 10), 3) ≈ ccdf(Binomial(10, 0.2), 2)
    @test prob_ao_menos_k([0.1, 0.5, 0.9], 0) ≈ 1
    @test prob_ao_menos_k([0.1, 0.5, 0.9], 3) ≈ 0.1 * 0.5 * 0.9
    @test prob_ao_menos_k([0.3, 0.6], 1) ≈ 1 - 0.7 * 0.4
end
