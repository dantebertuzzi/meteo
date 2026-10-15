# =============================================================================
# figuras.jl — figuras do artigo (CairoMakie → PDF vetorial no tamanho final)
#
# Convenções: largura = \textwidth (455 pt); fonte New Computer Modern 9 pt;
# paleta validada para daltonismo (azul #2a78d6, laranja #eb6834, verde-água
# #1baf7a; ΔE CVD ≥ 9). Seco = laranja, úmido = azul, neutro = cinza. Toda série
# com cor tem legenda ou rótulo direto; nenhum gráfico usa dois eixos y.
# =============================================================================

const MESES_ABREV = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
const COR_SECO = COR_LARANJA
const COR_UMIDO = COR_AZUL
fmt_br(x; d = 0) = replace(@sprintf("%.*f", d, x), "." => ",", "-" => "−")     # vírgula decimal e sinal de menos tipográfico
ticks_br(v; d = 0) = (v, [fmt_br(x; d) for x in v])

# ----------------------------------------------------------------------------- 1. área de estudo
function fig_area(P)
    est = [("Petrolina 82983 / A307", -40.5233, -9.3886, true)]
    for r in eachrow(P.V.info)
        c = only([k for (k, v) in VIZINHAS if v == r.estacao])
        f = filter(x -> startswith(x, "dados_$(c)_D_"), readdir(joinpath(RAIZ, "data", "raw", "bdmep")))
        isempty(f) && continue
        l = readlines(joinpath(RAIZ, "data", "raw", "bdmep", f[1]))
        lat = parse(Float64, split(l[3], ": ")[2]); lon = parse(Float64, split(l[4], ": ")[2])
        push!(est, (r.estacao * (r.usada ? "" : "*"), lon, lat, false))
    end
    fig = figura(250)
    ax = Axis(fig[1, 1]; xlabel = "Longitude (°)", ylabel = "Latitude (°)", aspect = DataAspect(),
              xtickformat = v -> fmt_br.(v), ytickformat = v -> fmt_br.(v))
    for raio in (100, 200)                                    # km → graus aproximados
        θ = range(0, 2π, length = 200); g = raio / 111.32
        lines!(ax, -40.5233 .+ g .* cos.(θ) ./ cosd(9.39), -9.3886 .+ g .* sin.(θ); color = (COR_CINZA, 0.6), linestyle = :dash, linewidth = 0.6)
        text!(ax, -40.5233, -9.3886 + g; text = "$(raio) km", fontsize = 7, color = COR_TEXTO2, align = (:center, :bottom))
    end
    for (nome, lon, lat, alvo) in est
        scatter!(ax, [lon], [lat]; color = alvo ? COR_LARANJA : COR_AZUL, markersize = alvo ? 9 : 6, strokecolor = :white, strokewidth = 0.8)
        # posição de cada rótulo escolhida à mão para não colidir (direção: :d direita, :e esquerda, :a acima, :b abaixo)
        pos = Dict("Petrolina" => (:b, (0, -6)), "Cabrobó" => (:e, (-5, 2)), "Floresta" => (:a, (0, 5)),
                   "Petrolândia" => (:e, (-5, 2)), "Paulo Afonso" => (:e, (-5, -12)))
        k = findfirst(c -> startswith(nome, c), collect(keys(pos)))
        dir, off = isnothing(k) ? (:d, (5, 2)) : pos[collect(keys(pos))[k]]
        al = dir == :d ? (:left, :bottom) : dir == :e ? (:right, :bottom) : dir == :a ? (:center, :bottom) : (:center, :top)
        text!(ax, lon, lat; text = nome, fontsize = 7, offset = off, align = al, color = COR_TEXTO)
    end
    xlims!(ax, -42.6, -37.9); ylims!(ax, -11.4, -7.4)
    Label(fig[2, 1], "* vizinha não usada na referência regional (ρ < 0,5 ou < 15 anos em comum)"; fontsize = 7, color = COR_TEXTO2, halign = :left)
    salva_pfig(fig, "fig01_area")
end

# ----------------------------------------------------------------------------- 2. ausências (MissingPatterns)
"Percentual de ausentes por ano e variável, calculado com `missingstats` (MissingPatterns.jl)."
function pct_ausente_por_ano(df, vars)
    anos = sort(unique(year.(df.data)))
    M = fill(NaN, length(anos), length(vars))
    for (i, a) in enumerate(anos)
        st = Dict(r.column => r.pct for r in missingstats(select(df[year.(df.data) .== a, :], vars)))
        for (j, v) in enumerate(vars)
            M[i, j] = st[v]
        end
    end
    return anos, M
end

function fig_ausencias(P)
    vars = [:prec, :tmax, :tmin, :ur, :vvel]
    nomes = ["Chuva", "Tmáx", "Tmín", "UR", "Vento"]
    a = CSV.read(joinpath(DIR_PROC, "diario.csv"), DataFrame)
    for v in vars                                           # A307 original: desfaz o preenchimento
        a[!, v] = [p ? missing : x for (x, p) in zip(a[!, v], a[!, Symbol(v, :_pre)])]
    end
    a = a[year.(a.data) .<= 2024, :]
    l = select(P.L.d, :data, vars)
    fig = figura(300)
    cm = cgrad([colorant"#f3f2ef", COR_LARANJA, colorant"#7a2a0c"])
    for (k, (df, tit)) in enumerate(((l, "82983 (convencional), 1961–2025"), (a, "A307 (automática), 2003–2024")))
        anos, M = pct_ausente_por_ano(df, vars)
        ax = Axis(fig[1, k]; title = tit, yticks = (1:length(vars), nomes), xticks = anos[1]:(k == 1 ? 10 : 5):anos[end],
                  xticklabelrotation = 0, ygridvisible = false)
        hm = heatmap!(ax, anos, 1:length(vars), M; colormap = cm, colorrange = (0, 100))
        k == 2 && Colorbar(fig[1, 3], hm; label = "% de dias ausentes", width = 7, ticklabelsize = 8)
    end
    colsize!(fig.layout, 1, Relative(0.58))
    salva_pfig(fig, "fig02_ausencias")
end

# ----------------------------------------------------------------------------- 3. climatologia
function fig_climatologia(P)
    c = P.Z.clim
    fig = figura(220)
    ax1 = Axis(fig[1, 1]; title = "(a) Chuva mensal (mm)", xticks = (1:12, MESES_ABREV), ytickformat = v -> fmt_br.(v))
    rangebars!(ax1, 1:12, c.p10_9120, c.p90_9120; color = (COR_AZUL, 0.45), linewidth = 6, whiskerwidth = 0)
    scatter!(ax1, 1:12, c.prec_9120; color = COR_AZUL, markersize = 7, label = "1991–2020 (média; barra P10–P90)")
    ok = .!ismissing.(c.prec_6190)
    scatter!(ax1, collect(1:12)[ok], Float64.(c.prec_6190[ok]); color = :white, strokecolor = COR_TEXTO, strokewidth = 0.8,
             marker = :diamond, markersize = 6, label = "1961–1990 (média)")
    axislegend(ax1; position = :ct, labelsize = 7)
    ylims!(ax1, 0, nothing)
    ax2 = Axis(fig[1, 2]; title = "(b) Temperatura, 1991–2020 (°C)", xticks = (1:12, MESES_ABREV), ytickformat = v -> fmt_br.(v))
    lines!(ax2, 1:12, c.tmax_9120; color = COR_LARANJA, linewidth = 1.5)
    lines!(ax2, 1:12, c.tmin_9120; color = COR_AZUL, linewidth = 1.5)
    text!(ax2, 12.2, c.tmax_9120[12]; text = "Tmáx", fontsize = 8, color = COR_TEXTO, align = (:left, :center))
    text!(ax2, 12.2, c.tmin_9120[12]; text = "Tmín", fontsize = 8, color = COR_TEXTO, align = (:left, :center))
    xlims!(ax2, 0.5, 13.6)
    salva_pfig(fig, "fig03_climatologia")
end

# ----------------------------------------------------------------------------- 4. chuva anual
function fig_chuva_anual(P)
    dh = P.L.dh; v = dh[dh.valido, :]
    nt = P.N.nt; m = mean(v.ano_hidro)
    fig = figura(210)
    ax = Axis(fig[1, 1]; ylabel = "Chuva no ano hidrológico (mm)", xticks = 1960:10:2020, ytickformat = v -> fmt_br.(v))
    seco = v.prec .<= P.E.x1213 .+ 0.5
    barplot!(ax, v.ano_hidro, Float64.(v.prec); color = [s ? COR_SECO : (COR_CINZA, 0.55) for s in seco], gap = 0.25)
    anos = 1963:2024
    lines!(ax, collect(anos), [exp(nt.a + nt.b * (a - m) / 10) for a in anos]; color = COR_TEXTO, linewidth = 1.3,
           label = "média da Gamma não estacionária")
    if !isnothing(P.N.nf)
        lines!(ax, collect(anos), [exp(P.N.nf.a + P.N.nf.b * (a - m) / 10) for a in anos]; color = COR_TEXTO2, linewidth = 1.0,
               linestyle = :dash, label = "inclinação regional (vizinhas)")
    end
    hlines!(ax, [P.E.x1213]; color = COR_SECO, linewidth = 0.8, linestyle = :dot)
    for r in eachrow(v[seco, :])
        text!(ax, r.ano_hidro, r.prec; text = r.rotulo, fontsize = 6.5, rotation = π / 2, align = (:left, :center), offset = (0, 3), color = COR_TEXTO)
    end
    axislegend(ax; position = :rt, labelsize = 7)
    ylims!(ax, 0, nothing)
    salva_pfig(fig, "fig04_chuva_anual")
end

# ----------------------------------------------------------------------------- 5. SPI/SPEI
function fig_spi_spei(P)
    dm = P.I.dm
    fig = figura(230)
    for (k, (v, tit)) in enumerate(((:spi12, "(a) SPI-12"), (:spei12, "(b) SPEI-12 (Hargreaves, log-logística)")))
        ax = Axis(fig[k, 1]; title = tit, xticks = 1960:10:2020, yticks = -3:1:3, ytickformat = v -> fmt_br.(v))
        t = year.(dm.data) .+ (month.(dm.data) .- 0.5) ./ 12
        z = dm[!, v]; ok = .!ismissing.(z); tt = t[ok]; zz = Float64.(z[ok])
        barplot!(ax, tt, zz; color = [x < 0 ? COR_SECO : COR_UMIDO for x in zz], gap = 0, width = 1 / 12)
        hlines!(ax, [-1.0]; color = COR_TEXTO2, linewidth = 0.6, linestyle = :dash)
        xlims!(ax, 1962, 2026); ylims!(ax, -3.2, 3.2)
        k == 1 && hidexdecorations!(ax; grid = false, ticks = false)
    end
    Label(fig[3, 1], "Barras em laranja: índice < 0 (seco); azul: > 0 (úmido). Tracejado: limiar de evento (−1). Lacunas = meses sem 12 meses válidos.";
          fontsize = 7, color = COR_TEXTO2, halign = :left)
    salva_pfig(fig, "fig05_spi_spei")
end

# ----------------------------------------------------------------------------- 6. período de retorno
function fig_retorno(P)
    x = Float64.(P.E.v.prec); n = length(x); g = fit_mle(Gamma, x)
    Ts = exp.(range(log(1.05), log(500), length = 80)); pnex = 1 ./ Ts
    rng = MersenneTwister(SEMENTE)
    bs = [fit_mle(Gamma, rand(rng, g, n)) for _ in 1:2000]
    Q = hcat([quantile.(Ref(b), pnex) for b in bs]...)
    lo = [quantile(Q[i, :], 0.025) for i in eachindex(Ts)]; hi = [quantile(Q[i, :], 0.975) for i in eachindex(Ts)]
    fig = figura(230)
    ax = Axis(fig[1, 1]; xscale = log10, xlabel = "Período de retorno T (anos) de um ano com chuva ≤ x",
              ylabel = "Chuva no ano hidrológico x (mm)", xticks = ([1, 2, 5, 10, 20, 50, 100, 200, 500], ["1", "2", "5", "10", "20", "50", "100", "200", "500"]),
              ytickformat = v -> fmt_br.(v))
    band!(ax, Ts, lo, hi; color = (COR_AZUL, 0.18))
    lines!(ax, Ts, quantile.(g, pnex); color = COR_AZUL, linewidth = 1.5, label = "Gamma estacionária (IC 95 % bootstrap)")
    nt = P.N.nt; t25 = (2025 - mean(P.E.v.ano_hidro)) / 10
    g25 = Gamma(nt.k, exp(nt.a + nt.b * t25) / nt.k)
    lines!(ax, Ts, quantile.(g25, pnex); color = COR_LARANJA, linewidth = 1.5, label = "Gamma não estacionária, condições de 2025")
    xs = sort(x); Temp = 1 ./ pos_weibull(n)
    scatter!(ax, Temp, xs; color = COR_TEXTO, markersize = 4, label = "observado (posição de Weibull)")
    hlines!(ax, [P.E.x1213]; color = COR_TEXTO2, linewidth = 0.6, linestyle = :dot)
    text!(ax, 1.1, P.E.x1213; text = "2012/13", fontsize = 7, color = COR_TEXTO2, align = (:left, :bottom))
    axislegend(ax; position = :rt, labelsize = 7)
    ylims!(ax, 0, nothing); xlims!(ax, 1, 500)
    salva_pfig(fig, "fig06_periodo_retorno")
end

# ----------------------------------------------------------------------------- 7. posterior de T
function fig_posterior(P)
    pe = P.B.draws["base"]
    po = P.N.po; tc = (2025 - mean(P.E.v.ano_hidro)) / 10
    pn = [cdf(Gamma(exp(po.logk[i]), exp(po.a[i] + po.b[i] * tc) / exp(po.logk[i])), P.E.x1213) for i in eachindex(po.b)]
    fig = figura(190)
    ax = Axis(fig[1, 1]; xlabel = "Período de retorno T (anos, escala log) de um ano com chuva ≤ 2012/13", ylabel = "densidade posterior",
              xticks = (log10.([5, 10, 20, 50, 100, 200, 500, 1000]), ["5", "10", "20", "50", "100", "200", "500", "1000"]))
    for (p, cor, nm) in ((pe, COR_AZUL, "estacionário"), (pn, COR_LARANJA, "não estacionário (2025)"))
        lt = log10.(1 ./ p)
        density!(ax, lt; color = (cor, 0.25), strokecolor = cor, strokewidth = 1.2, label = nm)
        vlines!(ax, [median(lt)]; color = cor, linewidth = 0.8, linestyle = :dash)
    end
    axislegend(ax; position = :rt, labelsize = 7)
    xlims!(ax, log10(3), log10(2000)); ylims!(ax, 0, nothing)
    salva_pfig(fig, "fig07_posterior_T")
end

# ----------------------------------------------------------------------------- 8. Markov
function fig_markov(P)
    t = P.M.tab
    fig = figura(165)
    cm = cgrad([colorant"#f3f2ef", colorant"#8db8ea"])
    for (k, i) in enumerate((1, 3))
        r = t[i, :]
        M = [r.pSS (1 - r.pSS); r.pNS (1 - r.pNS)]
        ax = Axis(fig[1, k]; title = i == 1 ? "(a) Tercil inferior (com tendência)" : "(b) Tercil sem tendência (PIT)",
                  xticks = (1:2, ["seco", "não seco"]), yticks = (1:2, ["seco", "não seco"]), yreversed = true,
                  xlabel = "ano t + 1", ylabel = "ano t", aspect = 1, ygridvisible = false)
        heatmap!(ax, 1:2, 1:2, permutedims(M); colormap = cm, colorrange = (0, 1))
        for a in 1:2, b in 1:2
            ic = a == 1 ? (b == 1 ? (r.pSS_inf, r.pSS_sup) : (1 - r.pSS_sup, 1 - r.pSS_inf)) :
                          (b == 1 ? (r.pNS_inf, r.pNS_sup) : (1 - r.pNS_sup, 1 - r.pNS_inf))
            txt = fmt_br(M[a, b]; d = 2) * "\n[" * fmt_br(ic[1]; d = 2) * "; " * fmt_br(ic[2]; d = 2) * "]"
            text!(ax, b, a; text = txt, align = (:center, :center), fontsize = 7.5, color = COR_TEXTO)
        end
    end
    salva_pfig(fig, "fig08_markov")
end

# ----------------------------------------------------------------------------- 9. tendências (forest)
function fig_tendencias(P)
    t = P.X.tend
    grupos = [("(a) Temperatura (°C/década)", ["Tmáx (°C)", "Tmín (°C)", "Tméd (°C)", "DTR (°C)", "TXx (°C)", "TNn (°C)"], false),
              ("(b) Frequência de extremos (p.p./década)", ["TX90p (% dias)", "TN10p (% dias)"], false),
              ("(c) Chuva: % da média/década", ["R1mm (dias)", "Rx1day (mm)", "SDII (mm/dia)", "R95p (mm)", "CWD (dias)", "CDD (dias)"], true)]
    fig = figura(300)
    for (k, (tit, vs, rel)) in enumerate(grupos)
        ax = Axis(fig[k, 1]; title = tit, yticks = (1:length(vs), [first(split(v, " (")) for v in vs]), yreversed = true,
                  xtickformat = v -> fmt_br.(v; d = rel ? 0 : 1), ygridvisible = false, xgridvisible = true, xgridcolor = (COR_CINZA, 0.25))
        for (i, v) in enumerate(vs)
            r = t[t.variavel .== v, :][1, :]
            esc = rel ? 100 / mean(skipmissing(v in ("CDD (dias)",) ? P.X.anu.CDD_civil :
                        P.X.ich[!, Symbol(replace(first(split(v, " (")), "R1mm" => "R1mm"))])) : 1.0
            sig = r.q_tfpw < 0.05
            cor = sig ? COR_AZUL : COR_CINZA
            rangebars!(ax, [i], [esc * r.sen_ic_inf], [esc * r.sen_ic_sup]; direction = :x, color = cor, linewidth = 1.2, whiskerwidth = 4)
            scatter!(ax, [esc * r.sen_decada], [i]; color = sig ? cor : :white, strokecolor = cor, strokewidth = 1.0, markersize = 7)
        end
        vlines!(ax, [0.0]; color = COR_TEXTO2, linewidth = 0.6)
    end
    Label(fig[4, 1], "Ponto: declive de Sen; barra: IC 95 %. Cheio (azul): q < 0,05 (MK-TFPW, Benjamini–Hochberg); vazio (cinza): não significativo.";
          fontsize = 7, color = COR_TEXTO2, halign = :left)
    salva_pfig(fig, "fig09_tendencias")
end

# ----------------------------------------------------------------------------- 10. estação chuvosa
function fig_estacao(P)
    s = P.Z.s
    fig = figura(200)
    ax = Axis(fig[1, 1]; ylabel = "Data (ano hidrológico)", xticks = 1960:10:2020,
              yticks = ([1, 32, 63, 93, 124, 154, 185, 216, 244, 275, 305, 336],
                        ["1/jul", "1/ago", "1/set", "1/out", "1/nov", "1/dez", "1/jan", "1/fev", "1/mar", "1/abr", "1/mai", "1/jun"]))
    s = s[s.estacao_definida, :]
    rangebars!(ax, s.ano_hidro, Float64.(s.inicio), Float64.(s.fim); color = (COR_CINZA, 0.5), linewidth = 1.5, whiskerwidth = 0)
    for (c, cor, nm) in ((:inicio, COR_AZUL, "início"), (:fim, COR_LARANJA, "fim"))
        scatter!(ax, s.ano_hidro, Float64.(s[!, c]); color = cor, markersize = 5, label = nm)
        ts = theil_sen(Float64.(s.ano_hidro), Float64.(s[!, c]))
        xx = [minimum(s.ano_hidro), maximum(s.ano_hidro)]
        lines!(ax, xx, ts.intercepto .+ ts.declive .* xx; color = cor, linewidth = 1.2)
    end
    Legend(fig[2, 1], ax; orientation = :horizontal, labelsize = 7, tellheight = true, padding = (0, 0, 0, 0))
    rowgap!(fig.layout, 2)
    Label(fig[3, 1], "Método de Liebmann et al. (2007); retas: declive de Sen; barras cinza: duração. Anos sem estação definida omitidos.";
          fontsize = 7, color = COR_TEXTO2, halign = :left)
    salva_pfig(fig, "fig10_estacao_chuvosa")
end

# ----------------------------------------------------------------------------- 11. quantílica
function fig_quantilica(P)
    q = P.QR
    fig = figura(170)
    ax = Axis(fig[1, 1]; xlabel = "Quantil τ", ylabel = "Tendência (°C/década)", xticks = ([0.1, 0.5, 0.9], ["0,1", "0,5", "0,9"]),
              ytickformat = v -> fmt_br.(v; d = 2))
    for (v, cor, nm, dx) in (("tmax", COR_LARANJA, "Tmáx", -0.015), ("tmin", COR_AZUL, "Tmín", 0.015))
        r = q[q.variavel .== v, :]
        rangebars!(ax, r.tau .+ dx, r.ic_inf, r.ic_sup; color = cor, linewidth = 1.2, whiskerwidth = 5)
        scatterlines!(ax, r.tau .+ dx, r.beta_decada; color = cor, markersize = 7, linewidth = 1.0, label = nm)
    end
    hlines!(ax, [0.0]; color = COR_TEXTO2, linewidth = 0.6)
    axislegend(ax; position = :rb, labelsize = 7, orientation = :horizontal)
    xlims!(ax, 0, 1)
    salva_pfig(fig, "fig11_quantilica")
end


# ----------------------------------------------------------------------------- 12. padrões de ausência (MissingPatterns)
"""
Padrões de ausência conjunta (à la mice::md.pattern), a partir de
`missingpatternstats` do MissingPatterns.jl: cada linha é uma combinação de
variáveis ausentes no mesmo dia; a barra mostra a fração de dias com o padrão.
"""
function fig_padroes(P)
    vars = [:prec, :tmax, :tmin, :ur, :vvel]; nomes = ["Chuva", "Tmáx", "Tmín", "UR", "Vento"]
    a = CSV.read(joinpath(DIR_PROC, "diario.csv"), DataFrame)
    for v in vars
        a[!, v] = [p ? missing : x for (x, p) in zip(a[!, v], a[!, Symbol(v, :_pre)])]
    end
    a = a[year.(a.data) .<= 2024, :]
    l = P.L.d[P.L.d.ano .>= INICIO_TEMP, :]
    fig = figura(190)
    for (k, (df, tit)) in enumerate(((l, "(a) 82983, 1971–2025"), (a, "(b) A307, 2003–2024")))
        ps = sort(missingpatternstats(select(df, vars)); by = r -> -r.n)[1:min(6, end)]
        np = length(ps)
        ax = Axis(fig[1, 2k-1]; title = tit, xticks = (1:length(vars), nomes), yticks = (1:np, ["" for _ in 1:np]),
                  yreversed = true, xgridvisible = false, ygridvisible = false, xticklabelsize = 8)
        for (i, r) in enumerate(ps), (j, v) in enumerate(vars)
            aus = getproperty(r.pattern, v)
            scatter!(ax, [j], [i]; marker = :rect, markersize = 11, color = aus ? COR_LARANJA : (COR_CINZA, 0.25))
        end
        ax2 = Axis(fig[1, 2k]; yticks = (1:np, ["" for _ in 1:np]), yreversed = true, xgridvisible = true,
                   ygridvisible = false, xlabel = "% dos dias", xtickformat = v -> fmt_br.(v))
        barh = [r.pct for r in ps]
        barplot!(ax2, 1:np, barh; direction = :x, color = COR_AZUL, gap = 0.35)
        for (i, b) in enumerate(barh)
            text!(ax2, b, i; text = fmt_br(b; d = 1), offset = (3, 0), align = (:left, :center), fontsize = 7, color = COR_TEXTO)
        end
        xlims!(ax2, 0, 115); hideydecorations!(ax2); linkyaxes!(ax, ax2)
        colsize!(fig.layout, 2k, Relative(0.17))
    end
    Label(fig[2, :], "Quadrado laranja: variável ausente no padrão; cinza: presente. Seis padrões mais frequentes (missingpatternstats, MissingPatterns.jl).";
          fontsize = 7, color = COR_TEXTO2, halign = :left)
    salva_pfig(fig, "fig12_padroes_ausencia")
end

function gera_figuras_paper(P)
    titulo("Figuras do artigo")
    set_theme!(tema_artigo())
    for f in (fig_area, fig_ausencias, fig_climatologia, fig_chuva_anual, fig_spi_spei, fig_retorno, fig_posterior,
              fig_markov, fig_tendencias, fig_estacao, fig_quantilica, fig_padroes)
        try
            f(P)
        catch e
            @warn "Falha em $(nameof(f))" exception = (e, catch_backtrace())
        end
    end
end
