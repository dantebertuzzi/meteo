# =============================================================================
# 09_figuras.jl — figuras para publicação (CairoMakie)
#
# Convenções: largura de página dupla = 504 pt (7 pol.), simples = 252 pt;
# PDF vetorial + PNG a 300 dpi. Paleta categórica validada para daltonismo
# (azul, laranja, aqua — os três primeiros slots passam todas as comparações par a par);
# divergente azul ↔ cinza ↔ vermelho para grandezas com sinal (SPI, correlação).
# Vírgula decimal nos rótulos (texto em português).
# =============================================================================

fmtnum(x, d = 2) = replace(replace(@sprintf("%.*f", d, x), "." => ","), "-" => "−")
fmt_ticks(vs) = [isinteger(v) && abs(v) < 1e5 ? replace(string(Int(v)), "-" => "−") :
                 replace(replace(string(round(v; sigdigits = 4)), "." => ","), "-" => "−") for v in vs]
fmtp(p) = p < 0.001 ? "p < 0,001" : "p = " * fmtnum(p, 3)
nanify(x) = [ismissing(v) ? NaN : Float64(v) for v in x]

const TEMA_ARTIGO = Theme(
    fontsize = 10, size = (504, 300), backgroundcolor = :white,
    Axis = (xgridvisible = false, ygridcolor = (:black, 0.07), topspinevisible = false,
            rightspinevisible = false, spinewidth = 0.7, xtickwidth = 0.7, ytickwidth = 0.7,
            xticksize = 3, yticksize = 3, titlealign = :left, titlesize = 10, titlefont = :bold,
            titlecolor = COR_TEXTO, xlabelcolor = COR_TEXTO2, ylabelcolor = COR_TEXTO2,
            xticklabelcolor = COR_TEXTO2, yticklabelcolor = COR_TEXTO2, xticklabelsize = 8.5,
            yticklabelsize = 8.5, xtickformat = fmt_ticks, ytickformat = fmt_ticks,
            subtitlecolor = COR_TEXTO2, subtitlesize = 8.5),
    Legend = (framevisible = false, labelsize = 8.5, patchsize = (14, 8), labelcolor = COR_TEXTO2),
    Colorbar = (ticklabelsize = 8.5, labelsize = 9, spinewidth = 0.5, tickformat = fmt_ticks),
    Lines = (linewidth = 1.6,), Scatter = (markersize = 7,))

function salva(fig, nome)
    save(joinpath(DIR_FIG, nome * ".pdf"), fig; pt_per_unit = 1)
    save(joinpath(DIR_FIG, nome * ".png"), fig; px_per_unit = 300 / 72)
    println("  figura: ", nome)
end

const SEQ_AZUL = cgrad([colorant"#fcfcfb", colorant"#cde2fb", colorant"#6da7ec", colorant"#2a78d6", colorant"#104281"])
const DIV_AZ_VM = cgrad([colorant"#104281", colorant"#6da7ec", COR_NEUTRO, colorant"#f08f8e", colorant"#b8201f"])

# ----------------------------------------------------------------------------- 01 completude
function fig_completude(R)
    fig = Figure(size = (504, 330))
    for (k, (tab, tit)) in enumerate([(R[:completude_prec], "Precipitação"), (R[:completude_t], "Temperatura do ar")])
        anos = sort(unique(tab.ano))
        M = [tab.pct_ausente[(tab.ano .== a) .& (tab.mes .== m)][1] for m in 1:12, a in anos]
        ax = Axis(fig[1, k]; title = tit, xtickformat = Makie.automatic, xticks = (1:12, ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]),
                  ytickformat = Makie.automatic, yticks = (anos[1:3:end], string.(anos[1:3:end])), yreversed = true, ylabel = k == 1 ? "Ano" : "")
        hm = heatmap!(ax, 1:12, anos, M; colormap = SEQ_AZUL, colorrange = (0, 100))
        k == 2 && Colorbar(fig[1, 3], hm; label = "Horas ausentes (%)", width = 8)
    end
    Label(fig[0, :], "Falhas de registro concentram-se na estação chuvosa (nov–abr) de 2010, 2013, 2018 e 2022–24";
          fontsize = 9.5, halign = :left, color = COR_TEXTO2)
    salva(fig, "fig01_completude")
end

# ----------------------------------------------------------------------------- 02 radiação
function fig_radiacao(R)
    r = filter(:n_pos => >(0), R[:radiacao])
    fig = Figure(size = (504, 250))
    ax = Axis(fig[1, 1]; yscale = log10, xlabel = "Ano", ylabel = "Radiação horária (kJ m⁻²)",
              title = "Radiação global: escala inconsistente em 2005–08 e queda de nível em 2014",
              xticks = 2003:3:2024, ytickformat = Makie.automatic, yticks = ([1e3, 3e3, 1e4, 3e4], ["1.000", "3.000", "10.000", "30.000"]))
    hlines!(ax, [LIMITES[:rad][2]]; color = COR_VERM, linestyle = :dash, linewidth = 1)
    text!(ax, 2016.5, LIMITES[:rad][2]; text = "teto físico ≈ 4.900", color = COR_VERM, fontsize = 8, align = (:left, :bottom))
    lines!(ax, r.ano, r.p95; color = COR_AZUL); scatter!(ax, r.ano, r.p95; color = COR_AZUL, strokecolor = :white, strokewidth = 0.8)
    lines!(ax, r.ano, r.p50; color = COR_LARANJA); scatter!(ax, r.ano, r.p50; color = COR_LARANJA, strokecolor = :white, strokewidth = 0.8)
    axislegend(ax, [LineElement(color = COR_AZUL), LineElement(color = COR_LARANJA)], ["percentil 95", "mediana"];
               position = :rt, orientation = :horizontal)
    salva(fig, "fig02_radiacao_inconsistencia")
end

# ----------------------------------------------------------------------------- 03 chuva anual
function cor_spi(z)
    ismissing(z) && return COR_CINZA
    z <= -1 ? COR_VERM : z >= 1 ? COR_AZUL : colorant"#b9b8b3"
end

function fig_prec_anual(R)
    dh = R[:dh]
    v = dh[(dh.ano_hidro .>= 2003) .& (dh.ano_hidro .<= 2023), :]
    x = v.ano_hidro; ok = .!ismissing.(v.prec)
    fig = Figure(size = (504, 290))
    ax = Axis(fig[1, 1]; ylabel = "Precipitação (mm ano⁻¹)", xlabel = "Ano hidrológico (jul–jun)",
              xtickformat = Makie.automatic, xticks = (x, [@sprintf("%02d/%02d", a % 100, (a + 1) % 100) for a in x]), xticklabelrotation = π / 2,
              title = "Precipitação anual e anos anômalos (SPI)")
    barplot!(ax, x[ok], Float64.(v.prec[ok]); color = cor_spi.(v.spi[ok]), gap = 0.25, strokewidth = 0)
    μ = mean(skipmissing(v.prec))
    hlines!(ax, [μ]; color = COR_TEXTO2, linestyle = :dash, linewidth = 0.9)
    text!(ax, x[end] + 0.6, μ; text = "média\n" * fmtnum(μ, 0) * " mm", fontsize = 7.5, color = COR_TEXTO2, align = (:left, :center))
    scatter!(ax, x[.!ok], fill(12.0, count(.!ok)); marker = :xcross, color = COR_CINZA, markersize = 8)
    mm = media_movel(x, v.prec)
    lines!(ax, x, nanify(mm); color = COR_TEXTO, linewidth = 1.3)
    for (xi, z) in zip(x[ok], v.spi[ok])
        abs(z) >= 1.5 && text!(ax, xi, Float64(v.prec[x .== xi][1]); text = "SPI " * fmtnum(z, 1),
                                fontsize = 7.5, align = (:center, :bottom), offset = (0, 2), color = COR_TEXTO)
    end
    xlims!(ax, x[1] - 0.8, x[end] + 2.6); ylims!(ax, 0, 1.12maximum(skipmissing(v.prec)))
    Legend(fig[2, 1], [PolyElement(color = COR_VERM), PolyElement(color = colorant"#b9b8b3"), PolyElement(color = COR_AZUL),
                       LineElement(color = COR_TEXTO), MarkerElement(marker = :xcross, color = COR_CINZA)],
           ["SPI ≤ −1 (seco)", "−1 < SPI < 1", "SPI ≥ 1 (chuvoso)", "média móvel 5 anos", "excluído (FCA > 5 %)"];
           orientation = :horizontal, nbanks = 2, tellheight = true)
    salva(fig, "fig03_precipitacao_anual")
end

# ----------------------------------------------------------------------------- 04 temperatura anual
function fig_temp_anual(R)
    a = R[:anual]; a = a[a.ano .<= 2024, :]
    fig = Figure(size = (504, 420))
    specs = [(:tmax, "Média das máximas diárias", COR_LARANJA), (:tmed, "Média diária", COR_TEXTO2),
             (:tmin, "Média das mínimas diárias", COR_AZUL)]
    for (k, (v, tit, c)) in enumerate(specs)
        ax = Axis(fig[k, 1]; title = tit, ylabel = "°C", xticks = 2003:3:2024,
                  xticklabelsvisible = k == 3, xlabel = k == 3 ? "Ano" : "")
        y = nanify(a[!, v])
        lines!(ax, a.ano, y; color = (c, 0.5), linewidth = 1)
        scatter!(ax, a.ano, y; color = c, strokecolor = :white, strokewidth = 0.8)
        lines!(ax, a.ano, nanify(media_movel(a.ano, a[!, v])); color = COR_TEXTO, linewidth = 1.3)
        xlims!(ax, 2002.5, 2024.5)
    end
    Legend(fig[4, 1], [LineElement(color = COR_TEXTO)], ["média móvel centrada de 5 anos (≥ 3 anos válidos)"];
           orientation = :horizontal)
    rowgap!(fig.layout, 6)
    salva(fig, "fig04_temperatura_anual")
end

# ----------------------------------------------------------------------------- 05 tendências
function painel_tendencia!(ax, mt, tend_row; cor = COR_AZUL, unidade = "")
    d = mt.dados
    tg = collect(range(minimum(d.ano), maximum(d.ano), length = 100))
    ŷ, lo, hi = faixa_confianca(mt, tg)
    band!(ax, tg, lo, hi; color = (cor, 0.18))
    lines!(ax, tg, ŷ; color = cor)
    ts = theil_sen(d.ano, d.y)
    lines!(ax, tg, ts.intercepto .+ ts.declive .* tg; color = COR_TEXTO, linestyle = :dash, linewidth = 1.1)
    scatter!(ax, d.ano, d.y; color = cor, strokecolor = :white, strokewidth = 0.8, markersize = 6)
    r = mt.resumo
    ax.subtitle = "MQO " * fmtnum(r.beta1_decada, unidade == "mm" ? 0 : 2) * " " * unidade * "/déc. (" * fmtp(r.p_beta1) *
                  ");  MK " * fmtp(tend_row.p_mk)
end

function fig_tendencias(R)
    mods, tend = R[:mods], R[:tend]
    lin(nome) = tend[startswith.(tend.variavel, nome * " "), :][1, :]
    specs = [("prec", "Precipitação (ano hidrol.)", "mm", COR_AZUL), ("tmed", "Temperatura média", "°C", COR_TEXTO2),
             ("tmax", "Temperatura máxima", "°C", COR_LARANJA), ("tmin", "Temperatura mínima", "°C", COR_AZUL),
             ("ur", "Umidade relativa", "%", COR_AQUA), ("dtr", "Amplitude térmica diária", "°C", COR_LARANJA)]
    fig = Figure(size = (504, 470))
    for (k, (v, tit, u, c)) in enumerate(specs)
        i, j = fldmod1(k, 2)
        ax = Axis(fig[i, j]; title = tit, ylabel = u, xticks = 2004:5:2024)
        painel_tendencia!(ax, mods[v], lin(v); cor = c, unidade = u)
    end
    Legend(fig[4, 1:2], [[PolyElement(color = (COR_CINZA, 0.25)), LineElement(color = COR_CINZA)],
                         LineElement(color = COR_TEXTO, linestyle = :dash)],
           ["reta MQO e IC 95 % da média", "Theil-Sen"]; orientation = :horizontal)
    salva(fig, "fig05_tendencias_ic")
end

# ----------------------------------------------------------------------------- 06 anomalias mensais
function fig_anomalias_mensais(R)
    dm = R[:dm]; tm = R[:tend_m]
    fig = Figure(size = (504, 360))
    for (k, (v, tit, c)) in enumerate([(:tmax_anom, "Anomalia mensal da temperatura máxima", COR_LARANJA),
                                       (:tmin_anom, "Anomalia mensal da temperatura mínima", COR_AZUL)])
        ok = .!ismissing.(dm[!, v]); t = dm.tempo[ok]; y = Float64.(dm[ok, v])
        linha = tm[tm.variavel .== String(v), :][1, :]
        ax = Axis(fig[k, 1]; title = tit, ylabel = "°C", xticks = 2003:3:2024,
                  subtitle = "tendência " * fmtnum(linha.tend_decada, 2) * " °C/década [IC 95 % HAC " *
                             fmtnum(linha.ic_inf, 2) * "; " * fmtnum(linha.ic_sup, 2) * "], " * fmtp(linha.p_hac))
        hlines!(ax, [0]; color = COR_TEXTO2, linewidth = 0.6)
        vspan!(ax, [2013.9], [2014.1]; color = (COR_CINZA, 0.25))
        scatter!(ax, t, y; color = (c, 0.55), markersize = 3.5)
        b = [ones(length(t)) t] \ y
        lines!(ax, [minimum(t), maximum(t)], b[1] .+ b[2] .* [minimum(t), maximum(t)]; color = COR_TEXTO, linewidth = 1.3)
        k == 1 && text!(ax, 2014.2, maximum(y); text = "possível troca instrumental\n(quebra na radiação)", fontsize = 7,
                        color = COR_TEXTO2, align = (:left, :top))
        xlims!(ax, 2002.8, 2025)
    end
    salva(fig, "fig06_anomalias_mensais_temperatura")
end

# ----------------------------------------------------------------------------- 07–09 distribuições
const CORES_DIST = Dict("Gamma" => COR_AZUL, "Weibull" => COR_LARANJA, "LogNormal" => COR_AQUA, "Normal" => colorant"#eda100")

function fig_distribuicoes(R)
    x = collect(Float64, skipmissing(R[:dh].prec)); dist = R[:dist]
    xs = range(0, 1.25maximum(x), length = 300)
    rot(r) = r.distribuicao * " (ΔAIC " * fmtnum(r.dAIC, 1) * ")"

    fig = Figure(size = (504, 280))
    ax = Axis(fig[1, 1]; xlabel = "Precipitação anual (mm)", ylabel = "Densidade",
              title = "Histograma e densidades ajustadas (MV)", subtitle = "n = $(length(x)) anos hidrológicos válidos")
    hist!(ax, x; bins = 0:125:1000, normalization = :pdf, color = (COR_CINZA, 0.3), strokecolor = :white, strokewidth = 1)
    for r in eachrow(dist)
        lines!(ax, xs, pdf.(r.modelo, xs); color = CORES_DIST[r.distribuicao], label = rot(r))
    end
    scatter!(ax, x, zeros(length(x)); marker = :vline, markersize = 8, color = COR_TEXTO)
    ax2 = Axis(fig[1, 2]; xlabel = "Precipitação anual (mm)", ylabel = "F(x)", title = "ECDF (preto) × modelos ajustados")
    ecdfplot!(ax2, x; color = COR_TEXTO, linewidth = 1.3, label = "empírica")
    for r in eachrow(dist)
        lines!(ax2, xs, cdf.(r.modelo, xs); color = CORES_DIST[r.distribuicao])
    end
    Legend(fig[2, 1:2], ax; orientation = :horizontal, nbanks = 2)
    salva(fig, "fig07_distribuicao_precipitacao")

    fig = Figure(size = (504, 430))
    for (k, r) in enumerate(eachrow(dist))
        i, j = fldmod1(k, 2); n = length(x)
        qt = quantile.(r.modelo, ((1:n) .- 0.5) ./ n)
        ax = Axis(fig[i, j]; title = r.distribuicao, xlabel = "Quantil teórico (mm)", ylabel = "Quantil observado (mm)",
                  subtitle = "A² = " * fmtnum(r.AD, 2) * ", p(bootstrap) = " * fmtnum(r.p_AD_boot, 2), aspect = 1)
        lim = (0, 1.1maximum([x; qt]))
        lines!(ax, collect(lim), collect(lim); color = COR_TEXTO2, linewidth = 0.8, linestyle = :dash)
        scatter!(ax, qt, sort(x); color = CORES_DIST[r.distribuicao], strokecolor = :white, strokewidth = 0.8)
        limits!(ax, lim, lim)
    end
    salva(fig, "fig08_qqplots_precipitacao")
end

# ----------------------------------------------------------------------------- 10 boxplots climatologia
function fig_boxplots(R)
    dm = R[:dm]
    specs = [(:prec, "Precipitação mensal", "mm", COR_AZUL), (:tmax, "Temperatura máxima (média mensal)", "°C", COR_LARANJA),
             (:tmin, "Temperatura mínima (média mensal)", "°C", COR_AZUL), (:ur, "Umidade relativa (média mensal)", "%", COR_AQUA)]
    fig = Figure(size = (504, 400))
    for (k, (v, tit, u, c)) in enumerate(specs)
        i, j = fldmod1(k, 2)
        d = dropmissing(dm[:, [:mes, v]])
        ax = Axis(fig[i, j]; title = tit, ylabel = u, xtickformat = Makie.automatic, xticks = (1:12, ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]))
        boxplot!(ax, d.mes, Float64.(d[!, v]); color = (c, 0.35), strokecolor = c, strokewidth = 0.8, whiskercolor = c,
                 mediancolor = c, medianlinewidth = 1.5, width = 0.65, markersize = 3.5, outliercolor = c)
    end
    Label(fig[0, :], "Ciclo anual (meses válidos, 2003–2024): estação chuvosa nov–abr concentra ≈ 85 % da chuva";
          fontsize = 9.5, halign = :left, color = COR_TEXTO2)
    salva(fig, "fig09_boxplots_climatologia_mensal")
end

# ----------------------------------------------------------------------------- 11 correlação
function fig_correlacao(R)
    Rm, N, P, vars = R[:cor_R], R[:cor_N], R[:cor_P], R[:cor_vars]
    k = length(vars)
    iu = [(i, j) for i in 1:k for j in 1:k if i < j]
    q = bh_ajuste([P[i, j] for (i, j) in iu]); Q = fill(NaN, k, k)
    for ((i, j), qq) in zip(iu, q); Q[i, j] = qq; Q[j, i] = qq; end
    nomes = Dict(:prec => "Chuva", :r1mm => "Dias ≥1 mm", :rx1day => "Rx1dia", :tmed => "Tméd", :tmax => "Tmáx",
                 :tmin => "Tmín", :dtr => "Amplit.", :ur => "UR", :td => "Td", :vvel => "Vento")
    M = [i > j ? Rm[i, j] : NaN for i in 1:k, j in 1:k]      # triângulo inferior
    fig = Figure(size = (420, 380))
    ax = Axis(fig[1, 1]; xtickformat = Makie.automatic, xticks = (1:k-1, [nomes[v] for v in vars[1:k-1]]),
              ytickformat = Makie.automatic, yticks = (2:k, [nomes[v] for v in vars[2:k]]),
              yreversed = true, xticklabelrotation = π / 4, aspect = 1, title = "Correlação de Spearman — anos hidrológicos",
              subtitle = "valores em negrito: q(BH) < 0,05; n = 14–16", xgridvisible = false, ygridvisible = false,
              leftspinevisible = false, bottomspinevisible = false)
    hm = heatmap!(ax, 1:k-1, 2:k, permutedims(M[2:k, 1:k-1]); colormap = DIV_AZ_VM, colorrange = (-1, 1), nan_color = :transparent)
    for i in 1:k, j in 1:k
        i > j || continue
        r = Rm[i, j]; isnan(r) && continue
        text!(ax, j, i; text = fmtnum(r, 2), align = (:center, :center), fontsize = 7,
              font = Q[i, j] < 0.05 ? :bold : :regular, color = abs(r) > 0.65 ? :white : COR_TEXTO)
    end
    Colorbar(fig[1, 2], hm; label = "ρ de Spearman", width = 8)
    salva(fig, "fig10_matriz_correlacao")
end

# ----------------------------------------------------------------------------- 12 dispersões
function painel_regressao!(ax, d, x, y; cor = COR_AZUL, markersize = 6)
    m = lm(term(y) ~ term(x), d)
    xg = collect(range(minimum(d[!, x]), maximum(d[!, x]), length = 100))
    nd = DataFrame(x => xg)
    p = predict(m, nd; interval = :confidence, level = 0.95)
    band!(ax, xg, Float64.(p.lower), Float64.(p.upper); color = (cor, 0.18))
    scatter!(ax, d[!, x], d[!, y]; color = (cor, markersize < 5 ? 0.5 : 1.0), strokecolor = :white,
             strokewidth = markersize < 5 ? 0 : 0.8, markersize = markersize)
    lines!(ax, xg, Float64.(p.prediction); color = COR_TEXTO, linewidth = 1.3)
    return m
end

function fig_dispersao(R)
    regs = R[:regs]; cm = R[:cor_mensal]
    fig = Figure(size = (504, 430))
    dh = dropmissing(R[:dh][:, [:prec, :tmax, :ur]])
    ax = Axis(fig[1, 1]; xlabel = "Precipitação anual (mm)", ylabel = "Tmáx média (°C)", title = "Chuva × Tmáx (anual)")
    m = painel_regressao!(ax, dh, :prec, :tmax; cor = COR_LARANJA)
    ax.subtitle = "r = " * fmtnum(cor(dh.prec, dh.tmax), 2) * ", R² = " * fmtnum(r2(m), 2) * ", n = $(nrow(dh))"
    ax = Axis(fig[1, 2]; xlabel = "Precipitação anual (mm)", ylabel = "UR média (%)", title = "Chuva × umidade relativa (anual)")
    m = painel_regressao!(ax, dh, :prec, :ur; cor = COR_AQUA)
    ax.subtitle = "r = " * fmtnum(cor(dh.prec, dh.ur), 2) * ", R² = " * fmtnum(r2(m), 2) * ", n = $(nrow(dh))"
    for (k, (x, y, xl, yl, tit, c)) in enumerate([(:ur_anom, :tmax_anom, "Anomalia de UR (p.p.)", "Anomalia de Tmáx (°C)", "UR × Tmáx (anomalias mensais)", COR_LARANJA),
                                                (:tmax_anom, :tmin_anom, "Anomalia de Tmáx (°C)", "Anomalia de Tmín (°C)", "Tmáx × Tmín (anomalias mensais)", COR_AZUL)])
        d = dropmissing(R[:dm][:, [x, y]]); d[!, x] = Float64.(d[!, x]); d[!, y] = Float64.(d[!, y])
        ax = Axis(fig[2, k]; xlabel = xl, ylabel = yl, title = tit)
        painel_regressao!(ax, d, x, y; cor = c, markersize = 4)
        l = cm[((cm.x .== String(x)) .& (cm.y .== String(y))) .| ((cm.x .== String(y)) .& (cm.y .== String(x))), :][1, :]
        ax.subtitle = "r = " * fmtnum(l.r, 2) * " [" * fmtnum(l.ic_inf, 2) * "; " * fmtnum(l.ic_sup, 2) * "], n_ef ≈ " * string(round(Int, l.n_efetivo))
    end
    Legend(fig[3, 1:2], [LineElement(color = COR_TEXTO), PolyElement(color = (COR_CINZA, 0.3))],
           ["reta MQO", "IC 95 % da média"]; orientation = :horizontal)
    salva(fig, "fig11_dispersao_regressoes")
end

# ----------------------------------------------------------------------------- 13 anos extremos
function fig_anos_extremos(R)
    dh = R[:dh]; dh = dh[(dh.ano_hidro .>= 2003) .& (dh.ano_hidro .<= 2023), :]; a = R[:anual]; a = a[a.ano .<= 2024, :]
    fig = Figure(size = (504, 400))
    ax = Axis(fig[1, 1]; title = "Índice de precipitação padronizado (SPI, ano hidrológico)", ylabel = "SPI",
              xtickformat = Makie.automatic, xticks = (dh.ano_hidro, [@sprintf("%02d/%02d", x % 100, (x + 1) % 100) for x in dh.ano_hidro]), xticklabelrotation = π / 2)
    ok = .!ismissing.(dh.spi)
    hspan!(ax, [-1], [1]; color = (COR_NEUTRO, 0.9))
    hlines!(ax, [-2, -1.5, 1.5, 2]; color = COR_TEXTO2, linestyle = :dot, linewidth = 0.7)
    barplot!(ax, dh.ano_hidro[ok], Float64.(dh.spi[ok]); color = [z < 0 ? COR_VERM : COR_AZUL for z in dh.spi[ok]], gap = 0.3)
    scatter!(ax, dh.ano_hidro[.!ok], zeros(count(.!ok)); marker = :xcross, color = COR_CINZA, markersize = 7)
    for (x, z) in zip(dh.ano_hidro[ok], dh.spi[ok])
        abs(z) >= 1.5 && text!(ax, x, z; text = fmtnum(z, 2), fontsize = 7.5, align = (:center, z > 0 ? :bottom : :top), offset = (0, z > 0 ? 2 : -2))
    end
    ylims!(ax, -2.6, 2.6)
    ax2 = Axis(fig[2, 1]; title = "Anomalia padronizada da temperatura média anual", ylabel = "z", xlabel = "Ano civil",
               xticks = 2003:2:2024)
    okt = .!ismissing.(a.z_tmed)
    hspan!(ax2, [-1], [1]; color = (COR_NEUTRO, 0.9))
    hlines!(ax2, [-1.5, 1.5]; color = COR_TEXTO2, linestyle = :dot, linewidth = 0.7)
    barplot!(ax2, a.ano[okt], Float64.(a.z_tmed[okt]); color = [z < 0 ? COR_AZUL : COR_VERM for z in a.z_tmed[okt]], gap = 0.3)
    scatter!(ax2, a.ano[.!okt], zeros(count(.!okt)); marker = :xcross, color = COR_CINZA, markersize = 7)
    for (x, z) in zip(a.ano[okt], a.z_tmed[okt])
        abs(z) >= 1.5 && text!(ax2, x, z; text = fmtnum(z, 2), fontsize = 7.5, align = (:center, z > 0 ? :bottom : :top), offset = (0, z > 0 ? 2 : -2))
    end
    ylims!(ax2, -2.3, 2.3)
    Legend(fig[3, 1], [PolyElement(color = COR_NEUTRO), MarkerElement(marker = :xcross, color = COR_CINZA)],
           ["faixa |z| < 1", "ano excluído (dados insuficientes)"]; orientation = :horizontal)
    salva(fig, "fig12_anos_extremos")
end

# ----------------------------------------------------------------------------- 14 períodos
function fig_periodos(R)
    per = R[:per]; ch, cc = R[:cortes]; dh, a = R[:dh], R[:anual]
    specs = [("prec", dh.ano_hidro, dh.prec, ch, "Precipitação (mm)"), ("tmed", a.ano, a.tmed, cc, "Tméd (°C)"),
             ("tmax", a.ano, a.tmax, cc, "Tmáx (°C)"), ("tmin", a.ano, a.tmin, cc, "Tmín (°C)")]
    fig = Figure(size = (504, 270))
    rng = MersenneTwister(1)
    for (k, (v, t, y, corte, yl)) in enumerate(specs)
        l = per[per.variavel .== v, :][1, :]
        ax = Axis(fig[1, k]; ylabel = yl, xtickformat = Makie.automatic, xticks = (1:2, ["antigo", "recente"]), title = Dict("prec" => "Chuva", "tmed" => "Tméd", "tmax" => "Tmáx", "tmin" => "Tmín")[v],
                  subtitle = "Δ = " * fmtnum(l.dif_medias, v == "prec" ? 0 : 2) * "\ng = " * fmtnum(l.hedges_g, 2) * "\n" * fmtp(l.p_welch))
        for (g, sel, c) in [(1, t .<= corte, COR_AZUL), (2, t .> corte, COR_LARANJA)]
            yy = collect(Float64, skipmissing(y[sel])); n = length(yy)
            scatter!(ax, g .+ 0.12 .* (rand(rng, n) .- 0.5), yy; color = (c, 0.6), markersize = 5)
            e = quantile(TDist(n - 1), 0.975) * std(yy) / sqrt(n)
            rangebars!(ax, [g + 0.28], [mean(yy) - e], [mean(yy) + e]; color = COR_TEXTO, whiskerwidth = 5, linewidth = 1.2)
            scatter!(ax, [g + 0.28], [mean(yy)]; color = COR_TEXTO, markersize = 6)
        end
        xlims!(ax, 0.5, 2.6)
    end
    Label(fig[2, :], "Antigo: 2003/04–2012/13 (chuva) e 2003–2013 (temperatura); recente: demais anos. Barras: média ± IC 95 % (t).";
          fontsize = 8, color = COR_TEXTO2, halign = :left)
    salva(fig, "fig13_comparacao_periodos")
end

# ----------------------------------------------------------------------------- 15 POT/GPD
function fig_pot(R)
    P = R[:pot]; dg = P.diag
    fig = Figure(size = (504, 440))
    ax = Axis(fig[1, 1]; title = "Vida residual média", xlabel = "Limiar u (mm)", ylabel = "E[X − u | X > u] (mm)")
    band!(ax, dg.u, dg.mrl .- 1.96dg.mrl_ep, dg.mrl .+ 1.96dg.mrl_ep; color = (COR_AZUL, 0.18))
    lines!(ax, dg.u, dg.mrl; color = COR_AZUL); vlines!(ax, [P.u]; color = COR_LARANJA, linestyle = :dash, linewidth = 1)
    ax = Axis(fig[1, 2]; title = "Estabilidade do parâmetro de forma", xlabel = "Limiar u (mm)", ylabel = "ξ̂")
    band!(ax, dg.u, dg.ξ .- 1.96dg.ξ_ep, dg.ξ .+ 1.96dg.ξ_ep; color = (COR_AZUL, 0.18))
    lines!(ax, dg.u, dg.ξ; color = COR_AZUL); vlines!(ax, [P.u]; color = COR_LARANJA, linestyle = :dash, linewidth = 1)
    hlines!(ax, [0]; color = COR_TEXTO2, linewidth = 0.6)
    y = sort(P.y); n = length(y); g = GeneralizedPareto(0.0, P.gpd.σ, P.gpd.ξ)
    qt = quantile.(g, ((1:n) .- 0.5) ./ n)
    ax = Axis(fig[2, 1]; title = "Q-Q: excessos × GPD", xlabel = "Quantil GPD (mm)", ylabel = "Excesso observado (mm)",
              subtitle = "σ̂ = " * fmtnum(P.gpd.σ, 1) * ", ξ̂ = " * fmtnum(P.gpd.ξ, 2) * " (EP " * fmtnum(P.gpd.ep_ξ, 2) * ")")
    lim = 1.05maximum([y; qt])
    lines!(ax, [0, lim], [0, lim]; color = COR_TEXTO2, linestyle = :dash, linewidth = 0.8)
    scatter!(ax, qt, y; color = COR_AZUL, strokecolor = :white, strokewidth = 0.6, markersize = 5)
    rl = P.rl
    Ts = exp.(range(log(1.1), log(100), length = 80))
    ax = Axis(fig[2, 2]; xscale = log10, title = "Nível de retorno (chuva diária)", xlabel = "Período de retorno (anos)",
              ylabel = "mm dia⁻¹", xtickformat = Makie.automatic, xticks = ([1, 2, 5, 10, 20, 50, 100], ["1", "2", "5", "10", "20", "50", "100"]),
              subtitle = "u = " * fmtnum(P.u, 1) * " mm; λ = " * fmtnum(P.λ, 1) * "/ano; " * fmtnum(P.anos, 1) * " anos efetivos")
    band!(ax, rl.T_anos, rl.ic_inf, rl.ic_sup; color = (COR_AZUL, 0.18))
    lines!(ax, Ts, nivel_retorno.(P.u, P.gpd.σ, P.gpd.ξ, P.λ, Ts); color = COR_AZUL)
    pk = sort(P.picos; rev = true); m = length(pk)
    scatter!(ax, (m + 1) ./ ((1:m) .* P.λ), pk; color = COR_TEXTO, markersize = 4.5)
    vspan!(ax, [2P.anos], [100]; color = (COR_CINZA, 0.15))
    text!(ax, 95, minimum(pk); text = "extrapolação\n(> 2× registro)", fontsize = 7, color = COR_TEXTO2, align = (:right, :bottom))
    xlims!(ax, 1, 100)
    salva(fig, "fig14_pot_gpd_retorno")
end

# ----------------------------------------------------------------------------- 16 resíduos
function fig_residuos(R, v, rotulo, nome)
    mt = R[:mods][v]; e = residuals(mt.m); f = fitted(mt.m); t = mt.dados.ano; n = length(e)
    fig = Figure(size = (504, 400))
    ax = Axis(fig[1, 1]; title = "Resíduos × tempo", xlabel = "Ano", ylabel = "Resíduo")
    hlines!(ax, [0]; color = COR_TEXTO2, linewidth = 0.6); lines!(ax, t, e; color = (COR_AZUL, 0.4), linewidth = 0.8)
    scatter!(ax, t, e; color = COR_AZUL, markersize = 6)
    ax = Axis(fig[1, 2]; title = "Resíduos × ajustados", xlabel = "Valor ajustado", ylabel = "Resíduo",
              subtitle = "Breusch-Pagan " * fmtp(mt.resumo.p_bp))
    hlines!(ax, [0]; color = COR_TEXTO2, linewidth = 0.6); scatter!(ax, f, e; color = COR_AZUL, markersize = 6)
    ax = Axis(fig[2, 1]; title = "Q-Q normal dos resíduos", xlabel = "Quantil N(0,1)", ylabel = "Resíduo padronizado",
              subtitle = "Shapiro-Wilk " * fmtp(mt.resumo.p_shapiro))
    z = sort(e ./ std(e)); qn = quantile.(Normal(), ((1:n) .- 0.5) ./ n)
    lines!(ax, [-2.5, 2.5], [-2.5, 2.5]; color = COR_TEXTO2, linestyle = :dash, linewidth = 0.8)
    scatter!(ax, qn, z; color = COR_AZUL, markersize = 6)
    ac = acf_irregular(t, e, 5)
    ax = Axis(fig[2, 2]; title = "Autocorrelação dos resíduos", xlabel = "Defasagem (anos)", ylabel = "r_k",
              subtitle = "Ljung-Box " * fmtp(mt.resumo.p_ljungbox) * "; DW = " * fmtnum(mt.resumo.DW, 2))
    hspan!(ax, [-1.96 / sqrt(n)], [1.96 / sqrt(n)]; color = (COR_NEUTRO, 0.9))
    stem!(ax, 1:5, ac; color = COR_AZUL, stemcolor = COR_AZUL, trunkcolor = COR_TEXTO2, markersize = 6)
    ylims!(ax, -1, 1)
    Label(fig[0, :], "Diagnóstico do modelo MQO de tendência — " * rotulo; fontsize = 9.5, halign = :left, color = COR_TEXTO2)
    salva(fig, nome)
end

# ----------------------------------------------------------------------------- 17 contagens
function fig_contagens(R)
    cont = R[:cont]; dh = R[:dh]; a = R[:anual]
    fig = Figure(size = (504, 250))
    d1 = dropmissing(a[:, [:ano, :tx90p_n, :tx_dias]])
    l = cont[startswith.(cont.indice, "TX90p"), :][1, :]
    ax = Axis(fig[1, 1]; title = "Dias quentes (TX90p)", ylabel = "% dos dias", xlabel = "Ano", xticks = 2004:5:2024,
              subtitle = "razão de chances/déc. " * fmtnum(l.razao_decada, 2) * " [" * fmtnum(l.ic_inf, 2) * "; " * fmtnum(l.ic_sup, 2) * "]")
    scatter!(ax, d1.ano, 100 .* d1.tx90p_n ./ d1.tx_dias; color = COR_LARANJA, strokecolor = :white, strokewidth = 0.8)
    tg = collect(range(minimum(d1.ano), maximum(d1.ano), length = 50))
    lines!(ax, tg, 100 .* predict(l.modelo, DataFrame(t = tg .- mean(d1.ano))); color = COR_TEXTO, linewidth = 1.2)
    hlines!(ax, [10]; color = COR_TEXTO2, linestyle = :dot, linewidth = 0.7)
    d2 = dropmissing(dh[:, [:ano_hidro, :r95_n, :dias_validos]])
    l = cont[startswith.(cont.indice, "R95"), :][1, :]
    ax = Axis(fig[1, 2]; title = "Dias de chuva intensa (> P95)", ylabel = "dias por ano hidrológico", xlabel = "Ano hidrológico (início)",
              xticks = 2004:5:2024,
              subtitle = "razão de taxas/déc. " * fmtnum(l.razao_decada, 2) * " [" * fmtnum(l.ic_inf, 2) * "; " * fmtnum(l.ic_sup, 2) * "]")
    barplot!(ax, d2.ano_hidro, Float64.(d2.r95_n); color = COR_AZUL, gap = 0.3)
    tg = collect(range(minimum(d2.ano_hidro), maximum(d2.ano_hidro), length = 50))
    lines!(ax, tg, predict(l.modelo, DataFrame(t = tg .- mean(d2.ano_hidro)); offset = fill(log(365.25), 50)); color = COR_TEXTO, linewidth = 1.2)
    salva(fig, "fig15_indices_extremos_contagem")
end

function gera_figuras(R)
    println("\nGerando figuras…")
    with_theme(TEMA_ARTIGO) do
        fig_completude(R); fig_radiacao(R); fig_prec_anual(R); fig_temp_anual(R); fig_tendencias(R)
        fig_anomalias_mensais(R); fig_distribuicoes(R); fig_boxplots(R); fig_correlacao(R); fig_dispersao(R)
        fig_anos_extremos(R); fig_periodos(R); fig_pot(R); fig_contagens(R)
        fig_residuos(R, "tmax", "temperatura máxima média anual", "fig16_residuos_tmax")
        fig_residuos(R, "prec", "precipitação do ano hidrológico", "fig17_residuos_precipitacao")
    end
end
