using StatsBase, Plots

# Gerar dados AR(1)
n = 1000
x = zeros(n)
x[1] = randn()
for i in 2:n
    x[i] = 0.8 * x[i-1] + randn()
end

max_lag = 20

# Calcular ACF e PACF
acf_values = autocor(x, 0:max_lag)
pacf_values = pacf(x, 1:max_lag)

# Calcular intervalo de confiança 95%
conf = 1.96 / sqrt(n)

# Função para adicionar o sombreamento (sem borda)
function add_confidence_band!(plot, lags, conf)
    x_min = minimum(lags)
    x_max = maximum(lags)
    plot!(
        plot,
        [x_min, x_max, x_max, x_min, x_min],
        [conf, conf, -conf, -conf, conf],
        seriestype=:shape,
        fillcolor=:lightblue,
        fillalpha=0.3,
        linewidth=0,  # Removemos a borda aqui
        label="95% CI"
    )
end

# Plot ACF
p1 = plot(0:max_lag, acf_values,
    seriestype=:stem,
    title="ACF",
    xlabel="Lag",
    ylabel="Correlation",
    markershape=:circle,
    ylim=(-0.5, 1.1),
    legend=:topright
)
add_confidence_band!(p1, 0:max_lag, conf)

# Plot PACF
p2 = plot(1:max_lag, pacf_values,
    seriestype=:stem,
    title="PACF",
    xlabel="Lag",
    ylabel="Partial Correlation",
    markershape=:circle,
    ylim=(-0.5, 1.1),
    legend=:topright
)
add_confidence_band!(p2, 1:max_lag, conf)

# Combinar plots
plot(p1, p2, layout=(2, 1), size=(800, 600))