# =============================================================================
# 00_config.jl — pacotes, caminhos e parâmetros metodológicos
# Todas as escolhas "arbitrárias" da análise ficam aqui, num único lugar,
# para que possam ser alteradas e a análise de sensibilidade seja trivial.
# =============================================================================
using CSV, DataFrames, Dates, Statistics, StatsBase, Distributions,
      HypothesisTests, GLM, StatsModels, CairoMakie, Printf, Random, LinearAlgebra

const RAIZ      = normpath(joinpath(@__DIR__, ".."))
const DIR_BRUTO = RAIZ                                   # CSVs originais do INMET
const DIR_PROC  = joinpath(RAIZ, "data", "processed")    # séries agregadas
const DIR_TAB   = joinpath(RAIZ, "results", "tabelas")   # tabelas estatísticas
const DIR_FIG   = joinpath(RAIZ, "figures")

const ESTACAO = "A307"          # Petrolina-PE. O arquivo A807 (Curitiba) é outra cidade.
const FUSO    = Hour(-3)        # INMET registra em UTC; dia climatológico em hora local

# --- critérios de completude --------------------------------------------------
const MIN_HORAS_DIA   = 20      # dia válido p/ temperatura, UR etc. (≥ 20 de 24 h)
const MIN_HORAS_PREC  = 23      # dia válido p/ índices diários de chuva
const MAX_DIAS_FALTA  = 10      # mês válido (regra OMM): ≤ 10 dias ausentes…
const MAX_DIAS_CONSEC = 4       # …e ≤ 4 dias ausentes consecutivos
const MAX_HORAS_FALTA_MES_PREC = 0.02  # mês de chuva válido: ≤ 2 % das horas ausentes
const MIN_MESES_ANO   = 10      # ano (civil/hidrológico) válido p/ médias: ≥ 10 meses
const MAX_FCA         = 0.05    # ano hidrológico de chuva válido: fração climatológica ausente ≤ 5 %
const MES_INICIO_HIDRO = 7      # ano hidrológico jul–jun (estação chuvosa nov–abr inteira)

# --- limites físicos (checagem de plausibilidade) ------------------------------
# rad: 4.900 kJ/m² ≈ constante solar (1.361 W/m²) × 3.600 s — teto no topo da atmosfera
const LIMITES = Dict(
    :prec => (0.0, 150.0), :t => (5.0, 48.0), :tmax => (5.0, 48.0), :tmin => (5.0, 48.0),
    :td => (-10.0, 35.0), :tdmax => (-10.0, 35.0), :tdmin => (-10.0, 35.0),
    :ur => (1.0, 100.0), :urmax => (1.0, 100.0), :urmin => (1.0, 100.0),
    :p => (900.0, 1050.0), :pmax => (900.0, 1050.0), :pmin => (900.0, 1050.0),
    :rad => (0.0, 4900.0), :vvel => (0.0, 40.0), :vraj => (0.0, 60.0), :vdir => (0.0, 360.0))

const Q_LIMIAR_POT = 0.90    # limiar POT = quantil dos dias úmidos (≥ 1 mm); ver diagnóstico

const ALFA   = 0.05
const NBOOT  = 4000
const SEMENTE = 20240928

# --- estilo gráfico (paleta validada p/ daltonismo; ver 09_figuras.jl) ----------
const COR_AZUL    = colorant"#2a78d6"
const COR_LARANJA = colorant"#eb6834"
const COR_AQUA    = colorant"#1baf7a"
const COR_VERM    = colorant"#e34948"
const COR_CINZA   = colorant"#8a8985"
const COR_TEXTO   = colorant"#0b0b0b"
const COR_TEXTO2  = colorant"#52514e"
const COR_NEUTRO  = colorant"#f0efec"

for d in (DIR_PROC, DIR_TAB, DIR_FIG)
    mkpath(d)
end
