# Carrega a configuração, os utilitários já existentes da análise e os módulos do artigo.
const _SRC = normpath(joinpath(@__DIR__, ".."))
include(joinpath(_SRC, "00_config.jl"))
using SpecialFunctions: gamma
using Optim, Turing, MissingPatterns
for f in ["01_importacao", "02_limpeza", "02b_preenchimento", "03_descritiva", "04_temporal",
          "05_distribuicoes", "06_extremos"]
    include(joinpath(_SRC, f * ".jl"))
end
for f in ["utilidades", "indices_seca", "dados_longos", "retorno", "markov", "homogeneidade", "p2", "figuras"]
    include(joinpath(@__DIR__, f * ".jl"))
end
