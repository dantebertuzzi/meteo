# =============================================================================
# scripts/baixa_ana.jl — baixa a chuva diária dos pluviômetros da rede da ANA
# (Hidroweb) em torno de Petrolina e grava um CSV por posto em data/raw/ana/.
#
# Fonte: serviço público da ANA, HidroSerieHistorica (tipoDados = 2, chuva).
# Para cada mês, usa o nível de consistência mais alto disponível
# (2 = consistido; 1 = bruto). A seleção dos postos usados no artigo
# (≥ 30 anos hidrológicos válidos em comum com a 82983) é feita em
# scripts/paper.jl; aqui se baixam todos os candidatos de POSTOS_ANA.
#
# Uso: julia --project=. scripts/baixa_ana.jl
# =============================================================================
using Downloads, Dates, Printf
include(joinpath(@__DIR__, "..", "src", "paper", "postos_ana.jl"))

const URL_ANA = "http://telemetriaws1.ana.gov.br/ServiceANA.asmx/HidroSerieHistorica"
const DIR_ANA = joinpath(@__DIR__, "..", "data", "raw", "ana")

campo(bloco, nome) = (m = match(Regex("<$(nome)>([^<]*)</$(nome)>"), bloco); isnothing(m) ? "" : strip(m.captures[1]))

"Lê o XML da ANA e devolve Dict(Date do 1º dia do mês => (nível, vetor de chuvas diárias ou missing))."
function le_xml_ana(xml)
    meses = Dict{Date, Tuple{Int, Vector{Union{Missing, Float64}}}}()
    for m in eachmatch(r"<SerieHistorica[^>]*>(.*?)</SerieHistorica>"s, xml)
        b = m.captures[1]
        dh = campo(b, "DataHora"); isempty(dh) && continue
        d0 = Date(dh[1:10]); nivel = something(tryparse(Int, campo(b, "NivelConsistencia")), 1)
        vals = Union{Missing, Float64}[]
        for i in 1:daysinmonth(d0)
            s = campo(b, @sprintf("Chuva%02d", i))
            push!(vals, isempty(s) ? missing : parse(Float64, replace(s, "," => ".")))
        end
        (!haskey(meses, d0) || nivel > meses[d0][1]) && (meses[d0] = (nivel, vals))
    end
    return meses
end

function baixa_posto(cod; inicio = "01/01/1955", fim = "31/12/2025")
    url = "$(URL_ANA)?codEstacao=$(lpad(cod, 8, '0'))&dataInicio=$(inicio)&dataFim=$(fim)&tipoDados=2&nivelConsistencia="
    xml = String(take!(Downloads.download(url, IOBuffer(); timeout = 300)))
    meses = le_xml_ana(xml)
    arq = joinpath(DIR_ANA, "ana_$(cod).csv")
    open(arq, "w") do io
        println(io, "data,prec,nivel_consistencia")
        for d0 in sort(collect(keys(meses)))
            nivel, vals = meses[d0]
            for (i, v) in enumerate(vals)
                println(io, Date(year(d0), month(d0), i), ",", ismissing(v) ? "" : v, ",", nivel)
            end
        end
    end
    @printf("%s %-26s %4d meses → %s\n", cod, POSTOS_ANA[cod].nome, length(meses), basename(arq))
end

if abspath(PROGRAM_FILE) == @__FILE__
    mkpath(DIR_ANA)
    for cod in sort(collect(keys(POSTOS_ANA)))
        baixa_posto(cod)
    end
end
