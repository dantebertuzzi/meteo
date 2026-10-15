# =============================================================================
# 01_importacao.jl — leitura dos CSVs horários do INMET
#
# Os arquivos têm dois layouts de cabeçalho (≤ 2018: "DATA (YYYY-MM-DD)"/"00:00";
# ≥ 2019: "Data"/"0000 UTC"), codificação ISO-8859-1, separador ';', vírgula
# decimal e ausentes codificados como -9999 ou campo vazio. O mapeamento de
# colunas é feito pelo *texto* do cabeçalho, não pela posição.
# =============================================================================

# Ordem importa: o primeiro prefixo que casar vence (ex.: "TEMPERATURA ORVALHO MAX"
# precisa ser testado antes de "TEMPERATURA DO PONTO DE ORVALHO").
const MAPA_COLUNAS = [
    "DATA"                                      => :data,
    "HORA"                                      => :hora,
    "PRECIPITAÇÃO TOTAL"                        => :prec,   # mm (acumulado na hora)
    "PRESSAO ATMOSFERICA AO NIVEL DA ESTACAO"   => :p,      # hPa (mB)
    "PRESSÃO ATMOSFERICA MAX"                   => :pmax,
    "PRESSÃO ATMOSFERICA MIN"                   => :pmin,
    "RADIACAO GLOBAL"                           => :rad,    # kJ/m² (acumulado na hora)
    "TEMPERATURA DO AR - BULBO SECO"            => :t,      # °C (instantânea)
    "TEMPERATURA DO PONTO DE ORVALHO"           => :td,     # °C
    "TEMPERATURA MÁXIMA NA HORA ANT"            => :tmax,   # °C (máx. da hora anterior)
    "TEMPERATURA MÍNIMA NA HORA ANT"            => :tmin,   # °C (mín. da hora anterior)
    "TEMPERATURA ORVALHO MAX"                   => :tdmax,
    "TEMPERATURA ORVALHO MIN"                   => :tdmin,
    "UMIDADE REL. MAX"                          => :urmax,  # %
    "UMIDADE REL. MIN"                          => :urmin,
    "UMIDADE RELATIVA DO AR"                    => :ur,
    "VENTO, DIREÇÃO"                            => :vdir,   # graus
    "VENTO, RAJADA"                             => :vraj,   # m/s
    "VENTO, VELOCIDADE"                         => :vvel,   # m/s
]
const VARS_NUM = [:prec, :p, :pmax, :pmin, :rad, :t, :td, :tmax, :tmin,
                  :tdmax, :tdmin, :urmax, :urmin, :ur, :vdir, :vraj, :vvel]

function mapeia_coluna(cab::AbstractString)
    h = uppercase(strip(cab))
    for (pref, nome) in MAPA_COLUNAS
        startswith(h, pref) && return nome
    end
    return nothing
end

latin1_para_utf8(bytes::Vector{UInt8}) = String(map(Char, bytes))  # Latin-1 ≡ U+0000–U+00FF

parse_num(x) = ismissing(x) ? missing : parse(Float64, replace(strip(x), ',' => '.'))

"""
    le_inmet(caminho) -> DataFrame

Lê um CSV horário do INMET e devolve um DataFrame com nomes padronizados,
`dt_utc`, `dt_local` e metadados do arquivo (layout do cabeçalho).
"""
function le_inmet(caminho::AbstractString)
    txt = latin1_para_utf8(read(caminho))
    df = CSV.read(IOBuffer(txt), DataFrame; delim=';', header=1, types=String,
                  missingstring=["", "-9999"], silencewarnings=true)
    nomes = mapeia_coluna.(names(df))
    manter = findall(!isnothing, nomes)                    # descarta coluna vazia final
    df = df[:, manter]
    rename!(df, Symbol.(nomes[manter]))
    faltando = setdiff([:data, :hora, VARS_NUM...], propertynames(df))
    isempty(faltando) || error("Colunas não encontradas em $(basename(caminho)): $faltando")

    layout = occursin("UTC", df.hora[1]) ? "2019+ (Data/0000 UTC)" : "≤2018 (DATA/00:00)"
    data = Date.(df.data, dateformat"dd/mm/yyyy")          # ambos os layouts usam dd/mm/aaaa
    hora = parse.(Int, first.(replace.(df.hora, ":" => ""), 2))
    out = DataFrame(dt_utc = DateTime.(data) .+ Hour.(hora))
    for v in VARS_NUM
        out[!, v] = parse_num.(df[!, v])
    end
    out.dt_local = out.dt_utc .+ FUSO
    out.arquivo .= basename(caminho)
    out.layout  .= layout
    return out
end

"""
    carrega_estacao(dir, codigo) -> (dados_horarios, inventario)

Lê todos os arquivos da estação `codigo`, concatena e devolve também um
inventário (um registro por arquivo).
"""
function carrega_estacao(dir::AbstractString, codigo::AbstractString)
    todos = filter(f -> endswith(uppercase(f), ".CSV") && startswith(f, "INMET_"), readdir(dir))
    arqs  = filter(f -> occursin("_$(codigo)_", f), todos)
    ignorados = setdiff(todos, arqs)
    isempty(ignorados) || @info "Arquivos de outras estações ignorados" ignorados
    partes = [le_inmet(joinpath(dir, f)) for f in arqs]
    inventario = DataFrame(
        arquivo  = arqs,
        layout   = [p.layout[1] for p in partes],
        linhas   = nrow.(partes),
        inicio   = [minimum(p.dt_utc) for p in partes],
        fim      = [maximum(p.dt_utc) for p in partes])
    d = sort!(vcat(partes...), :dt_utc)
    return d, inventario
end
