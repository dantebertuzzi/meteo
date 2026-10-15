# =============================================================================
# utilidades.jl — números do texto como macros LaTeX, tema gráfico e E/S
#
# Regra do artigo: nenhum número é digitado à mão no main.tex. Cada valor citado
# é registrado aqui com `registra!` e escrito em paper/resultados.tex como
# \newcommand{\Nome}{\num{...}} (siunitx cuida da vírgula decimal).
# =============================================================================

const DIR_PAPER  = joinpath(RAIZ, "paper")
const DIR_PFIG   = joinpath(DIR_PAPER, "figuras")
const DIR_PTAB   = joinpath(RAIZ, "results", "paper")
const LARG_TEXTO = 455.0          # \textwidth em pt (A4, margens de 2,5 cm)
const LARG_COL   = LARG_TEXTO     # artigo em coluna única

for d in (DIR_PAPER, DIR_PFIG, DIR_PTAB)
    mkpath(d)
end

# ----------------------------------------------------------------------------- macros
const MACROS = Dict{String, String}()

"Nome de macro LaTeX válido: só letras. Dígitos viram palavras (T20 → Tvinte)."
function nome_macro(s::AbstractString)
    dig = Dict('0' => "zero", '1' => "um", '2' => "dois", '3' => "tres", '4' => "quatro",
               '5' => "cinco", '6' => "seis", '7' => "sete", '8' => "oito", '9' => "nove")
    out = join(isdigit(c) ? dig[c] : string(c) for c in s if isletter(c) || isdigit(c))
    all(c -> 'a' <= c <= 'z' || 'A' <= c <= 'Z', out) || error("macro inválida: $s → $out")
    return out
end

"""
    registra!(nome, valor; dig = 2, sig = nothing)

Registra um número para o texto. `dig` casas decimais fixas, ou `sig` algarismos
significativos (útil para períodos de retorno, que variam em ordens de grandeza).
Inteiros saem sem casas. Texto (String) é gravado literalmente.
"""
function registra!(nome::AbstractString, valor; dig = 2, sig = nothing)
    k = nome_macro(nome)
    v = if valor isa AbstractString
        valor
    elseif valor isa Integer
        "\\num{$(valor)}"
    elseif ismissing(valor) || !isfinite(valor)
        "\\textemdash"
    else
        x = isnothing(sig) ? round(valor; digits = dig) : round(valor; sigdigits = sig)
        s = isnothing(sig) ? @sprintf("%.*f", dig, x) : fmt_sig(x, sig)
        "\\num{$s}"
    end
    MACROS[k] = v
    return v
end

function fmt_sig(x, sig)
    x == 0 && return "0"
    casas = max(0, sig - 1 - floor(Int, log10(abs(x))))
    @sprintf("%.*f", casas, x)
end

"Registra estimativa e IC 95 % de uma vez: \\<nome>, \\<nome>ICinf, \\<nome>ICsup."
function registra_ic!(nome, est, inf, sup; kw...)
    registra!(nome, est; kw...); registra!(nome * "ICinf", inf; kw...); registra!(nome * "ICsup", sup; kw...)
end

function escreve_macros(caminho = joinpath(DIR_PAPER, "resultados.tex"))
    open(caminho, "w") do io
        println(io, "% Gerado automaticamente por scripts/paper.jl — NÃO EDITAR À MÃO.")
        println(io, "% ", Dates.format(now(), "yyyy-mm-dd HH:MM"), " | ", length(MACROS), " macros")
        for k in sort(collect(keys(MACROS)))
            println(io, "\\newcommand{\\", k, "}{", MACROS[k], "}")
        end
    end
    return caminho
end

salva_ptab(df, nome) = CSV.write(joinpath(DIR_PTAB, nome * ".csv"), df)

# ----------------------------------------------------------------------------- figuras
"Tema único do artigo: fonte serifada do LaTeX (New Computer Modern), sem moldura superior/direita."
function tema_artigo()
    merge(theme_latexfonts(), Theme(
        fontsize = 9,
        Axis = (spinewidth = 0.6, xtickwidth = 0.6, ytickwidth = 0.6, xgridvisible = false,
                ygridvisible = true, ygridcolor = (COR_CINZA, 0.25), ygridwidth = 0.4,
                topspinevisible = false, rightspinevisible = false,
                titlefont = :regular, titlesize = 9, titlealign = :left),
        Legend = (framevisible = false, labelsize = 8, patchsize = (12, 8)),
        Lines = (linewidth = 1.0,), Scatter = (markersize = 4,)))
end

"Salva a figura em PDF vetorial no tamanho final (1 unidade = 1 pt)."
function salva_pfig(fig, nome)
    caminho = joinpath(DIR_PFIG, nome * ".pdf")
    save(caminho, fig; pt_per_unit = 1)
    println("  figura: ", nome)
    return caminho
end

figura(altura; largura = LARG_COL) = Figure(size = (largura, altura), figure_padding = (2, 6, 2, 2))

"p-valor para o texto: < 0,001 quando muito pequeno, 3 casas abaixo de 0,1 e 2 casas acima."
function registra_p!(nome, p)
    k = nome_macro(nome)
    # \PXVal: só o valor (para tabelas); \PX: com o sinal, para o texto escrever $p$\PX ("p = 0,03" ou "p < 0,001")
    MACROS[k * "Val"] = p < 0.001 ? "\\num{< 0.001}" : "\\num{" * @sprintf("%.*f", p < 0.1 ? 3 : 2, p) * "}"
    MACROS[k] = p < 0.001 ? "\\ensuremath{<}\\,\\num{0.001}" : "\\ensuremath{=}\\," * MACROS[k * "Val"]
end
