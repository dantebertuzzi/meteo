# Menos chuva e mais calor em Petrolina-PE (1961–2025)

Artigo estatístico reprodutível (LaTeX + Julia) sobre a recorrência de anos secos e as
mudanças no regime meteorológico de Petrolina, no semiárido de Pernambuco, a partir da
estação convencional do INMET 82983, validada contra a estação automática A307, estações
convencionais vizinhas, pluviômetros independentes da ANA e as normais oficiais do INMET, com
uma análise regional de frequência (índice de cheia com tendência) para o período de retorno.

**Autor:** Dante de Souza Bertuzzi (pesquisador independente) ·
ORCID [0009-0007-5223-2814](https://orcid.org/0009-0007-5223-2814)

## Perguntas

- **P1. Recorrência de seca.** Qual é o período de retorno de um ano hidrológico tão seco
  quanto 2012/13, como ele muda sob não estacionariedade, e existe persistência de anos secos?
- **P2. Mudança no regime.** Houve mudança na temperatura, na umidade, no vento, nos
  extremos (ETCCDI) e na sazonalidade da chuva?

## Reprodução

Requer Julia 1.13 e uma distribuição TeX com `biblatex`/`biber`.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. main.jl                    # A307: limpeza, agregação, preenchimento (~6 min)
julia --project=. scripts/baixa_ana.jl       # opcional: rebaixa os pluviômetros da ANA (já em data/raw/ana/)
python3 scripts/extrai_normais_inmet.py      # opcional: rebaixa as normais do INMET (requer openpyxl e xlrd)
julia --project=. -t 4 scripts/paper.jl      # análises do artigo → paper/resultados.tex + figuras (~3 min)
julia --project=. test/runtests.jl           # testes unitários dos métodos implementados à mão
cd paper && latexmk -pdf main.tex            # artigo
```

`julia --project=. scripts/paper.jl figuras` regenera só as figuras a partir do cache da última
execução.

### Verificação independente

`scripts/verifica.py` reimplementa os cálculos centrais em Python (SciPy, pymannkendall) a partir
dos dados brutos, sem usar o código Julia, e compara 70 números com as macros do artigo
(dados, ajustes de distribuição, período de retorno, tendências, modelo não estacionário, SPI,
Markov, teleconexão, referência regional, anos chuvosos recentes e pluviômetros da ANA). Também confere a consistência interna das macros.

```bash
python3 -m venv .venv-verif
.venv-verif/bin/pip install -r scripts/requirements-verificacao.txt
.venv-verif/bin/python scripts/verifica.py      # código de saída 1 se algo não conferir
```

**Nenhum número do texto é digitado à mão.** Todos vêm de macros em `paper/resultados.tex`,
escritas pelo script.

## Estrutura

```
data/raw/bdmep/        séries diárias do BDMEP (82983 e vizinhas)
data/raw/indices/      ONI, TNA, TSA (NOAA PSL)
data/raw/inmet_normais/ Normais Climatológicas do INMET e a linha da 82983 (normais_82983.csv)
INMET_*_A307_*.CSV     dados horários da A307 (INMET, pacotes anuais)
main.jl, src/0*.jl     pipeline da A307 (auditoria, agregação, preenchimento, diagnóstico de ausentes)
src/paper/             módulos do artigo (SPI/SPEI, retorno, Markov, homogeneidade, P2, figuras)
scripts/paper.jl       roda todas as análises do artigo
scripts/baixa_ana.jl   baixa a chuva diária dos pluviômetros da ANA (Hidroweb)
scripts/extrai_normais_inmet.py  baixa as normais do INMET e extrai a linha da 82983
scripts/verifica.py    verificação independente em Python dos números do artigo
test/runtests.jl       testes contra valores publicados (FAO-56) e propriedades conhecidas
paper/                 main.tex, refs.bib, resultados.tex (gerado), figuras/ (geradas)
results/paper/         tabelas CSV de todas as análises do artigo
```

## Dados

- INMET, Banco de Dados Meteorológicos (BDMEP): estação convencional 82983 (Petrolina) e
  vizinhas 82886, 82979, 82986, 82753, 82987, 82985, 83088, 82887; dados diários.
- INMET, dados históricos anuais das estações automáticas: A307 (Petrolina), 2003–2024.
- NOAA Physical Sciences Laboratory: índices ONI, TNA e TSA.

## Licenças

- Código: MIT (`LICENSE`).
- Texto, figuras e tabelas do artigo: CC-BY 4.0 (`paper/LICENSE`).
- Os dados brutos do INMET e da NOAA seguem os termos das instituições de origem.

## Como citar

Ver `CITATION.cff`. O DOI do Zenodo será adicionado após o depósito.
