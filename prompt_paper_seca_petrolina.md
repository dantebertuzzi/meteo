# Tarefa: artigo estatístico reprodutível (LaTeX + Julia) sobre seca e mudança climática em Petrolina-PE

## Contexto e papel
Você é um estatístico com experiência em climatologia aplicada, hidrologia estatística e teoria de valores extremos. Vai me ajudar a escrever meu primeiro artigo, que será depositado no Zenodo com DOI (licença CC-BY 4.0). O artigo precisa ter rigor de periódico de estatística aplicada: métodos justificados, incerteza sempre quantificada, limitações explícitas e referências reais.

Sou estudante de bacharelado em Estatística e desenvolvedor Julia. Sou autor do MissingPatterns.jl, que já usei no diagnóstico de ausentes. Todo o código deve ser em Julia, e as figuras em CairoMakie.

## Perguntas de pesquisa
**P1 — Recorrência de seca.** Com base nos dados, qual é a probabilidade de ocorrer outra seca comparável à de 2012/13 em Petrolina? Responda em termos de:
  - período de retorno (com IC) da chuva do ano hidrológico 2012/13 e do evento de seca medido por SPI/SPEI (duração, severidade, intensidade);
  - probabilidade de ≥ 1 evento igual ou pior nos próximos 5, 10, 20 e 30 anos, assumindo estacionariedade, e como isso muda sob o modelo não estacionário, se houver evidência de tendência;
  - persistência: dado um ano seco, qual a probabilidade de o próximo também ser seco? (A seca do Nordeste de 2012 a ~2017 foi plurianual.)
  - **Não** apresente previsão determinística de "ano da próxima seca". Se fizer análise espectral ou de ondaletas, diga explicitamente que 20–22 anos não permitem identificar ciclos com período maior que ~10 anos.

**P2 — Mudança e perfil climático.** O regime meteorológico de Petrolina mudou de 2003 a 2025? Caracterize:
  - sazonalidade: climatologia mensal, índice de sazonalidade (Walsh & Lawler), centroide e concentração da estação chuvosa (estatística circular), início e fim da estação chuvosa (método de Liebmann et al.) e se essas datas se deslocaram;
  - tendências em temperatura (Tmáx, Tmín, Tméd, amplitude térmica), UR, vento e índices de extremos ETCCDI (TXx, TNn, TX90p, TN10p, R1mm, Rx1day, SDII, R95p, CDD = dias secos consecutivos, CWD);
  - mudanças na distribuição, e não só na média: regressão quantílica da Tmáx/Tmín diária no tempo (quantis 0,1, 0,5 e 0,9);
  - comparação com as Normais Climatológicas do INMET 1991–2020 da estação 82983, quando disponível, como linha de base externa.

Primeiro, verifique nos dados se 2012/13 é de fato o ano hidrológico mais seco da série, ou se outro ano (ex.: 2011/12, 2015/16) é mais extremo. Defina o "evento de referência" com precisão, tanto por total anual quanto por SPI-12/SPEI-12, antes de qualquer cálculo.

## Dados disponíveis

### Fontes
1. **A307 – estação automática do INMET, Petrolina-PE** (série principal)
   - Dados horários, 21/02/2003 → 31/12/2024. São 22 CSVs anuais na raiz do projeto (`INMET_NE_PE_A307_PETROLINA_*.CSV`), com 191 640 registros.
   - Formato: ISO-8859-1, separador `;`, vírgula decimal, ausentes = `-9999` ou campo vazio. Há dois layouts de cabeçalho (≤ 2018: `DATA (YYYY-MM-DD)`/`00:00`; ≥ 2019: `Data`/`0000 UTC`), mas as datas estão sempre em `dd/mm/aaaa`. Os horários são UTC e foram convertidos para hora local (UTC−3) antes da agregação diária.
   - Ausência: cerca de 10,3 % das horas. As falhas são **em bloco**: quando o registrador para, todas as variáveis somem juntas (padrão confirmado com MissingPatterns.jl). Maiores lacunas: abr–jun/2010, abr/2013, mar–jun/2018, ago/2022–abr/2023 e a maior parte de 2024.
   - O CSV `INMET_S_PR_A807_CURITIBA_*` está na pasta, mas é de outra cidade e é **ignorado**.

2. **82983 – estação convencional do INMET, Petrolina-PE** (usada para preencher lacunas)
   - Lat −9,389, lon −40,523, altitude 372,5 m. Obtida no BDMEP, 01/01/2003 → 30/06/2025.
   - Diário: `data/raw/bdmep/dados_82983_D_2003-01-01_2025-06-30.csv`
   - "Horário": `data/raw/bdmep/dados_82983_H_2003-01-01_2025-06-30.csv`. São só 3 leituras por dia (00, 12 e 18 UTC). **Ainda não é usado no pipeline.**
   - Formato: UTF-8, 10 linhas de metadados, depois cabeçalho `Data Medicao;...`; separador `;`, vírgula decimal, ausentes = `null`, datas `AAAA-MM-DD`.
   - A chuva diária é lida às 12 UTC (9 h local), então o "dia" dela não coincide exatamente com o dia civil.

### Variáveis originais

**A307 (horária)**

| Código | Variável | Unidade |
|---|---|---|
| prec | precipitação acumulada na hora | mm |
| t | temperatura do ar (bulbo seco), instantânea | °C |
| tmax / tmin | temperatura máx./mín. na hora anterior | °C |
| td / tdmax / tdmin | ponto de orvalho (instantâneo, máx./mín. da hora) | °C |
| ur / urmax / urmin | umidade relativa (instantânea, máx./mín. da hora) | % |
| p / pmax / pmin | pressão ao nível da estação | hPa |
| vvel / vraj / vdir | vento: velocidade, rajada máxima, direção | m/s, m/s, graus |
| rad | radiação global | kJ/m² — **excluída da análise** (mudança de escala 2005–08, quebra ≈ 2014, sem dados 2023–24) |

**82983 (diária)**

| Variável | Unidade |
|---|---|
| precipitação total | mm |
| temperatura máxima, mínima e média compensada | °C |
| umidade relativa média e mínima | % |
| velocidade média do vento | m/s |
| insolação total | h (não usada ainda) |
| evaporação Piche | mm (não usada ainda) |

**82983 ("horária", 3×/dia, não usada):** chuva, bulbo seco, bulbo úmido, ponto de orvalho, UR, pressão (estação e nível do mar), vento (direção e velocidade), nebulosidade, nuvens (CH/CL/CM) e visibilidade.

### Tratamento
- Limites físicos: um valor fora da faixa vira `missing` (ex.: prec 0–150 mm/h, T 5–48 °C, UR 1–100 %, p 900–1050 hPa).
- Dia válido: ≥ 20 de 24 h para temperatura, UR etc.; ≥ 23 h para chuva.
- Mês válido (regra da OMM): ≤ 10 dias ausentes e ≤ 4 dias consecutivos. Para a chuva, ≤ 2 % das horas ausentes.
- Ano civil válido (temperatura, UR): ≥ 10 meses válidos. O valor anual é a climatologia mais a média das anomalias mensais.
- Ano hidrológico: julho a junho. É válido se a fração climatológica ausente (FCA) for ≤ 5 %, isto é, a parte da chuva esperada no ano que cairia em dias sem registro.
- **Preenchimento:** quando falta um valor diário da A307, entra o valor da 82983 corrigido pela diferença média entre as estações em cada mês. A correção soma uma constante em tmed, tmax, tmin, ur e urmin e multiplica por um fator em vvel. A chuva entra sem correção (o total difere < 2 % no período comum). Td e pressão não são preenchidos, porque a 82983 diária não tem essas variáveis.
- Concordância entre as estações nos dias em comum (correlação diária r): Tmáx 0,98; Tméd 0,98; Tmín 0,95; UR 0,83; chuva 0,74; vento 0,68. A UR mínima da 82983 é cerca de 12 pontos percentuais mais alta que a da A307 antes da correção.
- Resultado: 20 anos hidrológicos válidos para chuva (2003/04–2021/22 e 2023/24; ficam de fora 2002/03, 2022/23 e 2024/25) e 19–21 anos civis válidos para temperatura.

### Arquivos processados (`data/processed/`) — os mais importantes
- **`diario.csv`** (8 167 dias, 20/02/2003 → 30/06/2025): `data, nh` (horas registradas), `nh_prec, prec_obs` (soma bruta), `prec` (válida), `tmed, tmax, tmin, dtr, ur, urmin, td, p, vvel, ano, mes, ano_hidro`. Também tem colunas Bool que marcam o dia preenchido com a 82983: `prec_pre, tmax_pre, tmed_pre, tmin_pre, ur_pre, urmin_pre, vvel_pre`.
- **`mensal.csv`** (269 meses): totais e médias mensais, `frac_ausente_prec`, `tempo` (ano decimal) e anomalias `<var>_anom` em relação à climatologia do mês.
- **`anual_civil.csv`** (2003–2025): médias anuais com o número de meses válidos (`nmeses_<var>`) e índices de extremos de temperatura: `txx` (maior Tmáx), `tnn` (menor Tmín), `tx90p_n/tx_dias` (dias com Tmáx acima do P90 do mês), `tn10p_n/tn_dias` (dias com Tmín abaixo do P10).
- **`anual_hidrologico.csv`** (2002/03–2024/25): `prec_obs, fca, valido, prec, dias_validos` e índices de chuva: `r1mm` (dias ≥ 1 mm), `rx1day` (maior chuva diária), `sdii` (intensidade média nos dias úmidos), `r95_n` (dias acima do P95 dos dias úmidos), além das médias das demais variáveis.

### Análises já feitas (`results/tabelas/`, `figures/`)
- Auditoria de qualidade e diagnóstico de ausentes (MissingPatterns.jl).
- Estatística descritiva.
- Tendências: Mann-Kendall, TFPW, Theil-Sen, Pettitt, MQO, MQG-AR(1) e GLM Gamma, com correção de Benjamini-Hochberg.
- Ajuste de distribuições à chuva anual (Gamma, Weibull, LogNormal, Normal).
- SPI anual.
- Extremos diários de chuva por POT/GPD, com níveis de retorno e IC bootstrap.
- Correlações e regressões entre variáveis.
- Comparação de períodos (2003/04–2012/13 contra 2013/14–2023/24).
- Pipeline em Julia: `julia --project=. main.jl`, cerca de 6 min.

### Dados opcionais
Use se eu fornecer; senão, trate como limitação e sugira na discussão:
- série longa da 82983 desde ~1961 (BDMEP);
- índices ONI (ENSO) e TNA/TSA (Atlântico) mensais da NOAA PSL;
- dados em grade (Xavier et al., BR-DWGD; ou CHIRPS) para representar "a região" além do ponto da estação.

## Métodos esperados (justifique cada escolha no texto)
1. **Qualidade e homogeneidade**
   - Resuma a auditoria já feita. Inclua a figura de padrões de ausência gerada com MissingPatterns.jl.
   - Teste a homogeneidade da série combinada A307+82983 com a série de diferenças entre estações e um teste de quebra (Pettitt e/ou SNHT de Alexandersson).
   - Faça análise de sensibilidade: rode os resultados principais com e sem os dias preenchidos (`*_pre`).
2. **Índices de seca**
   - SPI em 3, 6 e 12 meses (McKee et al. 1993): Gamma com massa em zero, ajustada por mês calendário.
   - SPEI (Vicente-Serrano et al. 2010) com ETP de Hargreaves–Samani, a partir de Tmáx, Tmín e radiação extraterrestre calculada pela latitude. Justifique o uso de Hargreaves, já que a radiação medida foi excluída.
   - Compare as distribuições candidatas para SPEI (log-logística etc.; ver Stagge et al. 2015).
   - Teoria de corridas (Yevjevich 1967): eventos com duração, severidade e intensidade.
3. **Período de retorno**
   - Ajuste a chuva anual (já testei Gamma, Weibull, LogNormal e Normal).
   - Selecione o modelo por AIC/BIC, QQ-plot e teste de Anderson-Darling com p-valor por bootstrap paramétrico.
   - Calcule a probabilidade de não excedência do evento de referência e T = 1/p, com IC por bootstrap paramétrico e por verossimilhança perfilada.
   - Calcule P(≥1 evento em n anos) = 1 − (1 − p)^n, propagando a incerteza.
   - **Incluir uma abordagem bayesiana** (Turing.jl ou Extremes.jl), porque a amostra é pequena: priors fracamente informativos e justificados, e distribuição preditiva posterior do período de retorno.
   - Se houver tendência: modelo não estacionário com parâmetro dependente do tempo e/ou de ENSO. Compare com o estacionário por deviance/LOO.
4. **Persistência**
   - Cadeia de Markov de 1ª ordem com estados seco/normal/úmido (tercis ou limiares de SPI-12).
   - Matriz de transição com IC (Dirichlet ou bootstrap) e tempo médio de recorrência ao estado seco.
5. **Tendências e mudanças**
   - Mann-Kendall com correção para autocorrelação (TFPW de Yue et al. 2002 e/ou Hamed & Rao 1998), Theil-Sen e Pettitt.
   - Correção de Benjamini-Hochberg para o conjunto de testes.
   - Reporte a **magnitude com IC**, e não só o p-valor.
   - Inclua uma estimativa de poder ou tendência mínima detectável com n ≈ 20, para interpretar corretamente os resultados não significativos.
   - Decomposição STL (Cleveland et al. 1990) das séries mensais e GAM com suavizador sazonal cíclico + tendência suave como alternativa.
6. **Teleconexões** (se os índices forem fornecidos)
   - Correlação defasada e regressão da chuva da estação chuvosa (nov–abr) com ONI e o gradiente TNA–TSA.
   - Probabilidade de ano seco condicionada à fase do ENSO.

Nenhum número pode aparecer no texto sem vir de saída do código.

## Entregáveis
Estrutura de repositório pronta para Zenodo/GitHub:

```
paper/
  main.tex            % artigo
  refs.bib            % BibTeX, todas as entradas com DOI quando existir
  resultados.tex      % \newcommand gerados pelo Julia com cada número citado no texto
  figuras/*.pdf       % geradas por CairoMakie
src/                  % módulos Julia por etapa (spi.jl, retorno.jl, markov.jl, tendencias.jl, ...)
scripts/paper.jl      % roda tudo e regenera figuras + resultados.tex
Project.toml, Manifest.toml
README.md, LICENSE (MIT para código, CC-BY-4.0 para texto), CITATION.cff, .zenodo.json
```

- Os números do texto devem vir de macros em `resultados.tex` geradas pelo script (ex.: `\TretornoMed`, `\TretornoICinf`), nunca digitados à mão.
- Formato dos números: vírgula decimal e IC 95 % sempre que houver estimativa.

### LaTeX
- `article` (ou `elsarticle` em modo preprint), `biblatex` + `biber` em estilo autor-ano, `siunitx` (vírgula decimal, locale BR), `booktabs`, `hyperref`, `cleveref`, `amsmath`.
- Idioma: português, com *Abstract* em inglês além do *Resumo*.
- Seções:
  1. Introdução (contexto do semiárido, secas no Nordeste, lacuna e perguntas)
  2. Área de estudo
  3. Dados (fontes, tratamento, preenchimento, qualidade)
  4. Métodos (com as equações principais: SPI, SPEI/Hargreaves, Gamma, verossimilhança, período de retorno, cadeia de Markov, MK/TFPW, Sen)
  5. Resultados
  6. Discussão (comparar com a literatura sobre a seca de 2012–2017)
  7. Limitações
  8. Conclusões
  9. Disponibilidade de dados e código
  10. Referências
- Metadados do autor: Dante [SOBRENOME], [INSTITUIÇÃO], ORCID [XXXX].

### Figuras (CairoMakie, estilo limpo e consistente)
- Uma paleta só, fonte serifada compatível com o LaTeX, salvas em PDF vetorial no tamanho final da coluna.
- Mínimo esperado:
  - mapa ou esquema da área;
  - padrões de ausência (MissingPatterns.jl);
  - climatologia mensal com bandas de variabilidade;
  - série de chuva anual com o evento de referência destacado;
  - SPI/SPEI-12 no tempo;
  - ajuste de distribuição com curva de período de retorno e IC;
  - posterior do período de retorno;
  - matriz de transição;
  - painel de tendências (Sen ± IC) por variável;
  - mudança no início/fim da estação chuvosa;
  - regressão quantílica da temperatura.

## Referências — regra rígida
Cite apenas obras que você tem certeza de que existem, com autores, ano, título, veículo e DOI corretos. Se houver qualquer dúvida, marque a entrada com `note = {VERIFICAR}` em vez de inventar. Ponto de partida (confirme cada uma):

- McKee, Doesken & Kleist (1993) — SPI
- Guttman (1999)
- Vicente-Serrano, Beguería & López-Moreno (2010, J. Climate) — SPEI
- Stagge et al. (2015, Int. J. Climatol.)
- Hargreaves & Samani (1985)
- Yevjevich (1967)
- Mann (1945); Kendall (1975)
- Sen (1968, JASA)
- Pettitt (1979, Appl. Stat.)
- Yue et al. (2002, Hydrol. Process.)
- Hamed & Rao (1998, J. Hydrol.)
- Benjamini & Hochberg (1995, JRSS-B)
- Coles (2001)
- Wilks, *Statistical Methods in the Atmospheric Sciences*
- Cleveland et al. (1990) — STL
- Walsh & Lawler (1981)
- Liebmann et al. (2007, J. Climate)
- Zhang et al. (2011, WIREs Clim. Change) — ETCCDI
- Alexandersson (1986)
- Torrence & Compo (1998, BAMS)
- WMO (2017) — Guidelines on the Calculation of Climate Normals (WMO-No. 1203)
- Hastenrath & Heller (1977, QJRMS)
- Moura & Shukla (1981, JAS)
- Marengo, Torres & Alves (2017, Theor. Appl. Climatol.) — seca do Nordeste
- Xavier, King & Scanlon (2016, Int. J. Climatol.)
- Funk et al. (2015, Sci. Data) — CHIRPS
- Bezanson et al. (2017, SIAM Review) — Julia
- Danisch & Krumbiegel (2021, JOSS) — Makie
- referência do Extremes.jl e do Turing.jl, se usados
- MissingPatterns.jl, citado pelo DOI do Zenodo/registro

## Pacotes Julia sugeridos
CSV, DataFrames, Dates, Statistics, StatsBase, Distributions, HypothesisTests, GLM, Optim, Bootstrap, Extremes, Turing, MCMCChains, Loess, CairoMakie, MissingPatterns. Implemente à mão o que não existir de forma confiável (SPI, SPEI, MK modificado, Hargreaves), com testes unitários contra valores conhecidos.

## Como trabalhar
1. **Fase 1 — Plano.** Antes de escrever código, entregue o plano de análise e as decisões em aberto (definição do evento de referência, limiares, distribuições). Espere minha aprovação.
2. **Fase 2 — Código e resultados.** Escreva os módulos e o script. Se você tiver acesso aos arquivos, rode tudo e me mostre as tabelas principais. Se não tiver, escreva o código completo e deixe o LaTeX com as macros de `resultados.tex` como espaços reservados. **Não invente valores.**
3. **Fase 3 — Artigo.** Escreva `main.tex` completo a partir dos resultados reais.
4. **Fase 4 — Revisão.** Revise como um parecerista crítico: aponte afirmações fortes demais, testes sem correção, ICs ausentes e referências não verificadas.

Tom do texto: científico, direto, sem exagerar conclusões. Com n ≈ 20, "não detectamos" ≠ "não existe".
