# Análise estatística — Estação INMET A307, Petrolina-PE (2003–2024)

Reprodução: `julia --project=. main.jl` (≈ 1,5 min). Tabelas completas em `results/tabelas/`, séries agregadas em `data/processed/`, figuras (PDF vetorial + PNG 300 dpi) em `figures/`.

---

## 1. Diagnóstico inicial dos dados

| Item | Resultado |
|---|---|
| Arquivos | 22 CSVs horários da A307 (21/02/2003 → 31/12/2024). O CSV de Curitiba (A807, 2021) foi **excluído**: é outra cidade. |
| Registros | 191.640 horas; 0 timestamps duplicados; grade horária completa (as falhas são valores ausentes, não linhas faltando). |
| Formato | ISO-8859-1, `;`, vírgula decimal. Dois layouts de cabeçalho (≤ 2018 / ≥ 2019); ausentes = `-9999` ou vazio. O cabeçalho antigo diz `YYYY-MM-DD`, mas as datas estão em `dd/mm/aaaa`. |
| Variáveis | Chuva, temperatura (instantânea, máx./mín. horárias), ponto de orvalho, UR, pressão, radiação global, vento (velocidade, rajada, direção). |
| **Inexistentes** | Evapotranspiração, número de dias chuvosos e índices de extremos. Os dois últimos foram **derivados** da série diária. |
| Ausência | ≈ 10,3 % das horas em todas as variáveis (o mesmo registrador falha junto). Blocos longos: abr–jun/2010, abr/2013, mar–jun/2018, ago–dez/2022, jan–mar/2023 e a maior parte de 2024 (fig01). |
| Consistência | 0 violações de: Tmín > Tmáx, T fora de [Tmín; Tmáx], Td > T, URmín > URmáx, Pmín > Pmáx, rajada < média. Nenhuma sequência ≥ 12 h constante em T ou P (1 em UR). |
| Valores extremos | 2 horas com > 50 mm/h (68,2 mm em 14/01/2004; 64,4 mm em 18/03/2008). São plausíveis no semiárido, ficam sinalizadas e **não foram removidas**. |
| **Radiação** | **Inutilizável.** 2005–2007 têm valores ~10× maiores (máx. 50.086 kJ m⁻² h⁻¹, acima do teto físico de ≈ 4.900); 2008 é misto; o p95 cai de ≈ 3.900 para ≈ 2.650 a partir de 2014 (troca/recalibração de sensor); 2023–24 não têm dado (fig02). Excluída de todas as análises. |
| Vento | Queda de nível (≈ 3,9 → 3,4 m/s; 2,6 m/s em 2021), típica de mudança de exposição ou sensor. Reportado com ressalva. |

### Tratamento e agregação (nenhum dado inventado)

1. **Horário → diário**, em hora local (UTC−3).
   - Temperatura e UR: dia válido com ≥ 20 de 24 h.
   - Índices diários de chuva: dia válido com ≥ 23 h.
2. **Diário → mensal**: regra OMM, com ≤ 10 dias ausentes e ≤ 4 consecutivos. Total mensal de chuva só com ≤ 2 % das horas ausentes.
3. **Temperatura e UR anuais (ano civil)**: climatologia + média das **anomalias** mensais, exigindo ≥ 10 meses válidos. Médias brutas seriam enviesadas quando falta uma estação do ano inteira.
4. **Chuva no ano hidrológico (jul–jun)**, validado pela *fração climatológica ausente*:
   $$\mathrm{FCA}=\frac{\sum_{\text{dias}}(h_{\text{aus}}/24)\,\bar c_{\text{mês}}}{\bar P_{\text{anual}}}\le 5\%,$$
   onde $\bar c_{\text{mês}}$ é a chuva diária climatológica do mês. Um julho ausente (≈ 1 % da chuva anual) é tolerado; faltar março–abril (≈ 30 %) não é.

| Resultado | Anos |
|---|---|
| Chuva: **14 anos hidrológicos válidos** | 03/04, 04/05, 06/07, 08/09, 10/11, 11/12, 13/14–16/17, 18/19–21/22 |
| Chuva: excluídos | 05/06 (FCA = 5,2 %, limítrofe), 07/08 (8,8 %), 09/10 (23 %), 12/13 (16 %), 17/18 (25 %), 22/23 e 23/24 |
| Temperatura: **15 anos civis válidos** | 2003–05, 2007, 2009, 2011–17, 2019–21 |

> **Viés de exclusão importante:** 2012/13 registrou só 69 mm *com* abril ausente. Mesmo somando a climatologia de abril (≈ 60 mm), ficaria em torno de 130 mm, abaixo de todos os anos válidos. O ano mais seco da série (a seca de 2012–13 no semiárido) está, portanto, **fora** das estatísticas anuais. Isso tende a subestimar a variabilidade e a cauda seca.

---

## 2. Estratégia estatística (e por que cada método)

| Pergunta | Método | Justificativa |
|---|---|---|
| Tendência (chuva) | Mann-Kendall + Theil-Sen | Chuva é assimétrica; métodos de postos não exigem normalidade e resistem a anos extremos. |
| Tendência (temperatura, UR) | MQO com diagnóstico completo + MK | Resíduos aproximadamente normais; o MQO dá IC e R². O MK serve de checagem robusta. |
| Autocorrelação | r₁, Durbin-Watson, Ljung-Box; MK com TFPW; **MQG-AR(1) em tempo contínuo** | A série tem lacunas; a correlação ρ^\|Δt\| respeita a distância real entre os anos. |
| Mais poder estatístico | MQO nas **anomalias mensais** (n = 225) com EP **Newey-West (HAC)** | A série anual tem só n = 15. As anomalias removem o ciclo sazonal; o HAC corrige a autocorrelação. |
| Ponto de mudança | Pettitt | Localiza quebras e as cruza com trocas instrumentais conhecidas. |
| Distribuição da chuva | Normal, Lognormal, Gamma e Weibull por MV; AIC/BIC, pesos de Akaike, Anderson-Darling com **p por bootstrap paramétrico**, Q-Q | O p tabelado do AD é inválido quando os parâmetros são estimados dos próprios dados. A Exponencial foi descartada (moda em 0 e CV = 1 fixo). |
| Anos secos/chuvosos | SPI (Gamma → normal padrão) | Corrige a assimetria: desvios iguais em mm não são igualmente raros. |
| Extremos diários | POT/GPD com desagrupamento; exposição em "anos climatológicos" | Usa 59 eventos em vez de ~14 máximos anuais (GEV). A exposição corrige as falhas concentradas na estação chuvosa. |
| Contagens de extremos | GLM quase-Poisson / quase-binomial | Dados de contagem com sobredispersão. |
| Associação | Pearson e Spearman; anomalias mensais com n efetivo de Bretherton | Remove a correlação espúria do ciclo sazonal e corrige a dependência temporal. |
| Comparação de períodos | Welch, Mann-Whitney, g de Hedges, δ de Cliff, Hodges-Lehmann e IC bootstrap | Efeito + IC, não só p-valor. |
| Multiplicidade | Benjamini-Hochberg (q) em cada bateria de testes | São dezenas de testes simultâneos. |

---

## 3. Estatística descritiva (anos válidos)

| Variável | n | Média | Mediana | DP | CV % | Mín–Máx | Assim. (EP) | Curt. exc. | MAD |
|---|---|---|---|---|---|---|---|---|---|
| Chuva (mm/ano hidrol.) | 14 | 406,7 | 369,2 | 183,2 | 45,0 | 132,6–861,2 | 0,77 (0,60) | 1,00 | 145,2 |
| Dias ≥ 1 mm | 14 | 33,9 | 34,5 | 11,3 | 33,3 | 16–50 | −0,20 | −1,24 | 14,1 |
| Rx1dia (mm) | 14 | 69,1 | 63,9 | 31,8 | 46,0 | 29,2–135,6 | 0,64 | −0,52 | 33,5 |
| Tméd (°C) | 15 | 27,12 | 26,93 | 0,39 | — | 26,62–27,74 | 0,31 | −1,48 | 0,33 |
| Tmáx (°C) | 15 | 32,82 | 32,72 | 0,51 | — | 32,20–33,67 | 0,42 | −1,32 | 0,60 |
| Tmín (°C) | 15 | 22,29 | 22,30 | 0,27 | — | 21,79–22,73 | −0,14 | −0,87 | 0,26 |
| Amplitude diária (°C) | 15 | 10,53 | 10,47 | 0,35 | — | 9,84–11,10 | 0,17 | −0,38 | 0,26 |
| UR (%) | 15 | 53,4 | 54,6 | 3,3 | 6,1 | 48,3–58,1 | −0,19 | −1,45 | 3,6 |
| TXx (°C) | 15 | 38,5 | 38,5 | 0,79 | — | 37,4–40,1 | 0,31 | −0,79 | 0,89 |
| TNn (°C) | 15 | 17,8 | 17,8 | 0,46 | — | 16,9–18,5 | −0,33 | −0,71 | 0,59 |

- O CV é omitido para °C porque é escala intervalar: o CV mudaria com a origem da escala.
- A chuva anual tem assimetria positiva (g₁ = 0,77, Bowley = 0,67), puxada por 2003/04 (861 mm). Com n = 14, porém, o EP da assimetria (0,60) impede afirmar assimetria com segurança; a média aparada (392 mm) está próxima da média.
- Na escala **diária**, a chuva é extremamente assimétrica (g₁ = 10,9): a mediana diária é 0 mm.
- Ciclo anual (fig09): nov–abr concentra ≈ 85 % da chuva (jan 74, fev 73, mar 50, abr 60 mm; ago–set < 2 mm).

---

## 4. Análise temporal e tendências

### 4.1 Séries anuais (tabela `04_tendencias.csv`, figs 03–05)

| Variável | Theil-Sen /década [IC 95 %] | MK p | MK-TFPW p | MQO /década [IC 95 %] | q (BH) |
|---|---|---|---|---|---|
| Chuva (mm) | −79 [−308; +67] | 0,32 | 0,54 | −137 [−312; +38] | 0,45 |
| Tméd (°C) | +0,26 [−0,16; +0,72] | 0,14 | 0,47 | +0,27 [−0,09; +0,63] | 0,28 |
| Tmáx (°C) | +0,30 [−0,03; +0,92] | 0,075 | 0,59 | +0,43 [−0,03; +0,89] | 0,28 |
| Tmín (°C) | +0,17 [−0,11; +0,41] | 0,37 | 0,47 | +0,19 [−0,07; +0,44] | 0,47 |
| UR (%) | −3,5 [−6,0; +0,6] | 0,092 | 0,59 | −2,6 [−5,6; +0,4] | 0,28 |
| Vento (m/s) | −0,29 [−0,58; −0,08] | 0,023 | 0,049 | — | 0,28 |

**Resultado descritivo.** Todas as estimativas pontuais de temperatura são positivas (+0,2 a +0,4 °C/década), e as de chuva e UR são negativas.

**Evidência estatística.** Na escala anual, **nenhuma** tendência é significativa após a correção por multiplicidade (todos os q ≥ 0,28). Os IC incluem zero. Com n = 14–15, o poder é baixo: para detectar +0,3 °C/década seria preciso um DP interanual bem menor que o observado.

### 4.2 Anomalias mensais (n = 225 meses; EP HAC) — fig06

| Variável | Tendência/década [IC 95 % HAC] | p | q (BH) |
|---|---|---|---|
| Tmáx | **+0,45 [0,18; 0,72] °C** | 0,001 | 0,004 |
| Tméd | **+0,29 [0,06; 0,52] °C** | 0,014 | 0,020 |
| Tmín | **+0,21 [0,05; 0,37] °C** | 0,012 | 0,020 |
| Amplitude diária | **+0,24 [0,08; 0,41] °C** | 0,004 | 0,009 |
| UR | −1,5 [−3,4; +0,3] p.p. | 0,10 | 0,12 |
| Vento | −0,36 [−0,56; −0,17] m/s | < 0,001 | 0,002 |

**Evidência estatística.** Usando toda a informação mensal, o aquecimento é significativo e robusto à autocorrelação. É mais forte na máxima (+0,45 °C/década) que na mínima (+0,21), e por isso a amplitude diária cresce.

**Ressalva crítica: homogeneidade.** Pettitt aponta mudanças na pressão (após 2013; p = 0,043), e a radiação tem uma quebra instrumental documentada em 2014. Um degrau instrumental por volta de 2014 produziria exatamente esse padrão de "tendência" numa regressão linear. **Sem estação de referência para homogeneização relativa, não é possível separar sinal climático de artefato instrumental.** O resultado deve ser tratado como *compatível com* aquecimento, não como prova dele.

### 4.3 Pontos de mudança (Pettitt)
Nenhuma variável climática principal tem mudança significativa: chuva p = 0,36; Tmáx p = 0,11 (após 2011); UR p = 0,08. Pressão (p = 0,043, após 2013) e vento (p = 0,11) têm padrões sugestivos de efeito instrumental.

---

## 5. Distribuição da chuva anual (figs 07–08)

| Distribuição | Parâmetros (MV) | log L | AIC | BIC | ΔAIC | Peso Akaike | A² | p AD (bootstrap) |
|---|---|---|---|---|---|---|---|---|
| Gamma | k = 4,98; θ = 81,7 | −91,76 | 187,5 | 188,8 | 0,00 | 0,32 | 0,43 | 0,32 |
| Weibull | k = 2,44; λ = 458,8 | −91,84 | 187,7 | 189,0 | 0,17 | 0,30 | 0,40 | 0,40 |
| Lognormal | μ = 5,90; σ = 0,50 | −92,27 | 188,5 | 189,8 | 1,03 | 0,19 | 0,58 | 0,10 |
| Normal | μ = 406,7; σ = 176,5 | −92,29 | 188,6 | 189,9 | 1,08 | 0,19 | 0,41 | 0,34 |

**Interpretação.**
- Todos os ΔAIC são < 2 e os pesos de Akaike ficam entre 0,19 e 0,32, ou seja, **os dados não discriminam entre as quatro candidatas**. Nenhuma é rejeitada pelo AD.
- Nos Q-Q plots, a Lognormal superestima a cauda superior e a Normal atribui probabilidade positiva a totais negativos.
- **Escolha: Gamma.** O critério é físico e parcimonioso (suporte positivo, assimetria moderada, padrão do SPI), não uma superioridade estatística demonstrada. A Weibull é equivalente.
- Shapiro-Wilk não rejeita normalidade para nenhuma série anual (todos os q ≥ 0,35), mas com n ≈ 15 o teste tem pouco poder.

---

## 6. Variabilidade e extremos

### 6.1 Anos anômalos (fig12)

| Ano | Indicador | Classe |
|---|---|---|
| 2003/04 | SPI = **+2,05** (861 mm) | extremamente chuvoso (a única observação fora da cerca de Tukey, 726 mm) |
| 2016/17 | SPI = **−1,95** (133 mm) | severamente seco |
| 2011/12 | SPI = **−1,77** (148 mm) | severamente seco |
| 2012/13 | ≥ 69 mm observados, com abril ausente | provavelmente o mais seco, mas excluído |
| 2019 | z(Tméd) = +1,59; TXx = 40,1 °C (máximo da série) | ano mais quente |
| 2015 | z(Tmáx) = +1,67 | máxima média mais alta |
| 2007 | z(Tmín) = −1,84 | noites mais frias |

Os anos quentes (2012, 2013, 2015, 2016, 2019) coincidem com anos secos ou pós-secos. Isso é coerente com a forte relação chuva × temperatura (seção 7).

### 6.2 Chuva diária extrema: POT/GPD (fig14)
- **Limiar**: u = 29,8 mm (P90 dos dias úmidos). A vida residual média é aproximadamente linear até ~36 mm, e ξ̂ fica estável entre 0,1 e 0,2 na faixa de 13–37 mm.
- **Dados**: 59 picos desagrupados (corridas, r = 1 dia); 18,5 anos climatológicos efetivos; λ = 3,2 eventos/ano.
- **GPD**: σ̂ = 18,0 mm (EP 3,8); ξ̂ = 0,19 (EP 0,16). A cauda é possivelmente pesada, mas ξ = 0 não é rejeitado (IC ≈ [−0,13; 0,51]).

| T (anos) | Nível (mm/dia) | IC 95 % bootstrap |
|---|---|---|
| 2 | 70 | 57–85 |
| 5 | 95 | 73–120 |
| 10 | 118 | 84–157 |
| 20 | 144 | 95–207 |
| 50 | 184 | 105–302 (**extrapolação**) |

O maior evento diário, 159,8 mm em 18/03/2008, tem período de retorno empírico de ~20 anos. Os IC se alargam rapidamente. **Valores acima de T ≈ 20 anos são extrapolação** de um registro de 18,5 anos efetivos e não devem embasar projeto hidráulico sem regionalização.

### 6.3 Frequência de extremos (fig15)

| Índice | Razão/década [IC 95 %] | p | Dispersão φ̂ |
|---|---|---|---|
| R95 (dias > 40,4 mm) | 0,80 [0,32; 2,00] | 0,60 | 1,4 |
| Dias ≥ 1 mm | 0,88 [0,61; 1,26] | 0,45 | 4,3 |
| TX90p (dias quentes) | 1,66 [0,94; 2,94] | 0,075 | 9,4 |
| TN10p (noites frias) | 0,79 [0,48; 1,30] | 0,32 | 7,3 |

Os sinais têm a direção esperada num cenário de aquecimento (mais dias quentes, menos noites frias), mas não são significativos (q ≥ 0,30). A forte sobredispersão indica agrupamento interanual, e daí o uso de quase-verossimilhança.

---

## 7. Relações entre variáveis (figs 10–11)

### Escala anual (anos hidrológicos, n = 14–16, q de BH)
- **Chuva × Tmáx**: r = −0,82 [−0,94; −0,52], ρ = −0,88, q < 0,001. Anos mais chuvosos têm máximas mais baixas, cerca de −0,25 °C por 100 mm.
- **Chuva × UR**: r = +0,83 [0,53; 0,94]. **Chuva × Td**: r = +0,81.
- **Tméd × UR**: r = −0,94. **UR × Td**: r = 0,99 (UR e Td contêm praticamente a mesma informação interanual).
- **Rx1dia × chuva anual**: ρ = 0,59. Um único evento pesa no total anual.
- **Vento**: sem relação significativa com as demais variáveis.

Grupos por SPI: anos secos (SPI < 0) são **0,86 °C mais quentes** na Tmáx do que os úmidos [IC Welch 0,48; 1,23], com g = 2,3 e δ de Cliff = 1,0 (separação total). O teste é p = 0,0004 (q = 0,003), mas tem **só 6 vs. 8 anos**.

### Escala mensal (anomalias, n = 196–225, n efetivo de Bretherton)
- **UR × Tmáx**: r = −0,84 [−0,88; −0,79], n_ef ≈ 148.
- **Tmáx × Tmín**: r = 0,77 [0,70; 0,83]; inclinação de 0,38 °C por °C.
- **Chuva × UR**: ρ = 0,56.
- **Regressão múltipla** Tmáx ~ UR + chuva (EP HAC; VIF = 1,75): a UR domina (−0,17 °C por p.p., p < 10⁻³⁰). A chuva **deixa de ser significativa** quando a UR entra no modelo (p = 0,93); R² = 0,76.

**Interpretação.** As relações são fisicamente coerentes: menos chuva, menos umidade do solo e menos nebulosidade, o que significa mais calor sensível e máximas mais altas. A regressão múltipla sugere que o efeito da chuva sobre a Tmáx é mediado pela umidade atmosférica. **Isso é associação, não causalidade demonstrada.** O mecanismo é plausível, mas estes dados não excluem causas comuns, como a circulação de grande escala (ENOS, Atlântico).

---

## 8. Mudança entre períodos

| Variável | Antigo → recente | Δ [IC Welch 95 %] | g de Hedges [IC] | p Welch | p Mann-Whitney |
|---|---|---|---|---|---|
| Chuva | 473 → 357 mm | −117 [−367; +133] | −0,61 [−1,69; 0,48] | 0,31 | 0,28 |
| Tméd | 27,01 → 27,26 °C | +0,25 [−0,18; +0,68] | 0,63 [−0,41; 1,66] | 0,23 | 0,19 |
| Tmáx | 32,62 → 33,04 °C | +0,42 [−0,14; +0,97] | 0,81 [−0,25; 1,86] | 0,13 | 0,12 |
| Tmín | 22,22 → 22,37 °C | +0,15 [−0,15; +0,45] | 0,54 [−0,49; 1,57] | 0,29 | 0,54 |
| UR | 54,6 → 52,1 % | −2,5 [−6,0; +1,0] | −0,74 [−1,79; 0,30] | 0,15 | 0,15 |

Períodos: 2003/04–2012/13 vs. 2013/14–2023/24 para chuva; 2003–2013 vs. 2014–2024 para temperatura. Os efeitos são de tamanho médio a grande (|g| ≈ 0,5–0,8), mas **nenhuma diferença é estatisticamente distinguível de zero**: os IC de g são largos, com n = 6–8 por grupo. Além disso, o corte coincide com a quebra instrumental de 2014, o que torna esta comparação a mais vulnerável a artefatos.

---

## 9. Modelagem e diagnóstico (figs 16–17)

Modelo: $Y_t=\beta_0+\beta_1(t-\bar t)+\varepsilon_t$.

| Modelo | β₁/década (EP) | R² | AIC | Shapiro p | Breusch-Pagan p | DW | r₁ | Ljung-Box p | ρ̂ AR(1) | p LR(ρ = 0) |
|---|---|---|---|---|---|---|---|---|---|---|
| Chuva | −137 (80) | 0,20 | 187,5 | 0,98 | 0,27 | 1,84 | −0,20 | 0,90 | −0,24 | 0,53 |
| Tméd | +0,27 (0,17) | 0,17 | 16,2 | 0,06 | **0,011** | 1,89 | −0,14 | 0,94 | −0,15 | 0,65 |
| Tmáx | +0,43 (0,21) | 0,24 | 23,5 | 0,35 | 0,061 | 1,87 | −0,11 | 0,77 | −0,10 | 0,74 |
| Tmín | +0,19 (0,12) | 0,16 | 5,8 | 0,30 | 0,75 | 2,76 | −0,39 | 0,27 | −0,30 | 0,27 |
| UR | −2,6 (1,4) | 0,22 | 79,5 | 0,86 | 0,077 | 1,22 | +0,38 | 0,45 | +0,35 | 0,22 |

- **Normalidade**: aceitável em todos os modelos.
- **Heterocedasticidade**: indício em Tméd (BP p = 0,011), ligado à maior dispersão dos anos recentes, que também é compatível com uma quebra de nível. As conclusões não mudam: o IC já inclui zero.
- **Autocorrelação**: sem evidência na escala anual (LR para ρ = 0 com p ≥ 0,22). O MQG-AR(1) tem AIC **maior** que o MQO em todos os casos, então o MQO simples é adequado. A UR tem r₁ = 0,38 (DW = 1,22), mas não significativo com n = 15. A estimativa AR(1) (−2,5 %/década, p = 0,16) confirma a ausência de significância.
- **Chuva, GLM Gamma com ligação log**: −24,8 %/década [−50; +14], p = 0,16. O AIC (187,2) é praticamente igual ao do MQO; os dois modelos concordam.
- **R² baixo (0,16–0,24)**: a variabilidade interanual domina amplamente qualquer componente linear. Modelos mais complexos (ARIMA, polinômios) não se justificam com 14–15 pontos.

---

## 10. Síntese: descritivo → evidência → interpretação

| Achado | Robustez |
|---|---|
| Clima semiárido: 407 ± 183 mm/ano (CV 45 %), 85 % em nov–abr; Tméd 27,1 °C; UR 53 % | **Robusto** (descritivo) |
| Forte acoplamento interanual chuva ↔ UR ↔ Tmáx (\|r\| ≈ 0,8–0,9); anos secos ≈ 0,9 °C mais quentes | **Robusto**: sobrevive a BH, aparece em Pearson e Spearman, nas escalas anual e mensal, e é fisicamente coerente |
| Aquecimento de +0,2 a +0,45 °C/década nas anomalias mensais, maior na Tmáx | **Significativo, mas com incerteza de homogeneidade** (possível troca instrumental em 2014) |
| Tendências anuais de chuva, UR e temperatura | **Não significativas**: estimativas pontuais sugestivas, IC largos |
| Diferença entre as metades da série | **Não significativa**; efeitos médios com IC muito largos |
| Chuva anual ~ Gamma (ou Weibull) | **Incerto**: as quatro candidatas são indistinguíveis com n = 14 |
| Chuva diária de 10 anos ≈ 118 mm [84; 157] | **Moderado** até T ≈ 20 anos; incerteza alta acima disso |
| Queda do vento | **Provavelmente instrumental ou de exposição**, não climática |

**Interpretação climatológica, com cautela.** Os dados são compatíveis com aquecimento regional, especialmente nas máximas, e com uma ligação estreita entre seca e calor no semiárido. Contudo, (i) a série efetiva é curta, (ii) há indícios de inomogeneidade instrumental e (iii) o ano mais seco foi excluído. Com isso, **esta série isolada não permite atribuir as tendências a mudança climática**: uma associação temporal não é evidência de causa.

---

## 11. Limitações

1. **Série curta**: 14 anos válidos de chuva e 15 de temperatura. O poder para tendências é baixo e os extremos anuais não sustentam uma GEV.
2. **Falhas não aleatórias**: concentram-se na estação chuvosa. A exclusão de 2012/13 remove o ano provavelmente mais seco e enviesa a cauda seca e a variância.
3. **Homogeneidade não verificável**: sem estação vizinha, não há como testar homogeneidade relativa (SNHT, RHtest). Há indícios de trocas por volta de 2008 e 2014.
4. **Radiação descartada**: por isso não foi possível estimar a evapotranspiração (Penman-Monteith).
5. **Critérios de completude são escolhas**: FCA 5 %, regra OMM, ≥ 10 meses. 2005/06 (FCA 5,2 %) está no limite. Todos os limiares estão em `src/00_config.jl`.
6. **Percentis relativos ao próprio período** (TX90p, P95): não são comparáveis com uma normal climatológica oficial (1991–2020).
7. **Testes de grupos** (Welch, Mann-Whitney) supõem independência; com autocorrelação positiva, os p-valores seriam otimistas. Na escala anual, o r₁ observado é pequeno.

## 12. Análises adicionais justificadas pelos dados

1. **Homogeneização relativa** com estações próximas: a convencional do INMET em Petrolina, estações de Juazeiro-BA, ou reanálise ERA5-Land. É o passo mais importante antes de qualquer afirmação sobre tendência.
2. **Preencher 2012/13, 2017/18 e 2022–24** com pluviômetros vizinhos (ANA/APAC/CEMADEN), para recuperar anos secos críticos.
3. **Covariáveis de grande escala** (Niño 3.4, gradiente do Atlântico tropical) numa regressão da chuva e da Tmáx, para separar variabilidade natural de tendência.
4. **Análise sazonal** (Seasonal Kendall por mês ou estação): verificar se o aquecimento se concentra na estação seca (set–nov), onde a Tmáx é maior.
5. **Estimar ETo pela equação de Hargreaves**, que só exige Tmáx, Tmín e radiação extraterrestre calculada, contornando a radiação medida defeituosa.

---

### Índice de figuras
| Arquivo | Conteúdo |
|---|---|
| fig01_completude | Horas ausentes por ano × mês (chuva e temperatura) |
| fig02_radiacao_inconsistencia | Justificativa para excluir a radiação |
| fig03_precipitacao_anual | Série da chuva anual com SPI, média móvel e anos excluídos |
| fig04_temperatura_anual | Tmáx, Tméd e Tmín anuais com média móvel |
| fig05_tendencias_ic | Tendência MQO ± IC 95 % e Theil-Sen, 6 variáveis |
| fig06_anomalias_mensais_temperatura | Anomalias mensais com tendência HAC e marca de 2014 |
| fig07_distribuicao_precipitacao | Histograma, densidades e ECDF com os modelos ajustados |
| fig08_qqplots_precipitacao | Q-Q das 4 distribuições candidatas |
| fig09_boxplots_climatologia_mensal | Ciclo anual (boxplots mensais) |
| fig10_matriz_correlacao | Matriz de Spearman (anos hidrológicos) com marcação BH |
| fig11_dispersao_regressoes | Dispersões com retas de regressão e IC |
| fig12_anos_extremos | SPI e anomalias padronizadas de temperatura |
| fig13_comparacao_periodos | Período antigo × recente |
| fig14_pot_gpd_retorno | Escolha do limiar, Q-Q da GPD e curva de retorno |
| fig15_indices_extremos_contagem | TX90p e R95 com ajuste GLM |
| fig16/17_residuos_* | Diagnóstico dos modelos de tendência |
