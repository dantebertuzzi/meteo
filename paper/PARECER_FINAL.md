# Parecer final de revisão por pares (simulada) — 17/10/2026

**Manuscrito:** *Menos chuva e mais calor em Petrolina-PE (1961–2025)…*, 24 p., commit 547fe4e.
**Escopo:** leitura integral do texto; testes do código; recálculo independente dos números
principais e de números fora da verificação automática; análise visual das 12 figuras;
interpretação dos resultados; nota.

> **Aviso de independência.** Este parecer, como os quatro anteriores, foi feito pelo mesmo
> assistente que ajudou a escrever o artigo. As checagens numéricas são objetivas, mas o julgamento
> de mérito e a nota não substituem revisores humanos independentes.

---

## 1. Testes e verificação dos cálculos

| Verificação | Resultado |
|---|---|
| Testes unitários Julia (`test/runtests.jl`) | **73/73** passam, em 16 grupos: FAO-56, SPI, PWM, corridas, perfil e bootstrap, SNHT, Markov, Alexandersson–Moberg, ETCCDI, Hamed–Rao e TMD, STL, ESS, leitura do BDMEP, log-normal, não estacionário e Poisson-binomial |
| Verificação independente em Python (`verifica.py`) | **57/57** números conferem; nenhuma inconsistência interna (estimativas dentro dos IC, probabilidades crescentes com n) |
| Índices ETCCDI de chuva (fora do `verifica.py`) | Recalculados dos dados brutos: R1mm, Rx1day, SDII e R95p idênticos ao Julia (diferença máxima de 6e-14); declives de Sen e IC idênticos aos da Tabela 6 |
| Razão A307/82983 | Média 0,960, mínimo 0,59, máximo 1,12, como no texto |
| Normais oficiais do INMET | 419 / 483 / 577 mm conferidos contra as planilhas originais |
| ETo e Piche | Recalculadas por outro método (média simples): mesmas conclusões |
| Dado suspeito | O máximo diário do registro, 151,3 mm em 01/01/1978, coincide em valor com o total de 2012/13 (151,3 mm em 25 dias, sem falhas). Os dias anteriores têm chuva (35,5 mm em 28/12), então não parece valor acumulado. Sem esse ano, o Rx1day passa de −7,2 para −6,6 mm/década e o R95p de −41 para −40,5. **Robusto.** |

**Conclusão:** não encontrei erro de cálculo. A cadeia dados brutos → números do texto é
reprodutível e foi conferida por duas implementações independentes nos pontos centrais.

## 2. Análise das figuras

| Fig. | Avaliação | Pedido |
|---|---|---|
| 1 (área) | Clara; os dois painéis e o anel de distância ajudam. Juazeiro quase sobrepõe Petrolina (3 km), o que é inevitável. | — |
| 2 (ausências) | Boa leitura das lacunas (1966–70, 1985–91) e do corte de temperatura antes de 1971. | — |
| 3 (padrões de ausência) | Correta. No painel (b), três padrões com 0,0 % ocupam espaço sem informação. | **F1** (menor) |
| 4 (chuva anual) | Eficaz; destaca os três anos secos e as duas tendências. A linha não estacionária vai até cerca de 2024, além do último ano válido (2020/21), sem indicar que esse trecho é extrapolação. | **F2** (menor) |
| 5 (SPI/SPEI) | Clara; as lacunas são visíveis, e o SPEI começa em 1972 por causa do corte da temperatura. Mostra bem o argumento de que os eventos recentes refletem a cobertura. | — |
| 6 (período de retorno) | Os pontos observados (posição de Weibull) são a distribuição **marginal** de anos de uma série com tendência e por isso seguem a curva estacionária. Um leitor apressado conclui que o modelo estacionário se ajusta bem. A legenda não avisa. | **F3** (menor) |
| 7 (posterior de T) | Correta, mas mostra só a Gamma bayesiana, que é **sensibilidade**. O resultado principal (média de Akaike) não tem figura de incerteza. | **F4** (menor) |
| 8 (Markov) | Clara; a comparação (a) × (b) mostra o argumento central da persistência. | — |
| 9 (climatologia) | Boa; a legenda já explica a diferença com a normal oficial. | — |
| 10 (estação chuvosa) | Correta; a lacuna dos anos sem estação depois de 2010 é visível e está discutida no texto. Os inícios muito precoces de 1970–73 pesam no declive, o que é coerente com o p do TFPW não significativo. | — |
| 11 (tendências) | No painel (c), o **R95p aparece com −28 %/década**. Linearmente, isso daria −160 % da média ao longo do período, o que é impossível. O índice tem 9 zeros em 47 anos e é muito assimétrico: a média dos 10 primeiros anos é 229 mm e a dos 10 últimos, 59 mm (−74 %). A queda é real, mas a escala "% da média por década" exagera e não é interpretável para esse índice. | **R5.1** |
| 12 (quantílica) | Clara e honesta: os intervalos se sobrepõem, e o texto não força desigualdade entre os quantis. | — |

Numeração conforme o PDF: os arquivos `fig03_climatologia` e `fig12_padroes_ausencia` aparecem no
texto como Figuras 9 e 3.

**Qualidade gráfica geral:** boa. As figuras são vetoriais, com vírgula decimal, legendas
completas, paleta consistente, rótulos diretos e sem duplo eixo y. Estão acima da média da
literatura regional.

## 3. Interpretação dos resultados

**P2 (regime): bem sustentada.**
- A queda da chuva, −60 [−87; −33] mm/década, é:
  - confirmada pela referência regional, por 4 de 5 postos da ANA e pelas normais oficiais;
  - coerente com Bebedouro (−40 a −68) e com o Submédio;
  - estável em magnitude nas subamostras.
- A dependência da década de 2010 está declarada no resumo.
- O aquecimento da Tmín é o resultado de temperatura mais firme; a faixa da Tmáx (0,19–0,30) está
  honestamente apresentada.
- As ressalvas sobre dias chuvosos e evaporação (R4) estão aplicadas.

**P1 (frequência de anos secos): bem delimitada.** O texto não promete mais do que entrega:
- T = 12 [5,0; 88] anos é apresentado como ordem de grandeza;
- a sensibilidade à forma, à inclinação regional e ao fim da série é reportada;
- o viés de seleção de E1 está declarado;
- a validação fora da amostra sustenta a não estacionariedade com um limiar não viesado.

A conclusão "a incerteza é parte do resultado" é a leitura correta.

**Persistência e teleconexões:** os resultados nulos são interpretados como falta de poder, e não
como ausência de efeito, o que está correto.

**Ponto de atenção remanescente.** O valor "no clima de 2020/21" é a média ajustada no extremo da
série, onde a tendência é menos precisa. O IC largo reflete isso, e a sensibilidade com FCA ≤ 7,5 %
mitiga. Não exige mudança.

## 4. Pedidos desta rodada

| Ponto | Gravidade | Pedido |
|---|---|---|
| **R5.1** | Menor-moderada | R95p na Fig. 11c: trocar "% da média por década" por outra escala para índices com zeros (por exemplo, mm/década com eixo próprio, ou variação entre o primeiro e o último terço), ou anotar que a escala relativa não se aplica ao R95p. |
| **F3** | Menor | Legenda da Fig. 6: "os pontos são a distribuição marginal de todos os anos e só são comparáveis à curva estacionária". |
| **F2** | Menor | Fig. 4: terminar a linha não estacionária em 2020/21, ou tracejar o trecho extrapolado. |
| **F4** | Menor | Fig. 7: dizer na legenda que é a sensibilidade Gamma, ou levá-la para um apêndice. Opcional: histograma bootstrap da média de Akaike. |
| **F1** | Cosmético | Fig. 3b: mostrar só os padrões com fração > 0,05 %. |

Nenhum desses pontos muda números ou conclusões.

## 5. Limitações estruturais (não corrigíveis neste trabalho)

1. **Uma estação para a P1.** O período de retorno vale para um ponto. A extensão regional depende
   de dados em grade ou dos postos da ANA.
2. **Muitas escolhas de análise depois de ver os dados.** O estimador principal mudou durante a
   revisão (de Gamma para média de Akaike) e há muitas sensibilidades. Isso foi feito às claras e
   todas as alternativas estão na Tabela 4, mas um revisor externo pode levantar o risco de
   resultados dependentes de escolhas.
3. **Originalidade moderada na P2.** O sinal das tendências já era conhecido (Pimentel e Assis
   2022; Bebedouro; Assis et al. 2022). A contribuição é de rigor e quantificação, não de
   descoberta. A P1 é a parte nova.
4. **Revisão não independente**, como no aviso do início.

## 6. Veredito e nota

**Veredito: aceitar com revisão menor** (R5.1 e ajustes de figura). Como preprint no Zenodo, está
pronto para depósito mesmo sem esses ajustes.

| Critério | Peso | Nota (0–10) | Comentário |
|---|---|---|---|
| Rigor estatístico | 25 % | 9,0 | IC, TFPW, BH, TMD, média de Akaike, validação fora da amostra, viés de seleção declarado |
| Qualidade e tratamento dos dados | 15 % | 9,0 | Homogeneidade contra A307, vizinhas, ANA e normais; FCA; único instrumento |
| Reprodutibilidade | 15 % | 10 | Todos os números por macro; 73 testes; verificação independente 57/57; código aberto |
| Interpretação e honestidade | 15 % | 9,0 | Limites declarados no resumo; nulos tratados como falta de poder |
| Originalidade e contribuição | 15 % | 6,5 | P2 confirma o conhecido com mais rigor; P1 é nova, mas de uma estação |
| Clareza e redação | 10 % | 8,0 | Bem organizado; denso em partes (§ 5.2); 24 p. é longo para o escopo |
| Figuras | 5 % | 8,0 | Boas, com os ajustes F1–F4 e R5.1 pendentes |

**Nota final: 8,6 / 10** (média ponderada: 8,63).

**Leitura da nota.** Como trabalho de uma estação, o artigo está no topo do que se publica
regionalmente sobre Petrolina. Numa revista nacional da área (por exemplo, Revista Brasileira de
Meteorologia), a expectativa razoável é de aceitação com revisão menor ou moderada. Numa revista
internacional (Theoretical and Applied Climatology, International Journal of Climatology), o
limitante seria a originalidade e o escopo de um ponto: a resposta provável é de revisão maior,
com pedido de extensão regional (postos da ANA ou grade). O que eleva o trabalho é o rigor e a
transparência, e não a novidade do resultado.

---

## Resolução — 17/10/2026

- **R5.1:** o R95p saiu do painel relativo e ganhou o painel (d), em mm/década (−41 [−57; −21]).
  A legenda da Fig. 11 explica o motivo. Os painéis têm altura proporcional ao número de linhas.
- **F2:** na Fig. 4, as médias não estacionárias terminam no último ano válido (2020/21), e a
  legenda diz isso.
- **F3:** a legenda da Fig. 6 explica que os pontos são a distribuição marginal, comparável só à
  curva estacionária.
- **F4:** a legenda da Fig. 7 começa com "Sensibilidade:" e diz que são modelos Gamma bayesianos.
- **F1:** a Fig. 3 mostra só padrões com ≥ 0,05 % dos dias; o painel da A307 fica com três.

Nenhum número mudou: `verifica.py` 57/57, artigo com 24 p. Com os ajustes, o critério "Figuras"
passa de 8,0 para 8,5, e a nota ponderada vai de 8,63 para 8,65, o que **mantém a nota final em
8,6 / 10**.
