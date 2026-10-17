# Parecer final de revisão por pares (simulada), 2ª rodada — 17/10/2026

**Manuscrito:** versão com a extensão regional (período de retorno nos postos da ANA e análise
regional de frequência), os anos chuvosos recentes e o contexto internacional. 28 p., com apêndice.
**Escopo:**
- conteúdo novo;
- resquícios de conteúdo desatualizado;
- figuras: uma nova, distribuição entre texto e apêndice;
- auditoria completa das referências;
- testes e verificação independente;
- nota.

> **Aviso de independência.** Revisão feita pelo mesmo assistente que ajudou a escrever o artigo.
> Testes e recálculos são objetivos; o julgamento de mérito e a nota não substituem revisores
> humanos independentes.

---

## 1. Testes e verificação

| Verificação | Resultado |
|---|---|
| Testes unitários Julia | **81/81**, incluindo 8 novos do modelo regional (recuperação de índices e tendência em dados sintéticos; aninhamento com uma estação e com inclinação por estação) |
| Verificação independente em Python | **70/70**, antes 57/57. Os 13 números novos conferidos: tendência com 2021–23; tendências dos postos desde 1963/64, com e sem os anos recentes; limiares relativos; período de retorno estacionário de três postos |
| Pipeline completo e modo `figuras` (cache) | Sem erros. O modo cache falhava na figura nova, porque as funções desserializadas não podiam ser chamadas (*world age*); corrigido com `Base.invokelatest` |
| Compilação | 28 p., sem referências indefinidas nem *overfull boxes* |

A análise regional não tem verificação em Python, porque o modelo é complexo demais para
reimplementar com rapidez. Ela é coberta pelos testes unitários com dados sintéticos.

## 2. Conteúdo novo

**Período de retorno nos postos da ANA (Tabela 3).** Correto e bem justificado:
- o limiar é relativo porque o 2012/13 de Juazeiro não é válido;
- o 2012/13 do próprio posto, quando válido, entra como sensibilidade;
- todos os postos começam em 1963/64, o que exclui o segmento inomogêneo do Lagoa Grande.

O período de retorno cai nos 5 postos, de forma coerente com a 82983.

**Anos chuvosos recentes (Tabela 2, última linha).** É o achado mais importante da rodada, e o
texto o trata com honestidade: −43 [−70; −11] mm/década contra os −60 principais, a faixa aparece
no resumo e nas conclusões, e a Limitação 6 diz que a tendência de longo prazo pode ser menor.
Usar o total observado de 2022/23 como limite inferior é conservador, porque subestima esse ano.

**Análise regional (índice de cheia com tendência).** O modelo é adequado, e o bootstrap de anos
inteiros é a escolha certa para estações correlacionadas. Três pontos estão declarados no texto:
- a verossimilhança combinada superestima a evidência;
- os intervalos estreitam pouco;
- a tendência é heterogênea entre estações (−11,5 a −0,6 %/década).

Ressalva menor: os pesos de Akaike vêm da verossimilhança combinada e por isso são mais
concentrados do que deveriam. Como o IC vem do bootstrap, que reajusta os pesos em cada réplica,
a incerteza final está bem representada. Não exige mudança.

**Contexto internacional.**
- As afirmações foram conferidas nos resumos, e não há porcentagem sem fonte (Sahel).
- As três lições são pertinentes.
- A frase nos Métodos sobre 1/pₜ e o tempo de espera (Salas e Obeysekera 2014) melhora o rigor
  conceitual.

## 3. Resquícios desatualizados encontrados e corrigidos

| Onde | Antes | Agora |
|---|---|---|
| Resumo PT/EN | "com uma estação e cerca de meio século…" | "mesmo juntando estações vizinhas, cujos anos secos coincidem…" |
| Introdução (método) | comparação "com estações vizinhas" | inclui pluviômetros da ANA e a análise regional |
| Introdução (fim) | "o que uma única estação consegue dizer" | "uma estação, e mesmo um conjunto de estações vizinhas" |
| Discussão (seca 2011–17) | "passa à ordem de uma década" | "uma ou duas décadas, na 82983, nos postos e na análise regional" |
| Discussão (limites) | "mais poder exige mais estações ou dados em grade" | estações vizinhas pouco ajudam; é preciso estações distantes, séries longas ou grade |
| Limitação 1 | só o tamanho da amostra | acrescenta que a análise regional estreita pouco |
| README | "compara 34 números"; descrição antiga | 70 números; descrição com ANA, normais e análise regional |

Os trechos já atualizados nas rodadas anteriores (Limitações 6 e 9, Conclusões 1 e 4, Tabelas 2 e
3) conferem com os números atuais.

## 4. Figuras

- **Nova Figura 5:** período de retorno por estação e regional, estacionário (vazio) e não
  estacionário (cheio), com IC. Mostra numa imagem o principal resultado novo: a direção se repete
  em todas as estações e os intervalos se sobrepõem. A legenda foi levada para fora da área do
  gráfico.
- **Redistribuição:** três figuras de apoio foram para o **Apêndice A**:
  - padrões de ausência (A1);
  - posterior da Gamma bayesiana, que é sensibilidade (A2);
  - regressão quantílica, cujos números já estão no texto (A3).

  O texto principal passa de 12 para **10 figuras**, e as referências cruzadas se ajustam sozinhas.
- Para uma revista, as Figuras 2 (ausências) e 9 (estação chuvosa) também poderiam ir para o
  suplemento. No preprint, ficam.

## 5. Referências

- **60 entradas, 60 citadas.** `wilks2019` não era citada e foi removida. `alexandersson1986` não
  era citada, embora o SNHT seja usado várias vezes; agora é citada no primeiro uso.
- **46 DOIs conferidos** contra o Crossref (título e ano): todos corretos. Burnham e Anderson
  (2002) aparece como 2004 no Crossref por causa da data da versão eletrônica; 2002 é a 2ª edição
  impressa, e está correto.
- **2 DOIs do Zenodo** (DataCite) resolvem: MissingPatterns.jl e Gois et al. 2024. A versão citada
  do MissingPatterns.jl (0.6.0) é a usada no `Manifest.toml`.
- **12 entradas sem DOI** (livros, relatórios, anais, JMLR, PMLR, página do INMET) estão corretas,
  porque essas obras não têm DOI.
- Citação "Dijk et al." corrigida para "van Dijk et al." (`useprefix`).

## 6. Pontos remanescentes

| Ponto | Gravidade | Comentário |
|---|---|---|
| Tamanho: 28 p. com apêndice | Para revista | Adequado como preprint. Para submissão: 15–18 p., com suplemento |
| Estimativa principal −60 × −43 | Decisão do autor | −60 segue o critério definido antes de ver os dados; −43 está ao lado em todos os lugares relevantes |
| Análise regional sem verificação em Python | Menor | Coberta pelos testes unitários |
| Nota não independente | Estrutural | Ver aviso |

Nenhum erro de cálculo ou de referência pendente.

## 7. Veredito e nota

**Veredito: aceitar** como preprint, sem pendências. Para revista, aceitar com revisão menor
(encurtar).

| Critério | Peso | Antes (R1) | Agora | Por quê |
|---|---|---|---|---|
| Rigor estatístico | 25 % | 9,0 | 9,0 | Mantido; análise regional com bootstrap correto |
| Dados e homogeneidade | 15 % | 9,0 | 9,0 | Segmento inomogêneo do Lagoa Grande detectado e excluído |
| Reprodutibilidade | 15 % | 10 | 10 | 81 testes; verificação 70/70 |
| Interpretação e honestidade | 15 % | 9,0 | 9,5 | Achado dos anos recentes reportado contra o número de destaque |
| Originalidade e contribuição | 15 % | 6,5 | 7,5 | Período de retorno regional não estacionário; contexto internacional |
| Clareza e redação | 10 % | 8,0 | 7,5 | Mais conteúdo e mais denso; 28 p. |
| Figuras | 5 % | 8,5 | 9,0 | Figura nova informativa; texto principal com 10 figuras |

**Nota final: 8,9 / 10** (média ponderada: 8,85).

**Leitura:** o trabalho deixou de ser a análise de uma estação e passou a ter uma verificação
regional. O achado dos anos recentes tornou o resultado principal mais modesto e mais difícil de
contestar. Para uma revista internacional, a chance passa de moderada para boa, desde que o texto
seja encurtado.
