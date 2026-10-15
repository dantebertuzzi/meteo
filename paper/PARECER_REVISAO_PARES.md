# Parecer de revisão por pares (simulada) — 15/10/2026

**Manuscrito:** *Menos chuva e mais calor em Petrolina-PE (1961–2025): tendências do regime
meteorológico e o que elas implicam para a frequência de anos muito secos* (19 p.)

**Recomendação:** **revisão maior** (para periódico). Como preprint no Zenodo, é publicável
depois de corrigir M1 e M2, que afetam números do resumo.

Este parecer foi feito com o texto atual, com as macros de `paper/resultados.tex` e com
reanálises independentes em Python/SciPy (código abaixo de cada ponto, quando há).

---

## 1. Resumo do manuscrito

O trabalho usa a estação convencional INMET 82983 (47 anos hidrológicos válidos, 1963/64–2020/21)
para estimar tendências de chuva e temperatura em Petrolina e o período de retorno de um ano tão
seco quanto 2012/13, sob modelos estacionários e não estacionários. A homogeneidade é avaliada
contra a A307 (mesmo sítio) e contra quatro vizinhas.

## 2. Avaliação geral

**Pontos fortes**

- **Reprodutibilidade exemplar para o porte do trabalho.** Todos os números vêm de macros
  geradas pelo código; a verificação independente em Python confere **34/34** valores, e os testes
  unitários Julia passam. Isso é raro e deve ser mantido.
- Limitações declaradas com franqueza; o reenquadramento (regime primeiro, período de retorno
  como implicação, incerteza como resultado) é adequado ao que os dados permitem.
- Bons cuidados metodológicos: critério de completude ponderado pela climatologia (FCA),
  TFPW só com pares realmente consecutivos, tendência mínima detectável para resultados nulos,
  persistência testada com e sem tendência, sensibilidade à inclinação regional.

**Problemas principais:** a escolha da Gamma no modelo não estacionário determina o número
central do resumo (M1); há quebras na Tmáx frente às vizinhas que não aparecem nos Resultados (M2);
e a significância da tendência da chuva depende da década final (M3).

---

## 3. Pontos principais (exigem mudança)

### M1. O período de retorno não estacionário depende da distribuição, e a Gamma é a pior candidata

> **Resolvido em 15/10/2026.** `ns_media_akaike` (src/paper/retorno.jl) ajusta as quatro candidatas
> com posição log-linear; resultado principal = média por pesos de Akaike com bootstrap dos anos:
> **12 anos [6,3; 91]** em 2020/21 (antes: Gamma bayesiana, 31). Inclinação regional: 17 [11; 76].
> Gamma mantida como sensibilidade. Testes unitários e `verifica.py` (39/39) cobrem o novo cálculo.
> Obs.: a Gamma é a 2ª pior (a log-normal é a pior); a tabela acima foi feita sem a log-normal.

No caso estacionário, o artigo reconhece que nenhuma distribuição domina e usa a média por pesos
de Akaike (a Gamma tem peso 0,13). No não estacionário, porém, só a Gamma é ajustada. Reajustando
com a mesma estrutura (`log μ` ou média linear no tempo), na máxima verossimilhança:

| Modelo não estacionário | AIC | T de E1 em 2020/21 (anos) |
|---|---|---|
| Normal (média linear) | **609,5** | **8** |
| Weibull (log-escala linear) | 612,2 | 26 |
| Gamma (log-média linear; usado no artigo) | 614,5 | 36 |

A Gamma é a candidata com **pior** AIC, e o resultado varia de 8 a 36 anos conforme a
distribuição. Isso é coerente com o que o próprio texto admite: a cauda seca da Gamma é mais fina
que a dos dados (T empírico ≈ 24 contra 81 na Gamma estacionária).

Há também uma comparação imprópria no resumo e na conclusão 4: "de 38 anos (estacionário) para
31 (não estacionário)" compara uma média de quatro distribuições com uma Gamma. Na mesma
distribuição, a Gamma passa de **81 para 36 anos** (MV), uma mudança muito maior do que o texto
sugere.

**Pedido:** ajustar as mesmas candidatas no caso não estacionário e reportar a média por pesos de
Akaike (com bootstrap), como no estacionário; manter a Gamma bayesiana como sensibilidade.
Comparar estacionário e não estacionário dentro do mesmo esquema. A normal admite valores
negativos; isso é irrelevante na cauda de 151 mm, mas convém dizer.

### M2. Quebras na Tmáx frente às vizinhas no período analisado

`results/paper/homogeneidade_temperatura.csv` mostra, **depois de 1971**:

| Teste | Quebra | Salto | q (BH) |
|---|---|---|---|
| Tmáx × Cabrobó | 2000 | 0,27 °C | **0,043** |
| Tmáx × Remanso | 1974 | 1,93 °C | **0,001** |
| Tmín × Cabrobó | 1987 | — | 0,10 |
| Tmín × Remanso | 1991 | — | 0,10 |

O texto só reporta a comparação com vizinhas antes de 1971. Um salto de 0,27 °C em 2000 é da
mesma ordem da tendência total da Tmáx por década (0,28 °C) e fica **antes** do corte de 2008
usado como sensibilidade, que portanto não o elimina. A quebra de 1974 contra Remanso, no início
do período, pode ser da própria Remanso, mas precisa ser discutida.

> Nota: na rodada anterior, a limitação 4 foi alterada para dizer que a quebra frente às vizinhas
> se referia só ao período anterior a 1971. Isso estava errado e já foi corrigido em `main.tex`
> (agora usa as macros `\NquebrasTemp` e `\NtestesTemp`).

**Pedido:** reportar esses testes em Resultados; estimar a tendência da Tmáx com a quebra de 2000
ajustada (ou em 1971–1999 e 2001–2021 separadamente) e apresentar a Tmín, sem quebra
significativa frente às vizinhas, como o resultado de temperatura mais firme.

### M3. A significância da tendência da chuva depende da década de 2010

Médias por década dos anos válidos: 1970s ≈ 611 mm, 1990s ≈ 464, 2000s ≈ 480, **2010s ≈ 305**.
Reanálise de sensibilidade:

| Amostra | n | Sen (mm/déc.) [IC 95 %] | MK p |
|---|---|---|---|
| Todos | 47 | −60 [−87; −33] | < 0,001 |
| Sem 2011/12–2016/17 | 41 | −49 [−75; −21] | 0,007 |
| Até 2010/11 | 37 | −45 [−79; 2] | 0,073 |
| 1963–2008 | 35 | −45 [−81; 6] | 0,13 |

A **magnitude** é estável (−45 a −60) e a tendência não é só a seca de 2012–2017, o que é um bom
resultado. Mas a **significância** depende da última década: até 2010/11, o IC inclui zero. Ao
mesmo tempo, os anos chuvosos posteriores (2021/22, 2023/24) estão fora pelo critério FCA, e
2022/23 e 2024/25 faltam. A série termina, portanto, numa década seca e sem os anos que a
sucederam. Um modelo linear tem AIC menor que um degrau (476 contra 481), o que apoia a leitura
como tendência.

A frase do resumo "essas tendências resistem às análises de sensibilidade" precisa incluir esta
tabela, ou ao menos a linha "até 2010/11".

### M4. Uma queda de ~30 % na média exige confirmação independente

A média passa de ≈ 597 mm (antes de 1991) para ≈ 415 mm (de 1991 em diante), e o modelo implica
queda de ≈ 12 % por década. É um resultado forte para o semiárido, onde boa parte da literatura
não encontra tendência anual significativa. A confirmação atual vem de uma referência regional em
que duas das quatro vizinhas encerram o registro antes do fim do período, e a parte local tem
quebra marginal em 2003/04.

Dados que permitiriam essa confirmação são públicos: postos pluviométricos da ANA (Hidroweb) no
entorno de Petrolina e Juazeiro, a grade diária de Xavier et al. (já citada) e, possivelmente, a
estação agrometeorológica de Bebedouro (Embrapa Semiárido), em Petrolina (vale confirmar a
disponibilidade). Para o resultado central, a limitação 7 ("dados em grade permitiriam…") é
insuficiente: um revisor pediria ao menos uma dessas comparações.

---

## 4. Pontos menores

1. ✅ *(resolvido)* **Tabela 2 (probabilidades) não segue o resultado principal.** A linha "Estacionário" usa a
   Gamma (T = 81; 31 % em 30 anos), enquanto o resultado estacionário principal é a média por
   Akaike (T = 38, que daria ≈ 55 % em 30 anos). Alinhar ou explicar.
2. ✅ *(resolvido)* **Faixa "29 a 140 anos"** no resumo inclui a log-normal, rejeitada pelo Anderson–Darling
   (p = 0,012). Excluí-la ou dizer que foi rejeitada.
3. **Datas da estação chuvosa condicionais.** As tendências de início e duração excluem os 7 anos
   sem estação definida, que são anos secos e concentrados no período recente (1975/76, 1992/93,
   1997/98, 2011/12, 2012/13, 2015/16, 2016/17). A tendência descreve os anos *com* estação,
   provavelmente subestimando a mudança. Dizer isso.
4. **Contagem "11 de 16 indicadores significativos"** inclui o vento, que o próprio texto diz
   não interpretar climaticamente. Informar a contagem sem ele ou anotar.
5. **Período no título.** "1961–2025" é o período bruto; a chuva analisada vai de 1963/64 a
   2020/21 e a temperatura, de 1971 a 2023. Considerar o período efetivo no título ou deixá-lo claro
   na primeira frase do resumo.
6. **Critério FCA 5 %.** 2023/24 fica fora com FCA de 5,4 %. A sensibilidade a 7,5 % já cobre
   isso, mas convém citar o caso, porque é exatamente um ano chuvoso no fim da série.
7. **Disponibilidade de código.** O texto aponta para `github.com/dantebertuzzi/meteo`, mas o
   remoto está no commit antigo `23291c2`: o código do artigo não está lá. Fazer o push antes do
   depósito.
8. **Literatura local.** Comparar com estudos de tendência já publicados para Petrolina/Juazeiro e
   para o Submédio São Francisco, se houver, para situar a magnitude de −60 mm/década.

---

## 5. Reprodutibilidade (verificado)

- `scripts/verifica.py` (Python, independente): **34/34** valores conferem; nenhuma inconsistência
  interna nas macros.
- `test/runtests.jl`: todos os testes passam.
- `latexmk`: compila sem avisos.

## 6. Prioridade sugerida

1. M1 (muda o número central do resumo) e o ponto menor 1 (mesma família de números).
2. M2 (muda a afirmação sobre a Tmáx).
3. M3 (acrescentar a tabela de sensibilidade; não muda conclusões, muda a força delas).
4. M4 (dados externos; maior esforço, maior ganho de credibilidade).
5. Pontos menores 2–8.
