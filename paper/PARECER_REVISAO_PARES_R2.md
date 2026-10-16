# Parecer de revisão por pares (simulada) — 2ª rodada — 15/10/2026

**Manuscrito revisado:** *Menos chuva e mais calor em Petrolina-PE (1961–2025)…* (19 p.)
**Parecer anterior:** `PARECER_REVISAO_PARES.md` (M1–M4 e pontos menores 1–8).

**Recomendação:** **revisão maior** (para periódico). O M1 foi bem resolvido e, com o ponto R2.1
abaixo, o resultado de período de retorno fica mais forte do que antes. Mas M2 e M3 seguem
abertos, e o M2 afeta uma afirmação do resumo. Como preprint, publicável depois de R2.1 e M2.

Reanálises desta rodada: Python/SciPy, independentes do código Julia, sobre
`results/paper/anos_hidrologicos_82983.csv`.

---

## 1. Resposta aos pontos da 1ª rodada

| Ponto | Situação | Comentário |
|---|---|---|
| M1 — período de retorno depende da distribuição | **Resolvido** | Média por pesos de Akaike das 4 candidatas, com bootstrap; Gamma como sensibilidade. Implementação testada (8 testes unitários) e conferida em Python (39/39). Ver R2.1, que refina a conclusão. |
| M2 — quebras da Tmáx × vizinhas | **Resolvido (15/10)** | Testes reportados em Resultados; Tmáx sob 3 hipóteses: salto subtraído 0,19 [0,10; 0,29], degrau agnóstico 0,23 [0,06; 0,41], desde 1975 0,30 [0,19; 0,38]. Resumo, Discussão, Limitação e Conclusão citam a faixa; Tmín (sem quebra) destacada como resultado mais firme. `verifica.py` 47/47. |
| M3 — significância depende da década de 2010 | **Resolvido (15/10)** | Nova tabela de sensibilidade (todos, sem 2011/12–2016/17, até 2010/11, desde 1971/72, FCA ≤ 7,5 %) e comparação linear × degrau (AIC 476,1 × 481,4). Resumo PT/EN e Discussão dizem que a magnitude é estável, mas a detecção depende da década de 2010 (até 2010/11: p = 0,21). `verifica.py` 50/50. |
| M4 — confirmação com dados independentes | **Resolvido (15/10)** | 11 pluviômetros da ANA (Hidroweb) a até ~70 km, mesmo critério FCA; 5 com ≥ 30 anos comuns. Todos caem nos anos comuns, 4/5 com IC excluindo zero (−69 a −18 mm/déc.). Juazeiro (3 km): −69 contra −69 da 82983, razão sem tendência e sem mudança em 2003/04. Junco e Lagoa Grande sugerem excesso local. `verifica.py` 57/57. |
| Menores 1 (Tabela 2) e 2 (faixa com log-normal) | **Resolvidos** | — |
| Menores 3–8 | **3–6 resolvidos (15/10); 7–8 abertos** | Resolvidos: datas condicionais da estação, contagem com vento, período efetivo no resumo, FCA de 2023/24. Literatura local resolvida em 15/10 (Costa 2020, Carvalho 2020, Jatobá 2017, Teixeira 2010, conferidos no Crossref/PDF; parágrafo de comparação na Discussão). Aberto: push do código (ação do autor). |

---

## 2. Novos pontos principais

### R2.1. A dependência da distribuição vem da forma constante, não da família

> **Resolvido em 15/10/2026.** `ns_candidata(...; forma_var = true)` acrescenta log s = s₀ + c·t;
> a média principal usa 8 candidatas (4 famílias × forma fixa/variável): **12 anos [5,0; 88]**;
> com forma variável as famílias dão 10–12 anos; LR rejeita CV constante na Gamma (p = 0,034).
> Resumo, Métodos, Resultados, Tabela 1, Discussão, Limitação e Conclusão atualizados. R2.3
> (frase da checagem preditiva) também ajustado. Testes: 11/11; `verifica.py`: 44/44.

O texto revisado diz que, no caso não estacionário, o período de retorno "vai de 10 a 36 anos
conforme a distribuição" e que "a escolha da distribuição pesa tanto quanto a amostra". A
reanálise mostra que isso decorre de uma hipótese comum às quatro candidatas: **a forma
constante no tempo**. A Gamma, a Weibull e a log-normal mantêm o coeficiente de variação fixo; a
normal mantém o desvio-padrão fixo. Os dados respondem qual das duas hipóteses é razoável:

| Período | n | Média (mm) | DP (mm) | CV |
|---|---|---|---|---|
| 1963–1990 | 17 | 597 | 152 | 0,25 |
| 1991–2005 | 15 | 461 | 166 | 0,36 |
| 2006–2021 | 15 | 370 | 159 | 0,43 |

O desvio-padrão é estável e o CV quase dobra. Deixando também a forma variar no tempo
(log da forma linear no tempo, um parâmetro a mais):

| Família | AIC forma fixa | AIC forma variável | LR p | T fixa (anos) | T variável (anos) |
|---|---|---|---|---|---|
| Normal | 609,8 | 611,8 | 0,93 | 10,1 | 10,4 |
| Weibull | 612,2 | 611,5 | 0,10 | 25,8 | 11,5 |
| Gamma | 614,5 | 612,0 | **0,034** | 35,9 | 11,5 |
| Log-normal | 618,1 | 614,9 | **0,023** | 35,6 | 12,3 |

Com a forma livre, **as quatro famílias convergem para T ≈ 10–12 anos** em 2020/21, em linha com
a média por Akaike (12). A Gamma rejeita o CV constante (p = 0,034), e a normal não precisa do
parâmetro extra (p = 0,93). Isso é uma boa notícia para o artigo: o resultado principal não é
artefato da escolha da família.

**Pedido:**
1. Incluir as versões com forma variável no conjunto de candidatas (8 modelos) ou, no mínimo,
   como sensibilidade explícita na Tabela 1.
2. Reescrever "a escolha da distribuição pesa tanto quanto a amostra" (resumo PT/EN e
   Resultados) para algo como: "o resultado depende de a variabilidade acompanhar ou não a queda
   da média; os dados indicam que não, e com essa hipótese relaxada as quatro famílias dão
   10–12 anos".
3. A limitação "Forma da distribuição" e a Discussão ("manter o desvio-padrão ou o CV constante
   leva a 10 ou 36 anos") devem refletir que a hipótese de CV constante é rejeitada pelos dados.
4. A interpretação física merece uma frase: a chuva caiu sem que a amplitude da oscilação
   interanual em milímetros diminuísse, o que aumenta a frequência relativa de anos muito secos.

### R2.2. Validação fora da amostra: a Gamma com CV constante falha

> **Resolvido em 15/10/2026** (`bloco_validacao`, Tabela de validação). Com o limiar P10 do ajuste
> (sem viés de seleção): observados 7; estacionários esperam ≈ 2 (p = 0,002), não estacionários 5,6–5,9
> (p ≈ 0,3). Com E1: Gamma CV fixo p = 0,002; média de Akaike p = 0,12.

Ajustando os modelos só até 2005/06 e prevendo 2006/07–2020/21 (15 anos válidos, com a tendência
extrapolada), o número esperado de anos com chuva ≤ E1 é:

| Modelo (ajustado até 2005/06) | Esperado em 15 anos |
|---|---|
| Gamma estacionária | 0,02 |
| Normal estacionária | 0,17 |
| Gamma não estacionária (CV fixo) | 0,07 |
| Weibull não estacionária | 0,23 |
| Normal não estacionária (DP fixo) | 0,86 |
| Média por Akaike | 0,48 |
| **Observado** | **2** (2012/13 e 2016/17) |

Sob a Gamma não estacionária, observar 2 ou mais teria probabilidade da ordem de 0,2 %; sob a
normal não estacionária, cerca de 20 %. Mesmo a média por Akaike subestima.

**Ressalva importante:** o teste não é limpo, porque o limiar E1 foi escolhido justamente por
ser um ano do período de teste (viés de seleção a favor de "ver" o evento). O teste, portanto,
não valida valores absolutos, mas a **comparação entre modelos** continua informativa: todos
foram avaliados contra o mesmo evento selecionado. Isso reforça R2.1 e torna implausíveis os
valores da Gamma com CV constante (estacionária ou não).

**Pedido:** reportar este exercício com a ressalva do viés de seleção, de preferência também com
limiares não escolhidos a partir do período de teste (p.ex. o percentil 10 do período de ajuste).

### R2.3. A checagem preditiva posterior é um teste fraco

> **Resolvido em 15/10/2026** (frase reescrita em Resultados).

O texto diz que os modelos bayesianos (todos Gamma) "reproduzem a cauda seca observada", com
base na probabilidade preditiva do mínimo amostral (0,36 e 0,27). R2.2 mostra que a mesma Gamma
subestima muito a frequência de anos secos fora da amostra. A checagem do mínimo dentro da
amostra não detecta o problema de CV constante. Suavizar a frase ou, melhor, trocar a checagem por
uma que olhe a variância ao longo do tempo (p.ex. DP dos resíduos por período).

---

## 3. Verificações desta rodada que não exigem mudança

- **Persistência com estados definidos pela normal.** Os estados "sem tendência" ainda vêm da
  Gamma, que agora é só sensibilidade. Refazendo com a normal: P(seco|seco) = 0,24, P(seco|não
  seco) = 0,33, probabilidade de persistência 0,25 (Gamma: 0,19; 0,28; 0,26). A conclusão não
  muda. Convém dizer isso numa frase, já que o modelo principal mudou.
- **Tabelas e texto** estão consistentes com as macros novas; o PDF compila sem avisos.
- **Reprodutibilidade:** `verifica.py` 39/39; testes unitários passam.

## 4. Novos pontos menores

1. ✅ *(resolvido: Fig. 4 e 6 refeitas; Fig. 7 é a posterior Gamma, legenda explica)* **Figuras 4, 6 e 7 mostram só a Gamma**, que agora é sensibilidade e cuja hipótese de CV
   constante é rejeitada (R2.1). A Figura 7 mostra posterior centrada em 23–31 anos, enquanto o
   resultado principal é 12. Atualizar a Figura 4 (linha da média entre candidatas ou da normal) e
   a Figura 6 (curva de retorno não estacionária), ou ao menos chamar a atenção na legenda.
2. ✅ *(resolvido)* **"92 % em 30 anos"** (clima de 2020/21 mantido) é um número que chama atenção. Lembrar ao
   lado dele que supõe anos independentes e o clima congelado no nível de uma década excepcionalmente
   seca; o IC [28; 99] deve aparecer sempre junto.
3. ✅ *(resolvido: r com a normal no texto; persistência também)* **Teleconexões:** a correlação sem tendência também usa a PIT da Gamma; checar com a normal,
   como em §3 (não refeito nesta rodada).
4. **Frase "Os modelos bayesianos… servem como sensibilidade à forma de quantificar a incerteza,
   não à escolha da distribuição":** boa, manter, mas combinar com R2.3.

## 5. Prioridade sugerida

1. **R2.1** (muda a narrativa do resultado de período de retorno para melhor; esforço pequeno:
   a função `ns_candidata` já tem a estrutura, falta a forma log-linear).
2. **M2** (quebra da Tmáx em 2000; afeta o resumo).
3. **R2.2** e **R2.3** (validação fora da amostra; esforço pequeno).
4. **M3** (tabela de sensibilidade da tendência; já calculada na 1ª rodada).
5. Pontos menores e **M4**.
