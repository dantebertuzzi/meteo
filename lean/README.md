# MeteoLean — verificação formal das fórmulas fechadas do artigo

Provas em Lean 4 + Mathlib das identidades algébricas que o código Julia usa para transformar
probabilidades em períodos de retorno, riscos e contagens. Os enunciados são simples: o Lean **não**
verifica os números do artigo, que dependem dos dados, dos otimizadores e do bootstrap, nem as
escolhas de modelo. Ele verifica que essas fórmulas estão certas. O grupo de testes `Identidades
provadas em Lean` de `test/runtests.jl` confere que as funções Julia usadas no artigo
(`prob_ao_menos_um`, `prob_n_variavel`, `prob_ao_menos_k`, `markov_posterior`, `ns_media_akaike`)
satisfazem cada teorema. O tempo médio de espera é conferido por simulação, porque nenhuma função do
artigo o calcula diretamente.

| Arquivo | Teorema | O que garante | Código Julia |
|---|---|---|---|
| `Retorno.lean` | `soma_primeiro_evento` | $\sum_{k<n}(1-p)^k p = 1-(1-p)^n$ (Eq. 4) | `prob_ao_menos_um` |
| | `tempo_medio_espera` | tempo médio até o 1º evento $=1/p$ (justifica $T=1/p$) | — |
| | `probVar_const`, `probVar_mono` | $1-\prod(1-p_j)$ reduz-se à Eq. 4 e é monótona em cada $p_j$ | `prob_n_variavel` |
| | `continua_ge_congelado` | se $p_j$ só cresce, "tendência continua" ≥ "clima congelado" (Tab. 6) | `bloco_nao_estacionario` |
| `Akaike.lean` | `pesos_pos`, `pesos_soma_um` | pesos de Akaike positivos, somam 1 | `ns_media_akaike` |
| | `media_entre` | $\bar p$ fica entre o menor e o maior $p_i$ | |
| | `T_media_le_media_T` | $1/\sum w_i p_i \le \sum w_i/p_i$: o $T$ da média de $p$ nunca excede a média dos $T$ | |
| `Markov.lean` | `estacionaria`, `piSeco_mem` | $\pi_S=a/(a+b)$ é estacionária e está em $[0,1]$ | `markov_posterior` |
| | `duracao_media_seca` | duração média de uma sequência seca $=1/b$ | |
| | `sem_persistencia`, `persistencia_iff` | sem persistência, $\pi_S=P(S\mid N)$; persistência $\iff a+b<1$ | |
| `PoissonBinomial.lean` | `pb_nonneg`, `pb_suporte`, `pb_soma` | a recursão dá uma distribuição | `prob_ao_menos_k` |
| | `pb_media` | média $=\sum p_i$ ("esperado" da Tab. 4) | `bloco_validacao` |
| | `pb_ao_menos_um` | $P(N\ge1)$ da recursão $=1-\prod(1-p_i)$ | |
| | `pb_binomial` | com $p$ fixo, a recursão é a binomial | |

Nenhuma prova usa `sorry`. `#print axioms` mostra só `propext`, `Classical.choice` e `Quot.sound`.

## Uso

```bash
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
cd lean
lake exe cache get     # baixa o Mathlib pré-compilado (alguns GB)
lake build             # verifica todas as provas
```
