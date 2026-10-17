# Parecer de revisão por pares (simulada) — 4ª rodada — 17/10/2026

**Manuscrito:** *Menos chuva e mais calor em Petrolina-PE (1961–2025)…*, versão de 23 p.
(commit 01bb0a4), depois das normais oficiais do INMET e do resumo reduzido.
**Base:** leitura integral do PDF compilado, macros de `paper/resultados.tex`, `verifica.py`
(57/57) e reanálises rápidas em Python sobre os dados brutos da 82983 (código no fim).

**Recomendação:** **aceitar com revisão menor.** O núcleo científico (P2 forte, P1 em ordem de
grandeza, incerteza tratada como resultado) está sólido e bem defendido depois de quatro rodadas.
Os pontos abaixo são de redação e de uma afirmação sem suporte nos dados. Nenhum muda números
do resumo ou das conclusões.

---

## 1. Reavaliação das rodadas anteriores

| Ponto | Situação |
|---|---|
| M1 (período de retorno não estacionário dependia da Gamma) | Resolvido: média de Akaike das 8 candidatas, T = 12 [5,0; 88] |
| M2 (quebras na Tmáx) | Resolvido: 0,19–0,30 °C/déc. conforme a hipótese; Tmín é o resultado mais firme |
| M3 (detecção depende da década de 2010) | Resolvido e declarado no resumo |
| M4 (homogeneidade da chuva) | Resolvido: 5 postos da ANA; Juazeiro (3 km) reproduz −69 mm/déc. |
| R2.1–R2.3 (forma variável, validação fora da amostra) | Resolvidos |
| R3.1 (sensibilidade ao fim da série) | Resolvido: 11 anos com FCA ≤ 7,5 % |
| Linha de base externa (normais INMET) | Resolvido em 17/10: 419 mm nas duas; 577 → 483 → 419 mm |
| Resumo longo; mês mais chuvoso | Resolvidos em 17/10 (~250 palavras; jan./mar. empatam) |

Os números do texto continuam vindo todos de macros, e a verificação independente confere 57/57.

## 2. Pontos de conteúdo

### R4.1 — "Mais evaporação" não é sustentado pela ETo do próprio artigo

A Discussão diz: *"Com noites frias ficando mais raras, a demanda evaporativa aumenta; por isso o
SPEI […] é uma medida mais completa"* e *"Menos chuva e mais evaporação apontam na mesma direção"*.
O artigo não reporta tendência de ETo. Recalculando nos dados (1971+, anos com ≥ 10 meses válidos,
MK-TFPW e Sen):

| Série (82983) | n | Sen por década | p (TFPW) |
|---|---|---|---|
| ETo Hargreaves (a do SPEI) | 42 | +0,4 % | 0,28 |
| Evaporação Piche | 35 | +11 % | < 0,001 |
| Insolação | 43 | +3,2 % | < 0,001 |
| Vento | 35 | −8,9 % | < 0,001 |
| UR | 47 | +0,6 % | 0,31 |

(No pipeline Julia, com médias anuais por anomalias mensais, a ETo dá +0,56 %/déc. [−0,01; 1,09],
p = 0,069, e a Piche +8,7 %/déc. [4,8; 12,2], p < 0,001; a conclusão é a mesma.)
A ETo usada no SPEI **não tem tendência detectável**: Hargreaves depende de T̄ e de √DTR, e a DTR
não muda. A Piche sobe muito (+11 %/déc., cerca de +50 % em cinco décadas), mas o evaporímetro de
Piche é muito sensível à exposição e à manutenção. Além disso, o vento cai, e a UR não muda. Assim,
"mais evaporação" é, no máximo, uma hipótese física, e não um resultado. O raciocínio "noites frias
mais raras → mais demanda evaporativa" também é fraco: a demanda depende sobretudo de Tmáx,
radiação, déficit de vapor e vento.
*Pedido:* reescrever as duas frases sem afirmar aumento de evaporação. Por exemplo: dizer que o
aquecimento tende a aumentar a demanda evaporativa, mas que a ETo de Hargreaves não mostra
tendência detectável e que a Piche, que sobe, não é confiável para tendência. Convém reportar a
tendência da ETo (uma linha na Tabela 6 ou no texto, via macro) e acrescentar a Piche à lista de
"não interpretamos climaticamente", junto com o vento.

### R4.2 — O resumo trata as normais oficiais como uma confirmação a mais

No resumo, *"A queda aparece na referência regional, em 4 dos pluviômetros da ANA e nas normais
oficiais da estação"* põe as normais na mesma lista das verificações independentes. Elas são
calculadas com os mesmos dados da 82983, como a própria Limitação 10 reconhece. A frase sugere três
confirmações quando há duas.
*Pedido:* separar. Por exemplo: "…e em 4 dos 5 pluviômetros da ANA; as normais oficiais da
estação registram a mesma queda (577 mm em 1961–1990, 419 mm em 1991–2020)". O mesmo vale para o
abstract.

### R4.3 — Figura 9 × normal oficial de 1961–1990

A Fig. 9a mostra a média de 1961–1990 calculada com os meses válidos da série (março ≈ 149 mm,
janeiro ≈ 70 mm). O texto da Área de estudo cita a normal oficial (março = 136 mm). Um leitor que
compare os dois verá números diferentes para o "mesmo" período. A diferença é esperada: a série tem
20 a 22 meses válidos por mês do calendário em 1961–1990, e a normal do INMET tem outros critérios.
*Pedido:* uma frase na legenda ("média dos meses válidos da série; difere da normal oficial do
INMET, que usa outros critérios") ou, melhor, os pontos da normal oficial de 1961–1990 no painel.

### R4.4 — Dias chuvosos no resumo sem a ressalva do teste principal

O resumo afirma *"com menos dias chuvosos e chuva menos intensa"*. Na Tabela 6, o R1mm tem q = 0,059
no teste principal (MK-TFPW com BH), embora o IC de Sen exclua zero e o Hamed–Rao dê q = 0,012. Os
índices de intensidade são claros, mas a frequência é marginal no teste que o próprio artigo elegeu
como principal.
*Pedido:* "chuva menos intensa e, com evidência mais fraca, menos dias chuvosos", ou equivalente. Na
Discussão (§ "O regime mudou") e na Conclusão 1, a mesma ressalva em poucas palavras.

## 3. Pontos editoriais

1. **Referência da FAO-56:** sai "Food **e** Agriculture Organization". O biblatex trata `and` no
   campo `institution` como separador de lista. Basta envolver o nome em chaves:
   `institution = {{Food and Agriculture Organization of the United Nations}}`.
2. **Lista repetida:** os sete anos sem estação chuvosa (`\AnosSemEstacao`) aparecem por extenso em
   dois parágrafos seguidos (§ 5.5). Na segunda vez, basta "esses anos".
3. **Resumo:** "em 4 dos pluviômetros da ANA" → "em 4 dos 5" (já sai com R4.2).
4. **Disponibilidade de dados:** citar também as normais do INMET e o
   `scripts/extrai_normais_inmet.py`, ao lado do `baixa_ana.jl`.
5. **T empírico de 24 anos (§ 5.2):** a contagem 2/47 inclui o próprio 2012/13, que foi escolhido
   por ser extremo, e mistura anos de uma série com tendência. A comparação "cauda da Gamma mais fina
   que a dos dados" continua válida pela validação fora da amostra. Vale meia frase de ressalva ou
   uma remissão à Tabela 3.

## 4. Observações sem pedido de mudança

- A afirmação de que a 82983 "reproduz" a normal oficial é, a rigor, uma checagem de processamento
  (mesmos dados do BDMEP). A Limitação 10 já diz isso, e está adequado.
- A Fig. 6 avisa que a curva não estacionária chega a zero por causa da normal. Isso está honesto,
  mas para periódico um revisor pode pedir para cortar a curva abaixo de ~50 mm.
- A normal oficial de Tmáx de 1961–1990 (32,5 °C) acima da média de 1971–1990 (31,8 °C) reforça a
  exclusão da temperatura antes de 1971. Esse uso das normais é bom.

## 5. Veredito

**Aceitar com revisão menor.** Para o depósito no Zenodo, recomendo aplicar R4.1, porque é a única
afirmação do texto sem suporte nos dados, e também R4.2 e o item editorial 1. R4.3, R4.4 e os demais
itens editoriais são rápidos e melhoram o texto, mas não bloqueiam. Depois disso, o preprint está
pronto. Para submissão a periódico, a principal fragilidade remanescente é estrutural e já está
declarada: uma só estação para a P1.

---

### Código da checagem de R4.1

```python
# .venv-verif/bin/python -I  (a partir da raiz do repositório)
import pandas as pd, pymannkendall as mk
d = pd.read_csv('data/raw/bdmep/dados_82983_D_1961-01-01_2025-06-30.csv', sep=';', skiprows=10,
                na_values='null', decimal=',').iloc[:, :10]
d.columns = ['data','evap','inso','prec','tmax','tmed','tmin','ur','urmin','vento']
d['data'] = pd.to_datetime(d.data); d = d[d.data.dt.year >= 1971]
d['ano'], d['mes'] = d.data.dt.year, d.data.dt.month
def anual(v):                      # mês com ≤ 10 dias ausentes; ano com ≥ 10 meses
    g = d.groupby(['ano','mes'])[v]
    m = pd.DataFrame({'ok': g.count(), 'n': g.size(), 'val': g.mean()}); m = m[m.ok >= m.n - 10]
    a = m.reset_index().groupby('ano').agg(n=('mes','size'), val=('val','mean'))
    return a[a.n >= 10].val
e = pd.read_csv('results/paper/indices_mensais.csv'); e = e[e.ano >= 1971]
eto = e.groupby('ano').eto.agg(['count','mean']); eto = eto[eto['count'] >= 10]['mean']
for nome, s in [('Piche', anual('evap')), ('Insolação', anual('inso')), ('Vento', anual('vento')),
                ('UR', anual('ur')), ('ETo Hargreaves', eto)]:
    r = mk.trend_free_pre_whitening_modification_test(s.values); sl = mk.sens_slope(s.values)
    print(nome, len(s), f"{100*10*sl.slope/s.mean():.1f} %/déc.", f"p={r.p:.3f}")
```

---

## Resolução — 17/10/2026

- **R4.1:** tendências de ETo (Hargreaves) e da evaporação Piche calculadas em `scripts/paper.jl`
  (`SenEToPct`, `SenPichePct`), fora da família BH, para não mudar os q-valores da Tabela 6.
  ETo: +0,56 %/déc. [−0,01; 1,09], p = 0,069; Piche: +8,7 %/déc., tratada como o vento (sem
  interpretação climática). A Discussão deixou de afirmar "mais evaporação": a mudança no balanço
  hídrico vem sobretudo da queda da chuva, e o SPEI pouco acrescenta ao SPI nesta série.
- **R4.2:** no resumo e no abstract, as normais oficiais foram separadas das verificações
  independentes ("4 dos 5 pluviômetros"; as normais "registram a mesma queda").
- **R4.3:** a legenda da Fig. 9 explica que as médias usam os meses válidos e que a de 1961–1990
  difere da normal oficial.
- **R4.4:** "com evidência mais fraca, menos dias chuvosos" no resumo, no abstract, em Resultados,
  na Discussão e na Conclusão 1.
- **Editoriais:** FAO com chaves duplas; lista dos anos sem estação não se repete; disponibilidade
  de dados cita as normais e o script; T empírico com ressalva e remissão à Tabela 3.

Recompilado (23 p., sem referências indefinidas); `verifica.py` 57/57; as demais macros não mudaram.
