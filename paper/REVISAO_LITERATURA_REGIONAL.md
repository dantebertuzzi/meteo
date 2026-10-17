# Revisão crítica (simulada) da literatura regional — 17/10/2026

Pareceres, no formato de revisão por pares, dos estudos sobre tendências de chuva e temperatura
em Petrolina e no Submédio São Francisco encontrados na busca de 17/10/2026. O objetivo é decidir
**como** citá-los no artigo: como confirmação, como contraste ou só como contexto.

**Base:** textos completos de Pimentel e Assis (2022) (cópia do Internet Archive; o site da UFJF
estava fora do ar), Gois et al. (2024), Lacerda et al. (2010), Teixeira (2011), Assis et al. (2022)
e Nascimento et al. (2023). Rezende et al. (2023) foi lido por um resumo do texto completo (sem
acesso direto ao PDF). Onde possível, os números foram conferidos contra a série bruta da 82983.

Critérios usados em todos os pareceres:
1. **Homogeneidade e qualidade:** a série foi testada contra referências, ou só checada por
   consistência?
2. **Tendência:** o teste considera a autocorrelação? Há magnitude com intervalo de confiança?
3. **Multiplicidade:** com muitos índices ou estações, há controle de falsos positivos?
4. **Coerência interna:** texto, tabelas e conclusões dizem a mesma coisa?

---

## 1. Pimentel e Assis (2022)

*Análise da variabilidade climática no município de Petrolina-PE entre os anos de 1973–2021.*
Revista de Geografia (PPGEO-UFJF), v. 12, n. 2, pp. 281–303. DOI 10.34019/2236-837X.2022.v12.39026.

**Resumo.** O estudo usa chuva, Tmáx, Tmín e Tméd anuais da 82983 entre 1973/74 e 2021, com falhas
preenchidas por regressão linear com a estação de Bebedouro (Embrapa). Os anos são classificados em
"anos-padrão" (quantis), e os anômalos são comparados com ENSO e o Dipolo do Atlântico.
Mann–Kendall: chuva Z = −3,2; Tmáx 4,6; Tmín 5,9; Tméd 5,5.

**Pontos fortes.** É o estudo mais próximo do nosso: mesma estação, período longo e mesmo sinal nas
quatro variáveis. Também discute os moduladores (ENSO, Dipolo) e compara as normais 1961–1990 e
1991–2020.

**Problemas principais**
- **P1. Preenchimento entre instrumentos sem avaliação.** Bebedouro fica a cerca de 40 km. O ano de
  1985 inteiro e outras falhas foram preenchidos por regressão simples, sem quantificar o erro nem
  testar se a mistura cria tendência ou quebra artificiais. Para a chuva diária ou mensal no
  semiárido, a correlação entre postos a essa distância é modesta. Nós evitamos esse problema usando
  um único instrumento e deixando os anos incompletos de fora.
- **P2. Nenhum teste de homogeneidade da temperatura.** A série começa em 1973. Nós encontramos
  quebras na Tmáx frente a Remanso em 1975 e a Cabrobó em 2001, e indício de quebra frente à A307.
  Os valores Z de 4,6 a 5,9 podem estar inflados por essas quebras. Na Tmín, que no nosso caso não
  tem quebra, o sinal é robusto.
- **P3. Mann–Kendall sem correção para autocorrelação e sem magnitude.** O artigo só traz Z, sem
  declive de Sen, intervalo de confiança ou p exato. O leitor não sabe quanto a chuva caiu. A
  magnitude só aparece na reta do Excel, sem valor no texto.
- **P4. Uso das normais 1961–1990 como referência da temperatura.** O artigo conclui que as
  temperaturas subiram comparando as normais 1961–1990 e 1991–2020. A normal de Tmáx de 1961–1990
  inclui os anos anteriores a 1966, que diferem até 5 °C das vizinhas (nosso § 3.3). Parte da
  "mudança" entre normais é inomogeneidade.

**Problemas menores**
- A Tabela 2 rotula as classes como "P5–25, 25–50, 50–75" com frequências esperadas de 5/20/50/20/5 %,
  mas o texto descreve quartis, que dariam 25 % cada. Os sinais das desigualdades estão trocados
  ("Quente < 33,3" deveria ser "> 33,3").
- O período varia entre "1973–2021", "1974–2021", "48 anos" e "49 anos".
- O texto diz que as temperaturas "médias e mínimas" subiram em todos os meses e, logo depois, que
  as "médias" caíram no inverno; provavelmente o primeiro termo era "máximas".
- "Desde 2009 os totais ficaram abaixo da normal 1961–1990", mas 2009 teve 747,7 mm, acima de 577 mm.
- A associação com ENSO e o Dipolo é qualitativa, sem teste nem medida de efeito. O próprio texto
  mostra que 2012 e 2017, os anos mais secos, não tiveram El Niño.

**Veredito: revisão maior.** O sinal das tendências é crível e concorda com o nosso, mas as
magnitudes não são reportadas e a homogeneidade não é tratada.

**Uso no nosso artigo.** Citar na Introdução como estudo anterior com a mesma estação e o mesmo sinal.
Na Discussão, dizer o que acrescentamos: magnitudes com intervalo de confiança, homogeneidade
(Tmáx com quebras, Tmín sem), um único instrumento, verificação com a ANA e período de retorno.
Não é preciso apontar os defeitos no texto; basta a comparação.

---

## 2. Gois et al. (2024)

*Análise da tendência da precipitação pluviométrica em Bebedouro, Estado de Pernambuco (PE).*
III Simpósio de Ciências Ambientais na Amazônia Sul Ocidental (UFAC), 8 p.
DOI 10.5281/zenodo.13306022.

**Resumo.** Chuva mensal e anual de Bebedouro (Embrapa), 1975–2017, declarada sem falhas.
Shapiro–Wilk e Bartlett para "normalidade e homogeneidade", Mann–Kendall, Sen e Pettitt.
Resultado: Sen = −6,83 mm/ano e ruptura de Pettitt em 2010.

**Pontos fortes.** O dado é independente do nosso: outro instrumento, outro operador e outra rede, a
cerca de 20–40 km. A magnitude, **−68 mm/década**, é muito próxima da nossa (−60).

**Problemas principais**
- **P1. Confusão conceitual sobre homogeneidade.** O teste de Bartlett compara variâncias entre
  grupos; aqui, entre pares de meses do calendário. Isso não tem relação com homogeneidade
  climatológica, que é a ausência de quebras não climáticas na série. Nenhum teste de homogeneidade
  propriamente dito foi feito.
- **P2. Testes irrelevantes para a pergunta.** Que os dados mensais não sejam normais não afeta o
  Mann–Kendall, que é não paramétrico. A seção ocupa metade do texto sem consequência para as
  conclusões.
- **P3. Tendência e degrau não são distinguidos.** O Pettitt aponta ruptura em 2010, a 7 anos do fim.
  Com um degrau, o Mann–Kendall também dá "tendência", e os autores não discutem qual modelo descreve
  melhor a série. No nosso caso, o melhor degrau cai em 2009/10, mas a tendência linear vence por
  AIC (476,1 contra 481,4).
- **P4. Sem correção para autocorrelação e sem intervalo de confiança para o Sen.**

**Problemas menores**
- "Tendência de redução significativa **na maioria dos anos**": o Mann–Kendall dá um único resultado
  por série.
- "Diferença entre as médias antes e após a ruptura de 18,19 mm e 306,44 mm" não é interpretável:
  são duas diferenças? duas médias?
- Declarar ausência de falhas em 43 anos sem descrever o controle de qualidade é pouco crível.

**Veredito: rejeitar como artigo** (como resumo de simpósio, aceitável). O número central, −68
mm/década com ruptura por volta de 2010, é plausível e útil como indicação, mas não tem incerteza.

**Uso no nosso artigo.** Citar como confirmação independente da magnitude, com cautela ("resumo de
simpósio; sem intervalo de confiança"). A ruptura em 2010 conversa com o nosso resultado de que a
detecção depende da década de 2010.

---

## 3. Lacerda et al. (2010)

*Índices climáticos extremos para o município de Petrolina, PE.* XVI Congresso Brasileiro de
Meteorologia, Belém. (Repositório Alice/Embrapa.)

**Resumo.** Bebedouro (Embrapa), 1965–2007, chuva e temperatura diárias. RClimdex, com tendências
por mínimos quadrados. Chuva total: −4,2 mm/ano (p = 0,06). P10mm e P20mm caem (p ≈ 0,11). TXx
sobe 0,019 °C/ano (p = 0,049), TNx sobe 0,039 °C/ano (p = 0,001) e a amplitude térmica sobe
0,027 °C/ano (p = 0,004).

**Pontos fortes.** O procedimento segue o padrão ETCCDI, reporta declive, p e erro-padrão, e tem
dados independentes dos nossos.

**Problemas principais**
- **P1. Critério de significância frouxo e aplicado de forma inconsistente.** Os autores adotam
  p < 0,1 como "boa significância" e marcam como não significativos só os valores com p > 0,3. P10mm
  (p = 0,106) e P20mm (p = 0,110) ficam acima de 0,1 e, mesmo assim, entram como resultado. Com nove
  índices, não há controle de multiplicidade.
- **P2. Mínimos quadrados sem considerar autocorrelação** e sem teste de homogeneidade. O controle de
  qualidade do RClimdex detecta valores impossíveis, mas não quebras.
- **P3. Os rótulos "Tmáx" e "Tmin" da Tabela 2 são, pela definição da Tabela 1, TXx e TNx** (máximos
  mensais), e não médias. O texto fala em "aumento das temperaturas máximas e mínimas", o que leva o
  leitor a entender médias.

**Problema menor.** "55 % da chuva entre janeiro e abril" é coerente com os outros estudos, mas não
traz período nem método.

**Veredito: aceitável como trabalho de congresso; revisão maior para periódico.**

**Uso no nosso artigo.** É um contraste útil:
- A chuva cai −42 mm/década (p = 0,06) em Bebedouro até 2007, coerente com a nossa queda
  "não detectável até 2010/11" (−45, p = 0,21).
- A amplitude térmica **sobe** em Bebedouro, enquanto na 82983 ela não muda (−0,01 °C/década). Vale
  mencionar como diferença entre estações.

---

## 4. Teixeira (2011)

*Tendências climáticas no centro de desenvolvimento Petrolina-PE/Juazeiro-BA.* XVII Congresso
Brasileiro de Agrometeorologia, Guarapari. (Repositório Alice/Embrapa.)

**Resumo.** Estações agrometeorológicas de Bebedouro (Petrolina) e Mandacaru (Juazeiro), 1978–2008.
São analisadas tendências lineares de radiação, insolação, temperatura, UR, vento, chuva e
evaporação no tanque Classe A. Em Petrolina:
- chuva −4 mm/ano;
- Tmáx +0,04 °C/ano e Tmín −0,01 °C/ano;
- UR +0,01 %/ano e vento −0,02 m/s por ano;
- **evaporação do tanque −1 mm/ano**.

**Pontos fortes.** É o único trabalho regional que olha para a demanda evaporativa medida (tanque
Classe A, radiação, vento, UR). A observação de que o aumento da UR e a queda do vento compensam o
efeito da Tmáx sobre a evaporação é fisicamente correta.

**Problemas principais**
- **P1. Nenhum teste de significância nem intervalo de confiança.** Todas as tendências são
  coeficientes de regressão sem incerteza. Algumas são praticamente nulas (insolação +0,001 h/dia por
  ano; UR +0,01 % por ano) e mesmo assim são descritas como "elevações".
- **P2. Atribuição causal não testada.** A explicação pela irrigação (mais vapor, menos radiação,
  menos vento por rugosidade) é plausível, mas é conjectura. Falta comparar com estações fora das
  áreas irrigadas.
- **P3. Radiação e insolação divergem** (radiação cai, insolação sobe). Os autores atribuem isso ao
  vapor d'água, mas uma deriva do actinógrafo é igualmente provável e não é discutida.

**Problema menor.** A conclusão diz "maiores valores de deficiência hídrica em Juazeiro que em
Petrolina", mas o resultado mostra deficiência **aumentando** em Petrolina (+3 mm/ano) e
**diminuindo** em Juazeiro (−12 mm/ano). A frase confunde nível com tendência.

**Veredito: revisão maior.** Os resultados são descritivos e sem incerteza.

**Uso no nosso artigo.** É o mais relevante para a R4.1:
- Em Bebedouro, a **evaporação do tanque Classe A caiu** e o vento também caiu. Isso reforça que a
  alta da evaporação Piche na 82983 (+8,7 %/década) não deve ser lida como clima, e que a demanda
  evaporativa não aumentou de forma clara.
- A Tmín em queda contrasta com a nossa alta de +0,39 °C/década e com Lacerda et al. (2010), que
  acham TNx subindo na mesma estação. Sem testes, não dá para levar a divergência muito longe, mas
  vale registrar.

---

## 5. Assis et al. (2022)

*Analysis of climate indices and impacts on the rainfall regime in the Sub-medium stretch of the
São Francisco River Basin, Brazil.* Revista Principia, v. 59, n. 4. DOI 10.18265/1517-0306a2021id5570.

**Resumo.** 83 postos da APAC e da ANA foram triados, e 36 entraram na análise (26 em Pernambuco,
10 na Bahia), de 1964 a 2016. São calculados 11 índices ETCCDI de chuva com RClimdex, com tendências
por mínimos quadrados. PRCPTOT cai de 2,76 a 13,67 mm/ano. Em Petrolina, −5,41 mm/ano (cerca de
−54 mm/década). Também caem Juazeiro (−4,59), Remanso (−5,13) e Paulo Afonso (−3,51). CDD sobe;
R10mm, R20mm, SDII, R95p e R99p caem na maioria dos postos.

**Pontos fortes.** Tem cobertura espacial (36 postos), dados da ANA e da APAC independentes do INMET,
e índices padronizados. Os valores para Petrolina, Juazeiro, Remanso e Paulo Afonso concordam em
sinal e ordem de grandeza com os nossos (−60 na 82983; −42 na referência regional).

**Problemas principais**
- **P1. A conclusão contradiz a tabela.** O texto afirma que PRCPTOT é negativo "unanimemente, com
  significância estatística" em todos os postos. A Tabela 2, porém, tem valores **positivos** em
  Afogados da Ingazeira (+0,906), Arcoverde (+0,711) e Ibimirim (+1,081). O próprio texto reporta
  queda "de 2,76 a 13,67 mm/ano", o que exclui Afrânio (−1,58), Jacobina (−0,27) e Sento Sé (−0,63).
  A afirmação de unanimidade é falsa.
- **P2. Sem controle de multiplicidade nem significância de campo.** São 36 postos × 11 índices =
  396 testes, a p < 0,05 e p < 0,1, sem FDR nem teste de campo (por exemplo, Livezey–Chen). A
  correlação espacial entre postos vizinhos não é considerada.
- **P3. "Avaliação detalhada da qualidade e homogeneidade" sem método.** Não se sabe quais testes
  foram usados, nem por que 47 dos 83 postos saíram.
- **P4. Mínimos quadrados sem considerar autocorrelação**, o que é sensível à seca de 2012–2016 no
  fim da série.

**Problemas menores**
- A paginação diverge: o rodapé do PDF diz pp. 1224–1235, e o índice da revista diz pp. 1475–1486.
  Conferir antes de citar.
- Não informa qual posto é "Petrolina" (código da ANA ou da APAC); pode ser o próprio pluviômetro
  do INMET, o que tiraria a independência.
- A significância é indicada só por negrito ou itálico; na versão em texto, essa informação se perde.

**Veredito: revisão maior.** O padrão regional é crível e útil, mas a conclusão central exagera o
que a tabela mostra.

**Uso no nosso artigo.** Contexto regional: a queda é ampla no Submédio e mais intensa ao norte, em
Pernambuco. É útil para a Limitação 9 (representatividade espacial). Não citar como "unânime".

---

## 6. Rezende et al. (2023)

*Greening and Water Use Efficiency during a period of high frequency of droughts in the Brazilian
semi-arid.* Frontiers in Water, v. 5. DOI 10.3389/frwa.2023.1295286.

**Resumo.** SPI-12 da grade de Xavier no ponto da Embrapa (CPATSA), em 1961–1980, 1981–2000 e
2001–2020. Secas extremas passam de 7–9 % para 17–19 % dos meses; secas excepcionais surgem só em
2001–2020. O foco do artigo é vegetação (MODIS e GLEAM), não clima.

**Pontos fortes.** Tem revisão por pares em periódico internacional e usa um método padrão (SPI-12)
com períodos comparáveis.

**Problemas principais (para o uso que nos interessa)**
- **P1. Dado em grade, não estação.** A grade de Xavier interpola pluviômetros; no ponto de
  Petrolina, ela é fortemente influenciada pela própria 82983 e pelos postos da ANA, de modo que não
  é uma confirmação independente.
- **P2. Inconsistências internas.**
  - A frequência de secas excepcionais em 2001–2020 aparece como 7–9 % e, mais adiante, como 9–11 %.
  - Tendências são resumidas por "R²" com valores negativos (−0,33), o que é impossível para um
    coeficiente de determinação; provavelmente é um coeficiente de correlação ou um declive.
  - A Tabela 4 e a Discussão divergem sobre a eficiência no uso da água em Petrolina.
- **P3. Frequências por classes, sem teste** de diferença entre períodos.

**Veredito, só quanto à parte climática:** revisão menor para corrigir as inconsistências. A
conclusão qualitativa (mais secas extremas recentemente) é coerente com a nossa.

**Uso no nosso artigo.** Opcional, como contexto ao lado de Brito et al. (2018). Nós já mostramos que
a concentração recente dos eventos de SPI-12 reflete em parte a cobertura e a referência 1991–2020;
o mesmo cuidado vale para eles.

---

## 7. Nascimento et al. (2023)

*Analysis of annual precipitation variability and the rainy season period for the municipality of
Petrolina, Pernambuco.* Journal of Hyperspectral Remote Sensing, v. 13, n. 3, pp. 361–374.
DOI 10.29150/jhrs.v13.3.p361-374.

**Resumo.** Chuva mensal do INMET de 1991 a 2021. Os valores conferem com a 82983: 1995 = 702,0;
2004 = 691,4; 2009 = 747,7; 2012 = 141,1; 2017 = 147,4 mm. Aplica o índice de anomalia de chuva
(IAC) e regressão linear. A reta anual tem inclinação −5,08 mm/ano, com R² = 7,7 %.

**Pontos fortes.** Os dados conferem com a fonte, e o artigo descreve bem a década de 2010.

**Problemas principais**
- **P1. A tendência é afirmada sem teste, e não é significativa.** Com R² = 0,077 e n = 31, o p da
  regressão é **0,13**. O artigo conclui por "tendência de diminuição" e pela "tendência de aumento"
  em abril sem nenhum teste.
- **P2. Período curto e todo dentro da fase seca.** Em 1991–2021, a série começa depois da maior
  parte da queda (a normal de 1961–1990 era de 577 mm). Isso subestima a mudança de longo prazo.
- **P3. "Sem falhas" sem critério.** A série tem 1 a 5 dias ausentes em 7 desses 31 anos (1992,
  2000, 2001, 2011, 2019, 2020 e 2021). É pouco, mas os totais anuais são somas e o artigo não diz como
  tratou esses dias.

**Problema menor.** O IAC depende da média das dez maiores e das dez menores observações do próprio
período, o que torna as classes dependentes da amostra e não comparáveis entre estudos.

**Veredito: revisão maior.**

**Uso no nosso artigo.** Dispensável. Se citado, só como descrição de 2012 e 2017 como extremos
negativos.

---

## Síntese

| Estudo | Dado | Período | Chuva (mm/década) | Temperatura | Incerteza? | Veredito |
|---|---|---|---|---|---|---|
| **Este artigo** | 82983 (um instrumento) | 1963/64–2020/21 | **−60 [−87; −33]** | Tmín +0,39; Tmáx +0,19 a +0,30 | IC, TFPW, BH | — |
| Pimentel e Assis 2022 | 82983 + Bebedouro (preenchida) | 1973–2021 | Z = −3,2 (sem magnitude) | Z = 4,6 / 5,9 / 5,5 | Não | Revisão maior |
| Gois et al. 2024 | Bebedouro | 1975–2017 | −68 (ruptura 2010) | — | Não | Rejeitar (resumo) |
| Lacerda et al. 2010 | Bebedouro | 1965–2007 | −42 (p = 0,06) | TXx +0,19; TNx +0,39; DTR sobe | Só p | Congresso |
| Teixeira 2011 | Bebedouro | 1978–2008 | −40 | Tmáx +0,4; Tmín −0,1; tanque −10 mm/década | Não | Revisão maior |
| Assis et al. 2022 | 36 postos APAC/ANA | 1964–2016 | Petrolina −54; Juazeiro −46 | — | Só p, sem FDR | Revisão maior |
| Rezende et al. 2023 | Grade (Xavier) | 1961–2020 | SPI: secas extremas 2× | — | Não | Revisão menor |
| Nascimento et al. 2023 | 82983 | 1991–2021 | −51 (p = 0,13, n.s.) | — | Não | Revisão maior |

**Leitura de conjunto**
1. **O sinal da chuva é unânime, e a magnitude é coerente** entre −40 e −68 mm/década, em três
   fontes instrumentais diferentes (INMET, Embrapa e ANA/APAC). O nosso valor fica no meio do
   intervalo. Nenhum dos estudos reporta um intervalo de confiança comparável ao nosso.
2. **A temperatura mínima é o ponto de divergência.** Pimentel e Assis (Z = 5,9) e Lacerda (TNx)
   acham alta. Teixeira acha leve queda em Bebedouro até 2008, mas sem teste. A Tmáx sobe em todos.
3. **A demanda evaporativa não aumenta claramente**: o tanque Classe A cai em Bebedouro, e o vento
   cai e a UR sobe nas duas estações. Isso é coerente com a nossa ETo sem tendência significativa e
   reforça a leitura de que a alta da Piche é instrumental.
4. **O que nenhum deles faz:** testar homogeneidade contra referências, verificar com pluviômetros
   independentes, corrigir a autocorrelação, dar magnitudes com intervalo de confiança ou estimar
   período de retorno. É aí que o nosso artigo acrescenta.

## Recomendação para o artigo

| Onde | O que | Estudos |
|---|---|---|
| Introdução | Estudos anteriores com a estação: mesmo sinal, mas sem magnitude ou homogeneidade | Pimentel e Assis 2022 (junto de Jatobá 2017) |
| Discussão, comparação | Bebedouro (outro instrumento) dá −40 a −68 mm/década, com ruptura por volta de 2010 | Gois 2024, Lacerda 2010, Teixeira 2011 |
| Discussão, comparação | Queda ampla no Submédio, mais intensa em PE | Assis 2022 |
| Discussão, temperatura | Tmín: divergência com Teixeira 2011; DTR sobe em Bebedouro e não na 82983 | Teixeira 2011, Lacerda 2010 |
| Resultados ou Discussão, R4.1 | Tanque Classe A em queda em Bebedouro, coerente com a ETo sem tendência e com a Piche não climática | Teixeira 2011 |
| Não citar | Afirmações sem teste ou com tendência não significativa | Nascimento 2023 |

Rezende et al. (2023) é opcional, como contexto, ao lado de Brito et al. (2018).
