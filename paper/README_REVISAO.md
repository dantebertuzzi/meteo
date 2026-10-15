# Revisão crítica do artigo — pendências (15/10/2026)

> Itens marcados [x] aplicados em `main.tex` em 15/10/2026 (compila, 18 p.).
> McKee 1993: páginas 179–184 confirmadas pela citação padrão; nota VERIFICAR removida.

Avaliação de poder estatístico e consistência de `paper/main.tex`, feita com os valores
reais de `paper/resultados.tex`. Linhas citadas referem-se ao `main.tex` atual.

**Veredito:** para o Zenodo como preprint, está quase pronto. Para revista com revisão por
pares, ainda não. Corrigir os itens da seção 3 (≈ 2 h) deixa o preprint defensável.

---

## 1. Onde o poder estatístico é bom (manter)

- **Tendência da chuva:** −60 mm/década [−87; −33], MK e TFPW p<0,001. Resiste ao critério
  FCA ≤ 7,5 % (−47 [−75; −19]) e aparece na referência regional (−42). Resultado mais sólido.
- **Temperatura:** Tmín +0,39 e Tmáx +0,28 °C/década; resistem ao corte em 2008; regressão
  quantílica confirma em toda a distribuição.
- **Índices de intensidade:** Rx1day, SDII, R95p com q ≤ 0,003.

## 2. Onde o poder é fraco (ajustar o tom)

- [x] **Período de retorno:** IC de 19 a 220 anos é ordem de grandeza. A conclusão 3
  (L688–690) apresenta "38 anos" como número central — acrescentar que é ordem de grandeza.
- [x] **Persistência:** 43 transições; tendência e persistência praticamente não se separam
  com 47 anos, e o u_t vem de uma tendência ajustada aos mesmos dados (absorve qualquer
  agrupamento). Trocar "a persistência desaparece" (resumo L64–65 e inglês L89–90) por
  "não há evidência de persistência além da tendência".
- [x] **ENSO (L517–522):** a interpretação está invertida. r quase não muda sem tendência
  (−0,20 → −0,19), ou seja, o sinal do ENSO **não** vem da tendência. A queda de
  P(seco | El Niño) de 0,57 para 0,29 é 8/14 contra 4/14 anos — ruído de contagem pequena.
  Reescrever: correlação na direção esperada, não significativa, insensível à remoção da
  tendência; probabilidades condicionais instáveis com 14 anos por fase.
  Ajustar também a Discussão (L646–648).

## 3. Inconsistências a corrigir

- [x] **3.1 Datas da estação chuvosa (L531–535 e Discussão L640–641).** O texto diz que não
  há mudança, mas:
  - início: Sen +8,6 d/déc. [1,9; 15,5] (IC exclui zero), MK p = 0,014, HR p = 0,051,
    TFPW p = 0,16;
  - duração: Sen −9,7 d/déc. [−19,0; −0,7] (IC exclui zero), MK p = 0,040, HR p = 0,007,
    TFPW p = 0,31.
  Escrever "evidência mista": só o TFPW não é significativo, e a TMD do início
  (12,1 d/déc.) é maior que o efeito estimado → teste sem poder, não ausência de mudança.
- [x] **3.2 "A queda vem da intensidade" (L561).** R1mm também cai: −2,5 dias/década
  [−4,6; −0,6], q(HR) = 0,012, q(TFPW) = 0,059. A frequência de dias de chuva também diminui.
  Corrigir também a conclusão 5 (L693–694).
- [x] **3.3 Datas da seca, hoje com três versões:** "2011–2012" (intro L104), "2011–2016"
  (Discussão L625), "2011–2017" (L338, L685, L691, L501). Escolher uma; se a de Brito et al.
  for 2011–2016, dizer explicitamente que em Petrolina a sequência se estende até 2016/17.
- [x] **3.4 Período efetivo vs. título/resumo.** Título e resumo dizem 1961–2025, mas a
  análise da chuva vai até 2020/21 (índices até 2020). Os anos chuvosos 2021/22 e 2023/24
  ficam fora pelo critério FCA ≤ 5 % e atenuam a tendência. Dizer isso no resumo (PT e EN):
  é um viés possível na direção do resultado principal.
- [x] **3.5 Empate 2011/12 × 2012/13.** 151,4 mm contra 151,3 mm. Dizer que estão
  praticamente empatados (L58, L339, L685).
- [x] **3.6 Cauda da Gamma.** O texto (L374–375) admite cauda mais fina que os dados
  (empírico T ≈ 24; Gamma T = 81), mas o modelo não estacionário principal é Gamma.
  O T = 31 provavelmente está **superestimado**. Acrescentar uma frase em Resultados e uma
  em Limitações.

## 4. Risco científico maior (não resolvível com estes dados — enfatizar)

Cerca de um quarto da tendência local (−60 contra −42 da referência) é excesso local, com
quebra marginal em 2003/04 (SNHT p = 0,048; Pettitt p = 0,049) coincidente com a instalação
da A307. Um revisor pode argumentar que parte do resultado central é instrumental.

- [x] Dar mais destaque ao cenário com inclinação regional (T = 47 anos em 2025): citá-lo no
  resumo e na conclusão 3, não só nas tabelas.

## 5. Observações menores

- [x] SPI-12: todos os 6 eventos são de 1993 em diante (L355). Explicar se isso vem das
  lacunas (SPI-12 exige 12 meses seguidos) ou do período de referência 1991–2020 mais seco.
- [x] E3: T ≈ 85 anos (empírico 41) com um único evento — já ressalvado; manter assim.
- [x] Limitação 4 (L663) fala em quebra da Tmáx "frente a uma vizinha", mas os Resultados
  só mostram isso para antes de 1971. Alinhar.

## 6. Pendências objetivas

- [x] Páginas de McKee 1993 marcadas `VERIFICAR` em `refs.bib`.
- [ ] Normais oficiais INMET 1991–2020 (portal estava fora do ar) — ou manter como limitação 8.
- [ ] DOI do Zenodo (hoje "a ser atribuído", L703).
- [ ] Commitar o projeto (nada está commitado).

## Como regenerar

Números só via macros:

```bash
julia --project=. -t 4 scripts/paper.jl            # recalcula tudo e escreve paper/resultados.tex
julia --project=. -t 4 scripts/paper.jl figuras    # só figuras, usando o cache
cd paper && latexmk -pdf main.tex
```
