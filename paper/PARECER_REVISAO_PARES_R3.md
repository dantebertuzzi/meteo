# Parecer de revisão por pares (simulada) — 3ª rodada — 15/10/2026

**Manuscrito:** versão de 23 p., depois de M1–M4 e R2.1–R2.3.
**Recomendação:** **revisão menor**, já aplicada nesta rodada. Do ponto de vista científico, o
preprint está pronto para depósito; falta o push do código (ação do autor) e o DOI.

## Ponto de conteúdo

**R3.1 — sensibilidade do estimador principal ao fim da série.** O T não estacionário principal
(média de Akaike, 12 anos) era avaliado em 2020/21, no fim de uma década seca e sem os anos
chuvosos seguintes; a sensibilidade ao critério FCA ≤ 7,5 % existia só para a Gamma.
*Resolvido:* com FCA ≤ 7,5 %, a média de Akaike dá 11 anos [4,6; 69] em 2023/24 e 13 em 2020/21
(`TNSMediaFCASeteUlt`, Tabela 4 e Resultados). O resultado principal não depende da exclusão
desses anos.

## Pontos editoriais (todos resolvidos)

1. "duração diminui −9,7" → "varia −9,7" (duplo negativo).
2. Legenda da Tabela 4 dizia "quatro candidatas" também no caso não estacionário (são oito).
3. Resumo PT/EN: "quatro distribuições" → "cada uma com forma fixa ou variável no tempo".
4. Resumo PT/EN: "sem mudança relativa" → "sem mudança relativa significativa" (Q 1,15 → 1,34).
5. "mês mais chuvoso … é jan." → nome do mês por extenso (macro).
6. Figura 7: eixo y com ponto decimal → vírgula.
7. Limitação 3: "a queda é em parte local" (afirmado como fato) → "é regional, mas pode ter uma
   parte local", coerente com os pluviômetros da ANA.

## Observações sem mudança

- O resumo tem cerca de 380 palavras; para periódico, convém reduzir a ~250.
- A climatologia de Bebedouro (Teixeira 2010) indica março como mês mais chuvoso em 1963–2009,
  e a 82983 indica janeiro em 1991–2020; períodos diferentes, sem contradição, mas um revisor pode
  perguntar.

## Atualização — 17/10/2026

- **Normais oficiais do INMET** (portal de volta; exige user-agent de navegador): a climatologia
  1991–2020 calculada reproduz a normal oficial (419 mm; temperaturas mensais com diferença ≤ 0,1 °C).
  As normais anteriores registram a mesma queda (577 mm em 1961–1990, 483 mm em 1981–2010). Novo
  parágrafo "Normais oficiais" em Dados; a limitação "Linha de base externa" foi reescrita.
  Dados em `data/raw/inmet_normais/` via `scripts/extrai_normais_inmet.py`.
- **Mês mais chuvoso:** na normal oficial 1991–2020, janeiro (83,5 mm) e março (83,6 mm) empatam;
  em 1961–1990 março se destacava (136 mm), como em Bebedouro. O texto deixa de afirmar "janeiro".
- **Resumo** reduzido de ~380 para ~245 palavras (PT e EN), agora citando as normais oficiais.
