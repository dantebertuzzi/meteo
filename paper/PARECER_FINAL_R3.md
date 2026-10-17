# Parecer final de revisão por pares (simulada), 3ª rodada — versão de depósito — 17/10/2026

**Manuscrito:** versão final com o DOI 10.5281/zenodo.23438824. 28 p.
**Escopo:**
- releitura das partes alteradas;
- reprodução completa a partir do **pacote de depósito limpo**, sem o repositório;
- conferência do conteúdo dos arquivos do Zenodo;
- nota.

> **Aviso de independência.** Revisão feita pelo mesmo assistente que ajudou a escrever o artigo.

## 1. Reprodução a partir do pacote de depósito

Os zips de código e de dados foram extraídos numa pasta vazia e rodados sem nada do repositório
original.

| Verificação | Resultado |
|---|---|
| Testes unitários (`test/runtests.jl`) | 81/81 |
| Verificação independente (`verifica.py`) | 70/70 |
| Pipeline completo (`scripts/paper.jl`) | Roda sem erro; o `paper/resultados.tex` gerado é **idêntico** ao do artigo, exceto a linha de data |

**Conteúdo do pacote:** só o que o artigo usa. Ficaram de fora:
- o CSV de Curitiba;
- o arquivo horário não usado;
- o prompt de trabalho;
- os pareceres simulados;
- o relatório, as figuras e as tabelas da análise exploratória antiga da A307.

Nem o README nem o artigo citam os arquivos excluídos.

## 2. Problemas encontrados e corrigidos

1. **Duplo negativo** em 5 lugares: "a chuva anual **caiu −60** mm por década" no resumo, na
   Discussão e na Conclusão 1, "declined by −60" no abstract, e "ele cai −69" (Juazeiro).
   Com o sinal, a frase dizia o contrário do pretendido. Agora: "teve tendência de −60" e
   "had a trend of −60". Quatro passagens parecidas também foram ajustadas ("a queda passa a
   −43" → "a tendência passa a −43"; "caem de −57 para −37" → "passam de −57 para −37"; e
   "reduzem"/"reduz" → "atenuam"/"atenua"). Nenhuma rodada anterior tinha pegado isso.
2. **Ambiguidade** em "a análise regional dá 31 e 16 anos" (resumo, abstract e Conclusão 4): agora
   diz "31 anos sob estacionariedade e 16 com a média de 2020/21".
3. **Ordem das seções:** o apêndice vinha antes da Disponibilidade de dados e das Referências.
   Agora é Conclusões → Disponibilidade → Referências → Apêndice A.

## 3. Observações sem pedido de mudança

- O resumo tem cerca de 300 palavras, acima das ~250 usuais em periódicos. Para o preprint,
  está adequado.
- Números menores que dez aparecem em algarismos vindos das macros ("4 estações vizinhas"), o
  que é aceitável num texto técnico.

## 4. Veredito e nota

**Veredito: aceitar.** Pronto para publicação no Zenodo, sem pendências.

Nota mantida em **8,9 / 10**. As correções de redação compensam a densidade do texto, e a
reprodução a partir do pacote de depósito confirma o 10 em reprodutibilidade.
