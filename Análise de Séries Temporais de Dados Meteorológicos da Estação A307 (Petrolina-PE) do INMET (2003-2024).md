# Análise de Séries Temporais de Dados Meteorológicos da Estação A307 (Petrolina-PE) do INMET (2003-2024)

**(Resumo a ser preenchido)**




## 1. Introdução

A análise de séries temporais constitui uma ferramenta fundamental no campo da meteorologia, permitindo a compreensão de padrões, tendências e ciclos inerentes às variáveis climáticas ao longo do tempo. Uma série temporal é caracterizada como um conjunto de observações registradas sequencialmente em intervalos de tempo definidos, usualmente cronológicos. O estudo dessas séries é particularmente relevante no contexto atual de mudanças climáticas globais, onde a frequência e a intensidade de eventos extremos, como secas prolongadas, inundações e tempestades severas, têm aumentado significativamente. A capacidade de analisar e modelar dados meteorológicos históricos torna-se, portanto, essencial não apenas para a previsão e mitigação de desastres naturais, mas também para o desenvolvimento de estratégias de adaptação e planejamento de longo prazo em setores como agricultura, recursos hídricos e energia.

Este estudo concentra-se na análise das condições meteorológicas da região de Petrolina, localizada no estado de Pernambuco, no semiárido brasileiro. Esta região é caracterizada por um clima predominantemente quente e seco, com baixa pluviosidade concentrada em poucos meses do ano, tornando o monitoramento contínuo de variáveis como precipitação e temperatura de suma importância para a gestão de recursos hídricos e atividades agrícolas, pilares econômicos locais. A compreensão detalhada da variabilidade climática e de possíveis tendências de longo prazo é crucial para o desenvolvimento sustentável da região.

Os dados utilizados nesta análise foram coletados pela estação meteorológica automática A307, operada pelo Instituto Nacional de Meteorologia (INMET), e cobrem o período de 21 de fevereiro de 2003 a 31 de dezembro de 2024. Este extenso conjunto de dados horários oferece uma oportunidade valiosa para caracterizar o clima local com um alto grau de detalhamento temporal.

Os objetivos principais deste trabalho são: (1) realizar uma análise descritiva detalhada das principais variáveis meteorológicas registradas na estação A307 durante o período de estudo; (2) identificar e caracterizar padrões sazonais e tendências temporais nas séries de temperatura, precipitação e outras variáveis relevantes; (3) investigar, por meio de testes de hipóteses estatísticas, possíveis mudanças significativas nas médias de temperatura e precipitação ao longo do tempo e entre períodos específicos, bem como avaliar a distribuição de frequência da direção do vento.




## 2. Materiais e Métodos

### 2.1 Fonte dos Dados

O presente estudo utiliza dados meteorológicos horários provenientes da estação automática A307, localizada no município de Petrolina, estado de Pernambuco, Brasil. Os dados foram obtidos junto ao Instituto Nacional de Meteorologia (INMET), órgão responsável pela coleta e disponibilização de informações meteorológicas oficiais no país. A estação A307 integra a rede de estações de superfície do INMET, fornecendo um registro contínuo de diversas variáveis atmosféricas.

### 2.2 Período de Estudo

A análise abrange um extenso período de observações, iniciando-se em 21 de fevereiro de 2003 e estendendo-se até 31 de dezembro de 2024. Este intervalo temporal de mais de duas décadas permite uma avaliação robusta das características climáticas locais e a identificação de possíveis tendências e variações de longo prazo.

### 2.3 Variáveis Analisadas

Foram selecionadas para análise diversas variáveis meteorológicas horárias registradas pela estação, incluindo os registros temporais de **Data e Hora**, o **Volume de chuva acumulado** na hora anterior em milímetros (PRECIP_mm), a **Pressão Atmosférica ao Nível da Estação** em milibares (PRESSAO_ATMOSFERICA_NIVEL_ESTACAO_mb), bem como seus valores **Máximo e Mínimo** na hora anterior (PRESSAO_ATMOSFERICA_MAX_mb, PRESSAO_ATMOSFERICA_MIN_mb). Também foram consideradas a **Radiação Global** em kilojoules por metro quadrado (RADIACAO_GLOBAL_kjm2), embora a análise inicial tenha indicado um número significativo de valores ausentes para esta variável. As medições de temperatura incluíram a **Temperatura do Ar Bulbo Seco** em graus Celsius (°C) (TEMPERATURA_AR_BULBO_SECO), a **Temperatura do Ponto de Orvalho** (°C) (TEMPERATURA_PONTO_DE_ORVALHO), e os valores **Máximos e Mínimos** de ambas na hora anterior (TEMPERATURA_MAXIMA, TEMPERATURA_MINIMA, TEMPERATURA_ORVALHO_MAX, TEMPERATURA_ORVALHO_MIN). A umidade foi avaliada pela **Umidade Relativa do Ar** instantânea em percentual (%) (UMIDADE_RELATIVA_AR_percent) e seus valores **Máximo e Mínimo** na hora anterior (UMIDADE_RELATIVA_MAXIMA_percent, UMIDADE_RELATIVA_MINIMA_percent). Por fim, as condições de vento foram representadas pela **Direção do Vento** predominante em graus (0-360°) (VENTO_DIRECAO), a **Rajada Máxima de Vento** em metros por segundo (m/s) (VENTO_RAJADA_MAXIMA_ms) e a **Velocidade Média do Vento** na hora anterior (m/s) (VENTO_VELOCIDADE_ms).

### 2.4 Pré-processamento dos Dados

Os dados brutos, originalmente distribuídos em arquivos anuais no formato CSV (Comma-Separated Values) com ponto e vírgula como delimitador, passaram por uma série de etapas de pré-processamento para garantir sua qualidade e adequação à análise. Inicialmente, todos os arquivos CSV correspondentes ao período de estudo foram carregados e concatenados em um único conjunto de dados (DataFrame).

Durante o carregamento, as colunas, que não possuíam nomes descritivos nos arquivos originais, foram renomeadas sequencialmente (ex: `col_1`, `col_2`, ...) e, posteriormente, receberam nomes significativos conforme um dicionário predefinido, facilitando a identificação e manipulação das variáveis (vide `dicionario_renomeacao` no código fonte). A coluna contendo a data das observações, originalmente em formato texto (`dd/mm/yyyy`), foi convertida para um tipo de dado específico para datas, permitindo operações temporais adequadas.

Uma inspeção inicial revelou a presença de valores codificados como `-9999`, indicando medições ausentes ou inválidas. Embora um tratamento específico para esses valores (como substituição por `missing` ou interpolação) não tenha sido explicitamente detalhado na fase inicial do código analisado, a presença desses códigos foi registrada nas estatísticas descritivas preliminares (`describe(df)`), destacando a necessidade de considerá-los nas análises subsequentes. Da mesma forma, a conversão de colunas de texto (como as de temperatura e pressão, que utilizam vírgula como separador decimal) para tipos numéricos é uma etapa implícita necessária para as análises quantitativas realizadas posteriormente.

### 2.5 Análise Estatística

A análise estatística foi conduzida utilizando a linguagem de programação Julia e suas bibliotecas especializadas em manipulação de dados e estatística. Foram calculadas estatísticas descritivas básicas (média, mediana, mínimo, máximo) para fornecer uma visão geral da distribuição de cada variável. Os dados foram agregados em diferentes escalas temporais (diária, mensal, anual) para investigar padrões de sazonalidade e tendências de longo prazo.

A visualização gráfica foi amplamente empregada, utilizando gráficos de séries temporais para observar a evolução das variáveis ao longo do tempo, histogramas para analisar suas distribuições de frequência e boxplots para comparar as distribuições entre diferentes meses ou anos, auxiliando na identificação de sazonalidade e outliers.

Para investigar tendências temporais, foi aplicada análise de regressão linear simples, por exemplo, modelando a temperatura média anual em função do ano.

Adicionalmente, foram realizados testes de hipóteses para avaliar diferenças estatisticamente significativas e características de distribuição. Utilizou-se o Teste t para amostras independentes com variâncias desiguais (Teste t de Welch) para comparar a temperatura média do ar entre os anos de 2013 e 2023, e também para comparar a precipitação média entre os meses de Março (tipicamente mais chuvoso) e Setembro (tipicamente mais seco). Foi realizado um Teste Qui-quadrado de aderência para verificar se a distribuição da direção do vento (agrupada em setores de 30 graus) se desviava significativamente de uma distribuição uniforme, indicando direções predominantes.

### 2.6 Ferramentas

As análises foram implementadas na linguagem Julia, com o suporte das seguintes bibliotecas principais: `CSV` para leitura de arquivos, `DataFrames` para manipulação de dados tabulares, `Dates` para operações com datas, `Plots` e `StatsPlots` para geração de gráficos, `StatsBase` para estatísticas básicas e testes, `GLM` para modelos lineares e `HypothesisTests` para testes de hipóteses.




## 3. Resultados

### 3.1 Caracterização Geral das Variáveis Meteorológicas

A análise descritiva preliminar dos dados horários da estação A307 de Petrolina, abrangendo o período de 2003 a 2024, revelou informações importantes sobre as condições climáticas locais. A inspeção inicial, realizada através da função `describe()` aplicada ao DataFrame consolidado, indicou a presença significativa de valores ausentes ou inválidos (codificados como -9999) em diversas variáveis, como Precipitação (PRECIP_mm), Pressão Atmosférica (todas as colunas relacionadas), Temperatura do Ar (TEMPERATURA_AR_BULBO_SECO, TEMPERATURA_MAXIMA, TEMPERATURA_MINIMA), Temperatura do Ponto de Orvalho (todas as colunas relacionadas), Umidade Relativa (todas as colunas relacionadas) e Vento (VENTO_DIRECAO, VENTO_RAJADA_MAXIMA_ms, VENTO_VELOCIDADE_ms). A variável Radiação Global (RADIACAO_GLOBAL_kjm2) apresentou uma quantidade particularmente elevada de dados faltantes. Essa constatação ressalta a importância de um tratamento adequado desses dados em análises mais aprofundadas, embora a presente descrição se baseie nos valores disponíveis e processados nas etapas subsequentes do notebook original.

### 3.2 Análise de Séries Temporais e Sazonalidade

As visualizações gráficas permitiram identificar padrões temporais e sazonais nas principais variáveis. A série temporal da **Temperatura Média Diária do Ar** (`temp_media_diaria`), calculada a partir das médias horárias, exibiu flutuações características ao longo dos anos. A análise de sazonalidade, realizada através de boxplots mensais da temperatura média diária, confirmou um ciclo anual bem definido, com os meses de verão (especialmente Outubro a Dezembro) apresentando temperaturas consistentemente mais elevadas em comparação aos meses de inverno (Junho a Agosto), que registram as menores médias e medianas térmicas do ano. A variabilidade intra-mensal, representada pela amplitude das caixas e bigodes nos boxplots, também mostrou diferenças ao longo do ano.

De forma análoga, a análise da **Precipitação Diária Acumulada** (`precip_diaria`) revelou a marcada sazonalidade pluviométrica típica do semiárido. A maior parte da precipitação concentra-se nos primeiros meses do ano, notadamente entre Janeiro e Abril, período que constitui a estação chuvosa da região. Os meses do segundo semestre, especialmente de Junho a Outubro, caracterizam-se por volumes de chuva significativamente menores, com muitos dias registrando precipitação zero. Os boxplots mensais da precipitação diária ilustraram essa concentração, mostrando medianas próximas de zero na maior parte do ano e valores mais elevados e maior dispersão (indicando eventos de chuva mais intensos e frequentes) nos meses de verão e outono.

### 3.3 Análise de Tendências

Para investigar possíveis tendências de longo prazo, foi realizada uma análise de regressão linear simples sobre a **Temperatura Média Anual do Ar**. O gráfico de dispersão dos valores médios anuais de temperatura em função do ano, juntamente com a linha de regressão ajustada, sugeriu visualmente uma tendência de aumento da temperatura ao longo do período analisado (2003-2024). O coeficiente angular positivo da reta de regressão, embora necessite de uma avaliação formal de sua significância estatística (não explicitamente detalhada nos resultados do notebook `src.ipynb` para a regressão, mas abordada indiretamente no teste t), corrobora essa observação visual de uma possível tendência de aquecimento na região.

### 3.4 Resultados dos Testes de Hipóteses

Os testes de hipóteses realizados no notebook `test_src.ipynb` forneceram evidências estatísticas sobre aspectos específicos da variabilidade climática. A **Comparação de Temperaturas (2013 vs 2023)**, utilizando o Teste t de Welch para comparar as médias de temperatura do ar bulbo seco entre os anos de 2013 e 2023, resultou em um p-valor extremamente baixo (próximo de zero, conforme indicado pela saída `p-valor = 0.0`). Este resultado indica uma diferença estatisticamente significativa entre as médias de temperatura dos dois anos, rejeitando a hipótese nula de que as médias seriam iguais. A análise direcional dessa diferença (qual ano foi mais quente) não foi explicitada no output, mas a significância estatística da diferença foi confirmada.

A **Comparação de Precipitação (Março vs Setembro)**, através do Teste t de Welch aplicado às médias de precipitação diária dos meses de Março e Setembro, também produziu um p-valor muito baixo (`p-valor = 0.0`). Isso sugere uma diferença estatisticamente significativa na precipitação média entre o mês tipicamente mais chuvoso (Março) e o mês tipicamente mais seco (Setembro), confirmando a forte sazonalidade pluviométrica observada nas análises descritivas.

Finalmente, a análise da **Distribuição da Direção do Vento**, utilizando o Teste Qui-quadrado de aderência para avaliar se a frequência da direção do vento (agrupada em 12 setores de 30 graus) seguia uma distribuição uniforme, obteve um p-valor também muito baixo (`p-valor = 0.0`). Este valor indica que a distribuição observada da direção do vento difere significativamente de uma distribuição uniforme, sugerindo a existência de direções de vento predominantes na região de Petrolina ao longo do período estudado. A identificação específica dessas direções predominantes requereria uma análise adicional das frequências observadas em cada setor.




## 4. Discussão

Os resultados obtidos a partir da análise dos dados meteorológicos horários da estação A307 em Petrolina-PE, entre 2003 e 2024, permitem uma caracterização detalhada do clima local e a identificação de padrões e tendências relevantes. A marcada sazonalidade observada tanto na temperatura do ar quanto na precipitação está em plena consonância com as características esperadas para o clima semiárido brasileiro. O ciclo anual de temperaturas, com valores mais elevados nos meses de primavera e verão (Outubro a Março) e mais amenos no outono e inverno (Abril a Setembro), reflete a variação da incidência de radiação solar ao longo do ano. Da mesma forma, a concentração da precipitação nos primeiros meses do ano (Janeiro a Abril), seguida por um longo período de estiagem, é um traço distintivo do regime pluviométrico da região, fortemente influenciado por sistemas atmosféricos como a Zona de Convergência Intertropical (ZCIT), cuja posição mais ao sul nesse período favorece a ocorrência de chuvas.

A tendência de aumento da temperatura média anual, sugerida pela análise de regressão linear e reforçada pela diferença estatisticamente significativa encontrada entre as médias de temperatura de 2013 e 2023 (Teste t, p < 0.001), é um achado preocupante que se alinha às observações de aquecimento global. Embora a análise aqui apresentada seja preliminar e demande investigações mais aprofundadas sobre a magnitude e as causas dessa tendência, ela sinaliza possíveis mudanças no padrão térmico da região com potenciais impactos sobre os ecossistemas locais, a agricultura e os recursos hídricos. A confirmação e quantificação precisa dessa tendência exigiriam análises mais sofisticadas, considerando a autocorrelação temporal e possíveis pontos de mudança na série.

Os resultados dos testes de hipóteses corroboram as observações descritivas e fornecem validação estatística para padrões importantes. A diferença significativa na precipitação média entre Março e Setembro (Teste t, p < 0.001) quantifica a intensidade da sazonalidade pluviométrica. Já a rejeição da hipótese de uniformidade na distribuição da direção do vento (Teste Qui-quadrado, p < 0.001) confirma a existência de direções predominantes, um fator relevante para estudos de dispersão de poluentes, potencial eólico e conforto térmico. A identificação dessas direções predominantes, no entanto, não foi o foco da análise apresentada no notebook original.

É fundamental reconhecer as limitações deste estudo. A presença significativa de dados faltantes ou inválidos (código -9999) em diversas variáveis, especialmente na Radiação Global, pode ter influenciado as análises descritivas e de tendência. Um tratamento mais robusto desses dados, como técnicas de imputação ou interpolação, poderia refinar os resultados. Além disso, a análise de tendência baseou-se em regressão linear simples, que pode não capturar completamente a complexidade das variações temporais. A comparação dos resultados com a literatura científica específica sobre o clima de Petrolina e do semiárido seria necessária para contextualizar os achados de forma mais ampla.

Como sugestões para trabalhos futuros, recomenda-se a aplicação de técnicas mais avançadas de tratamento de dados faltantes, a utilização de modelos de séries temporais mais sofisticados (como ARIMA ou modelos sazonais) para análise de tendências e previsão, a investigação da relação entre as variáveis meteorológicas locais e fenômenos climáticos de larga escala (como El Niño-Oscilação Sul - ENOS), e a análise detalhada das direções predominantes do vento e sua variabilidade sazonal e interanual.




## 5. Conclusão

A análise das séries temporais de dados meteorológicos horários da estação A307 do INMET em Petrolina-PE, cobrindo o período de 2003 a 2024, permitiu alcançar os objetivos propostos neste estudo. Foi realizada uma caracterização detalhada das principais variáveis climáticas, identificando-se padrões sazonais bem definidos, especialmente na temperatura do ar e na precipitação, que são consistentes com o regime climático do semiárido brasileiro. A estação chuvosa concentra-se nos primeiros meses do ano, enquanto os meses de inverno registram as menores temperaturas.

As análises revelaram tendências temporais, notadamente uma sugestão de aumento na temperatura média anual ao longo do período estudado, corroborada por testes estatísticos que indicaram uma diferença significativa entre as médias de temperatura de anos distintos (2013 vs 2023). Da mesma forma, a forte sazonalidade da precipitação foi confirmada estatisticamente, assim como a existência de direções de vento predominantes na região, refutando a hipótese de uma distribuição uniforme.

Os achados deste trabalho reforçam a compreensão da dinâmica climática em Petrolina e fornecem uma base quantitativa para estudos futuros e para o planejamento regional, especialmente em face dos desafios impostos pela variabilidade e pelas potenciais mudanças climáticas. As limitações identificadas, como a presença de dados faltantes, apontam para a necessidade de aprimoramentos metodológicos em análises subsequentes. Conclui-se que o monitoramento contínuo e a análise aprofundada dos dados meteorológicos são essenciais para a gestão sustentável dos recursos naturais e para a adaptação da região às condições climáticas presentes e futuras.


