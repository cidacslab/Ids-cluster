# Titulo: Clusterizando curvas de covid19 por UFs do Brasil
# Autor: Jonatas Silva do Espirito Santo
# Criacao: 2021-07-31
# Ajustado: 2024-01-31

library(dplyr)
library(ggplot2)
library(dtw)
library(dtwclust)
library(purrr)
library(lubridate)
library(geobr)
library(furrr)
library(data.table)
library(readr)

# Configurações iniciais --------------------------------------------------

# definindo variavel cujas series serao clusterizadas
target <- 'totalCases_per_100k_inhabitants'
#target <- 'totaldeaths_per_100k_inhabitants'
#target <- 'newCases_per_100k_inhabitants'
#target <- 'newDeaths_per_100k_inhabitants'

# definindo o numero de cluster
n_clusters = as.numeric(6)

# Definindo tamanho da janela para cálculo do DBA
length_window <- 10

# fixando a semente
set.seed(16032018)

# Definindo a estrategia de processamento como sendo multisessao
future::plan(strategy = 'multisession', workers = 3)

# ETL ---------------------------------------------------------------------

# Pegando shapefile pór UF para construir mapa
shape_UF <- geobr::read_state()

# pegando serie historica do github de Wesley Cota
ts_UF_raw <- data.table::fread(
  input =  "https://raw.githubusercontent.com/wcota/covid19br/master/cases-brazil-states.csv"
)

# # Guardando a serie historica em disco
# data.table::fwrite(
#   x = ts_UF_raw,
#   file = "./data/ts_UF.csv",
#   sep = ";"
# )

# #Lendo arquivo do disco
# ts_UF_raw <- data.table::fread(
#   file = './data/ts_UF.csv'
# )

# Retirando a totalização do Brasil
# filtrando somente os 18 meses de pandemia no Brasil 2020-02-25 a 2021-08-25
ts_UF <- ts_UF_raw %>% 
  dplyr::filter(state!='TOTAL') %>% 
  dplyr::mutate(date = lubridate::ymd(date)) %>% 
  dplyr::mutate(
    epiweek = ifelse(
      test = epi_week < 100,
      yes = epi_week,
      no = epi_week - 47
      
    )) %>% 
  dplyr::relocate(epiweek)  %>% 
  dplyr::filter(date < "2022-01-01")

# Calculando dados da população
ts_UF <- ts_UF %>% 
  mutate(populacao = totalCases/totalCases_per_100k_inhabitants*100000)

# Agrupando dados por semana epidemiologica
ts_UF <- ts_UF %>% 
  dplyr::group_by(epiweek,state) %>% 
  summarise(
    newCases = sum(newCases),
    newDeaths = sum(newDeaths),
    populacao = dplyr::last(populacao),
    deaths = dplyr::last(deaths),
    deaths_by_totalCases = dplyr::last(deaths_by_totalCases),
    totalCases_per_100k_inhabitants = dplyr::last(totalCases_per_100k_inhabitants),
    totaldeaths_per_100k_inhabitants = dplyr::last(deaths_per_100k_inhabitants),
    deaths_by_totalCases = dplyr::last(deaths_by_totalCases)
  )

# Calculado newCases_per_100k_inhabitants e newDeaths_per_100k_inhabitants
ts_UF <- ts_UF %>% 
  mutate(
    newCases_per_100k_inhabitants = round((newCases/populacao*100000),5),
    newDeaths_per_100k_inhabitants = round((newDeaths/populacao*100000),5)
  )
  

# lista auxiliar para calcular quantis
q_aux <- list(
  c(0.33,0.66),
  c(0.25,0.5,0.75),
  c(0.2,0.4,0.6,0.8),
  c(0.1,0.3,0.5,0.7,0.9),
  c(0.05,0.25,0.45,0.65,0.85,1.0),
  c(0.05,0.20,0.35,0.60,0.75,0.9,1.0),
  c(0.05,0.15,0.25,0.40,0.60,0.75,0.9,1.0),
  c(0.05,0.15,0.25,0.35,0.50,0.65,0.75,0.85,1.0),
  c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,1.0)
) %>% setNames(2:10)

#Criando um vetor com as siglas das UF
vetor_UF <- ts_UF %>% 
  as.data.frame() %>% 
  dplyr::select(state) %>% 
  unique() %>% 
  dplyr::pull() %>% 
  sort() 

# transformando df_9 em lista para calcular DBA
list_ts_UF <- furrr::future_map(
  .x = vetor_UF,
  .f = ~{ts_UF %>% 
      dplyr::filter(state==.x) %>% 
      dplyr::pull(target) %>% 
      as.numeric()
  }  
) %>% setNames(vetor_UF)


# definindo quais UF serao os centroides iniciais
obs_cent_ini <- ts_UF %>%  
  dplyr::group_by(state) %>% 
  #dplyr::summarise(l_obs = dplyr::last(get(target))) %>% 
  dplyr::summarise(l_obs = max(get(target))) %>% 
  dplyr::arrange(l_obs) %>% 
  #dplyr::slice(order_by=round(q_aux[[as.character(n_clusters)]]*length(vetor_UF))) %>% 
  dplyr::slice(round(q_aux[[as.character(n_clusters)]]*length(vetor_UF))) %>% 
  dplyr::pull(state)

# gerando a lista com as series de centroides iniciais
centroide <- furrr::future_map(
  .x = obs_cent_ini,
  .f = ~{list_ts_UF[[.x]]
  }  
) 

n_interacoes <- 100

i=1

vec_WGSS <- numeric(0) 

list_index_cluster <- NULL

tictoc::tic()  
while (i <= n_interacoes) {
  
# Calculando dtw entre as series temporais e os centroides
index_cluster <- purrr::map_dbl(
  .x = vetor_UF,
  .f = ~{
    temp_UF <- .x
    temp_dtw <- furrr::future_map_dbl(
      .x = 1:n_clusters,
      .f = ~{
        dtw::dtw(
          list_ts_UF[[temp_UF]],
          centroide[[.x]]
        )$normalizedDistance
      }  
    )
    cluster <- match(min(temp_dtw),temp_dtw)
  }  
) %>% setNames(vetor_UF) 

# Calculando o centroide
centroide <- purrr::map(
  .x = 1:n_clusters,
  .f = ~{
    
    # nomes das UF no cluster
    names_instancias <- names(index_cluster[index_cluster==.x])
    
    # Definindo a serie referência para o dba
    state_ref <- furrr::future_map_dbl(
      .x = names_instancias,
      .f = ~{length(list_ts_UF[[.x]])}  
    ) %>% 
      purrr::set_names(names_instancias) %>% 
      base::as.data.frame() %>%
      tibble::rownames_to_column() %>%
      purrr::set_names(c('state','length')) %>% 
      dplyr::filter(length == max(length)) %>% 
      dplyr::slice(1L) %>% 
      dplyr::pull(state)
    
    # Calculando o centroido do grupo
    dtwclust::dba(
      X = list_ts_UF[names(index_cluster[index_cluster==.x])],
      #centroid = NULL,
      window.size = as.numeric(length_window),
      centroid = list_ts_UF[[state_ref]]
    )
  }  
)

# Guardando os cluster gerados em cada iteração
list_index_cluster[[i]] <- index_cluster

# Calculando distâncias de cada ponto até o centroide
# d2cc - distancia de cada serie ao centroide do seu cluster
d2cc <- purrr::map(
  .x = 1:n_clusters,
  .f = ~{
    temp_cluster <- .x
    furrr::future_map_dbl(
      .x = names(index_cluster[index_cluster==.x]),
      .f = ~{
        dtw::dtw(
          list_ts_UF[[.x]],
          centroide[[temp_cluster]]
        )$normalizedDistance^2
      }  
    )
  }  
)

# WGSS (within-group sum of square distance)
WGSS = purrr::reduce(d2cc,sum)

# Guardando o WGSS no vetor
vec_WGSS[i] <- WGSS

# Parando se houver convergência do algoritmo
if(i>1){

  # Critério de parada
  var_WGSS <- (sqrt((vec_WGSS[i]-vec_WGSS[i-1])^2))/vec_WGSS[i-1]
  
  # decidindo sobre a parada da interação
  if(var_WGSS<0.00000000000000001){break}
}

# Atualizando indice da interação
i = i + 1

}

tictoc::toc()

# Espiando o vetor de WGSS
vec_WGSS

#Espiando a lista de Index
print(list_index_cluster)

# gerando mapa
mapa_cluster <-
  dplyr::left_join(
    x = shape_UF,
    y = data.frame(state = names(index_cluster), cluster = as.character(index_cluster)),
    by = c("abbrev_state" = "state")
  )

ggplot() +
  geom_sf(
    data = mapa_cluster,
    aes(fill = cluster)
  ) +
  geom_sf_text(
    data = mapa_cluster, 
    aes(label = abbrev_state),
    color = 'gray10',
    size = 5
  )+
  labs(x='',y='')+
  theme_bw()+
  theme(
    legend.position = c(0.1, 0.2),
    legend.title = element_text(colour="Black", size=24, face="bold"),
    legend.text = element_text (size = 18),
    legend.background = element_rect(size=0.5, linetype="blank")
  )+
  scale_fill_manual(values = c("#fef0d9","#fdd49e","#fdbb84","#fc8d59","#e34a33","#b30000", "#8A0000", "#420000"))


#Gerando os gráficos dos centroides
purrr::map(
  .x = 1:n_clusters,
  .f = ~{

    ts_UF_temp <-  ts_UF %>%
      dplyr::select(c('epiweek','state',target)) %>%
      dplyr::filter(state %in% names(index_cluster[index_cluster==.x]))


    ggplot() +
      geom_line(
        data = (data.frame(
          (1:length(centroide[[.x]]))+8,
          centroide[[.x]]) %>%
          set_names(c('epiweek',target))),
        aes_string(x='epiweek',y=target)#,size = 1
      )+
      geom_vline(xintercept=c(26, 65), linetype='dashed', color=c('gray', 'black'))+
      ylim(c(0,max(ts_UF %>% pull(c(target)))))+
      xlim(c(0,80))+
      ggplot2::labs(
        y= (
          switch(
            target,
            'totalCases_per_100k_inhabitants' = 'Total Cases / 100k habitants',
            'totaldeaths_per_100k_inhabitants' = 'Total Deaths / 100k habitants',
            'newCases_per_100k_inhabitants' = 'New Cases / 100k habitants',
            'newDeaths_per_100k_inhabitants' = 'New Deaths / 100k habitants',
          )
        ),
        x = "Epiweek"
      )+
      theme_bw()+
      theme(
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_text(colour="Black", size=12),
        legend.text = element_text (size = 12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        legend.background = element_rect(fill="white", size=0.3, linetype="blank")
      )
  }
) %>% gridExtra::grid.arrange(grobs=., ncol = 3)


# Gerando o gráfico das series temporias clusterizadas
purrr::map(
  .x = 1:n_clusters,
  .f = ~{ts_UF %>% 
      dplyr::select(c('epiweek','state',target)) %>% 
      dplyr::filter(state %in% names(index_cluster[index_cluster==.x])) %>% 
      ggplot(aes_string(x='epiweek',y=target,color='state'))+
      geom_line(aes(color=state), linewidth = 1)+
      #geom_vline(xintercept=c(26, 65), linetype='dashed', color=c('gray', 'black'))+
      ylim(c(0,max(ts_UF %>% pull(c(target)))))+
      ggplot2::labs(
        y= (
          switch(
            target,
            'totalCases_per_100k_inhabitants' = 'Total Cases / 100k habitants',
            'totaldeaths_per_100k_inhabitants' = 'Total Deaths / 100k habitants',
            'newCases_per_100k_inhabitants' = 'New Cases / 100k habitants',
            'newDeaths_per_100k_inhabitants' = 'New Deaths / 100k habitants',
          )
        ), 
        x = "Epiweek"
      )+
      theme_bw()+
      theme(
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_text(colour="Black", size=12),
        legend.text = element_text (size = 12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        legend.background = element_rect(fill="white", size=0.3, linetype="blank")
      )
  }  
) %>% gridExtra::grid.arrange(grobs=., ncol = 2,)

# # Gerando o gráfico das series temporias clusterizadas com centroides
purrr::map(
  .x = 1:n_clusters,
  .f = ~{

    ts_UF_temp <-  ts_UF %>%
      dplyr::select(c('epiweek','state',target)) %>%
      dplyr::filter(state %in% names(index_cluster[index_cluster==.x]))

      ggplot(data = ts_UF_temp)+
      geom_line(aes_string(x='epiweek',y=target,color='state'))+
      geom_line(
        data = (data.frame(
          (1:length(centroide[[.x]]))+(min(ts_UF_temp$epiweek)-1),
          centroide[[.x]]) %>%
            set_names(c('epiweek',target))),
        aes_string(x='epiweek',y=target)
      )+
      geom_vline(xintercept=c(26, 65), linetype='dashed', color=c('gray', 'black'))+
      ylim(c(0,max(ts_UF %>% pull(c(target)))))+
      theme_bw()+
      theme(
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_text(colour="Black", size=8),
        legend.text = element_text (size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        legend.background = element_rect(fill="white", size=0.3, linetype="blank")
      )
  }
) %>% gridExtra::grid.arrange(grobs=., ncol = 2)

# Calculando distâncias de cada ponto até o centroide
# d2cc - distancia de cada serie ao centroide do seu cluster

d2cc <- purrr::map(
  .x = 1:n_clusters,
  .f = ~{
    temp_cluster <- .x
    purrr::map_dbl(
      .x = names(index_cluster[index_cluster==.x]),
      .f = ~{
        dtw::dtw(
          list_ts_UF[[.x]],
          centroide[[temp_cluster]]
        )$normalizedDistance^2
      }  
    )
  }  
)

# WGSS (within-group sum of square distance)
WGSS = purrr::reduce(d2cc,sum)

# Definindo a serie referência para o centroide geral
state_ref <- purrr::map_dbl(
  .x = names(list_ts_UF),
  .f = ~{length(list_ts_UF[[.x]])}  
) %>% 
  purrr::set_names(names(list_ts_UF)) %>% 
  base::as.data.frame() %>%
  tibble::rownames_to_column() %>%
  purrr::set_names(c('state','length')) %>% 
  dplyr::filter(length == max(length)) %>% 
  dplyr::slice(1L) %>% 
  dplyr::pull(state)

# Gerando o centroid geral
centroide_global <- dtwclust::dba(
  X = list_ts_UF,
  window.size = as.numeric(length_window),
  centroid = list_ts_UF[[state_ref]]
)

# d2cg - distancia dos centroides de cada cluster ao centroide global
d2cg <- purrr::map_dbl(
  .x = 1:n_clusters,
  .f = ~{
    dtw::dtw(
      centroide[[.x]],
      centroide_global
    )$normalizedDistance^2
  }  
)

# número de series temporais no cluster k
n_k <- base::table(index_cluster) %>% as.numeric()

# BGSS (between group sum of square distances)
BGSS <- purrr::map2(
  .x = d2cg,
  .y = n_k,
  .f = ~{
    .x * .y
  }  
) %>% purrr::reduce(sum)  


# Indice de Calinski-Harabasz  
CH <- (BGSS/(n_clusters-1))/(WGSS/(sum(n_k)-n_clusters))


#Visualizando as metricas
BGSS 
WGSS
CH  
