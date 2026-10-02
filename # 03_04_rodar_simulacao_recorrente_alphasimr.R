# ============================================================
# SCRIPT 04 v2
# Rodar simulação de seleção recorrente com AlphaSimR
#
# Versão correta para resultado final do artigo
#
# Correções principais:

# - n_cycles = 6
# - MT e Culling usam índice multicaracteres com escala fixa
#   baseada no Grupo B inicial
# - Genótipos dos elites são salvos no RDS final
# - Arquivo final salvo como:
#   resultados_simulacao_recorrente_alphasimr_v2.rds
#
# Importante:
# Maior score = melhor.
# Os efeitos já estão orientados para redução.
# Não inverter sinal.
# ============================================================



library(AlphaSimR)
library(dplyr)

set.seed(123)

# ============================================================
# 1. Parâmetros gerais
# ============================================================

n_rep <- 10
n_cycles <- 6

cenarios <- c(
  "AP",
  "AE",
  "FF",
  "FM",
  "MT",
  "Culling"
)

n_cruzamentos_por_ciclo <- 50
n_DH_por_cruzamento <- 20

n_TC1 <- 800
n_TC2 <- 200
n_TC3 <- 100
n_elite <- 2

prop_sel_uc <- 0.05

i_uc <- dnorm(qnorm(1 - prop_sel_uc)) / prop_sel_uc

max_candidatos_uc <- 15000

h2_traits <- c(
  AP = 0.90,
  AE = 0.91,
  FF = 0.89,
  FM = 0.89
)

fator_h2_estagio <- c(
  TC1 = 0.60,
  TC2 = 0.80,
  TC3 = 0.95
)

cat("\nParâmetros principais:\n")
cat("Repetições:", n_rep, "\n")
cat("Ciclos:", n_cycles, "\n")
cat("Cenários:", paste(cenarios, collapse = ", "), "\n")
cat("Cruzamentos por ciclo:", n_cruzamentos_por_ciclo, "\n")
cat("DHs por cruzamento:", n_DH_por_cruzamento, "\n")
cat("DHs por ciclo:", n_cruzamentos_por_ciclo * n_DH_por_cruzamento, "\n")
cat("Elites por ciclo:", n_elite, "\n")
cat("Intensidade UC:", i_uc, "\n")

# ============================================================
# 2. Carregar entradas
# ============================================================

entrada <- readRDS("sim_inputs_real.rds")

matriz_marcadores <- entrada$matriz_marcadores
mapa_AlphaSimR <- entrada$mapa_AlphaSimR
matriz_efeitos <- entrada$matriz_efeitos

entradas_cruzamentos <- readRDS(
  "entradas_cruzamentos_iniciais_simulacao.rds"
)

cruzamentos_iniciais_por_cenario <- entradas_cruzamentos$cruzamentos_iniciais_por_cenario
testadores <- entradas_cruzamentos$testadores
grupos_resumo <- entradas_cruzamentos$grupos_resumo

testador_ids <- testadores$Linhagem

parentais_grupo_B <- grupos_resumo$Linhagem[
  grupos_resumo$Grupo == "B"
]

cat("\nDimensão da matriz de marcadores:\n")
print(dim(matriz_marcadores))

cat("\nDimensão da matriz de efeitos:\n")
print(dim(matriz_efeitos))

cat("\nNúmero de parentais no Grupo B inicial:\n")
print(length(parentais_grupo_B))

cat("\nTestadores do Grupo A:\n")
print(testador_ids)

# ============================================================
# 3. Conferências iniciais
# ============================================================

if (!all(colnames(matriz_marcadores) == mapa_AlphaSimR$marker)) {
  stop("Matriz de marcadores e mapa AlphaSimR não estão alinhados.")
}

if (!all(colnames(matriz_marcadores) == rownames(matriz_efeitos))) {
  stop("Matriz de marcadores e matriz de efeitos não estão alinhadas.")
}

if (!all(parentais_grupo_B %in% rownames(matriz_marcadores))) {
  stop("Alguns parentais do Grupo B não estão na matriz de marcadores.")
}

if (!all(testador_ids %in% rownames(matriz_marcadores))) {
  stop("Alguns testadores não estão na matriz de marcadores.")
}

valores_marcadores <- sort(unique(as.vector(matriz_marcadores)))

cat("\nValores encontrados na matriz de marcadores:\n")
print(valores_marcadores)

if (!all(valores_marcadores %in% c(0, 2))) {
  stop("A matriz de marcadores deve conter apenas 0 e 2.")
}

# ============================================================
# 4. Criar escala fixa para índice multicaracteres
# ============================================================

cat("\nCriando escala fixa para o índice MT...\n")

traits_indice <- c("AP", "AE", "FF", "FM")

gv_grupo_B_inicial <- as.matrix(
  matriz_marcadores[parentais_grupo_B, , drop = FALSE] %*% matriz_efeitos
)

colnames(gv_grupo_B_inicial) <- colnames(matriz_efeitos)

media_base_indice <- colMeans(
  gv_grupo_B_inicial[, traits_indice, drop = FALSE],
  na.rm = TRUE
)

sd_base_indice <- apply(
  gv_grupo_B_inicial[, traits_indice, drop = FALSE],
  2,
  sd,
  na.rm = TRUE
)

sd_base_indice[
  sd_base_indice == 0 | !is.finite(sd_base_indice)
] <- 1

cat("\nMédias fixas do Grupo B inicial:\n")
print(media_base_indice)

cat("\nDesvios-padrão fixos do Grupo B inicial:\n")
print(sd_base_indice)

# ============================================================
# 5. Preparar mapa para AlphaSimR
# ============================================================

mapa_importacao <- mapa_AlphaSimR %>%
  dplyr::select(
    marker,
    chr,
    pos
  ) %>%
  as.data.frame()

mapa_importacao$marker <- as.character(mapa_importacao$marker)
mapa_importacao$chr <- as.integer(mapa_importacao$chr)
mapa_importacao$pos <- as.numeric(mapa_importacao$pos)

# ============================================================
# 6. Importar população base no AlphaSimR
# ============================================================

cat("\nImportando população base no AlphaSimR...\n")

genoma_fundador <- importInbredGeno(
  geno = matriz_marcadores,
  genMap = mapa_importacao
)

parametros_simulacao <- SimParam$new(genoma_fundador)

parametros_simulacao$setTrackPed(TRUE)

parametros_simulacao$addSnpChipByName(
  markers = colnames(matriz_marcadores)
)

pop_base <- newPop(
  rawPop = genoma_fundador,
  simParam = parametros_simulacao
)

pop_base@id <- rownames(matriz_marcadores)

cat("\nPopulação base criada:\n")
cat("Indivíduos:", pop_base@nInd, "\n")
cat("Cromossomos:", pop_base@nChr, "\n")

# ============================================================
# 7. Funções auxiliares
# ============================================================

calc_indice_MT <- function(gv) {
  
  traits <- c("AP", "AE", "FF", "FM")
  
  gv_traits <- gv[, traits, drop = FALSE]
  
  gv_padronizado <- sweep(
    gv_traits,
    2,
    media_base_indice[traits],
    "-"
  )
  
  gv_padronizado <- sweep(
    gv_padronizado,
    2,
    sd_base_indice[traits],
    "/"
  )
  
  indice <- rowMeans(
    gv_padronizado,
    na.rm = TRUE
  )
  
  return(as.numeric(indice))
}

score_por_cenario <- function(gv, cenario) {
  
  if (cenario %in% c("AP", "AE", "FF", "FM")) {
    return(as.numeric(gv[, cenario]))
  }
  
  if (cenario == "MT") {
    return(calc_indice_MT(gv))
  }
  
  if (cenario == "Culling") {
    return(calc_indice_MT(gv))
  }
  
  stop("Cenário não reconhecido.")
}

adicionar_erro <- function(gv, h2_estagio) {
  
  pheno <- gv
  
  for (tr in colnames(gv)) {
    
    varG <- stats::var(gv[, tr], na.rm = TRUE)
    
    if (!is.finite(varG) || varG == 0) {
      varE <- 0
    } else {
      varE <- varG * (1 - h2_estagio[tr]) / h2_estagio[tr]
    }
    
    pheno[, tr] <- gv[, tr] + stats::rnorm(
      n = nrow(gv),
      mean = 0,
      sd = sqrt(varE)
    )
  }
  
  return(pheno)
}

selecionar_por_cenario <- function(df, n_sel, cenario) {
  
  if (cenario %in% c("AP", "AE", "FF", "FM", "MT")) {
    
    df_sel <- df %>%
      dplyr::arrange(dplyr::desc(Score_selecao)) %>%
      dplyr::slice_head(n = n_sel)
    
    return(df_sel)
  }
  
  if (cenario == "Culling") {
    
    n1 <- max(n_sel, ceiling(nrow(df) * 0.70))
    n2 <- max(n_sel, ceiling(n1 * 0.70))
    n3 <- max(n_sel, ceiling(n2 * 0.70))
    
    df_sel <- df %>%
      dplyr::arrange(dplyr::desc(AP)) %>%
      dplyr::slice_head(n = n1) %>%
      dplyr::arrange(dplyr::desc(AE)) %>%
      dplyr::slice_head(n = n2) %>%
      dplyr::arrange(dplyr::desc(FF)) %>%
      dplyr::slice_head(n = n3) %>%
      dplyr::arrange(dplyr::desc(FM)) %>%
      dplyr::slice_head(n = n_sel)
    
    return(df_sel)
  }
  
  stop("Cenário não reconhecido.")
}

calcular_gv_testcross_media <- function(
    genotipos_DH,
    testador_ids,
    matriz_pool,
    matriz_efeitos
) {
  
  lista_gv <- vector("list", length(testador_ids))
  
  for (i in seq_along(testador_ids)) {
    
    id_testador <- testador_ids[i]
    
    geno_testador <- matriz_pool[id_testador, , drop = TRUE]
    
    geno_hibrido <- sweep(
      genotipos_DH,
      2,
      geno_testador,
      "+"
    ) / 2
    
    gv_testador <- as.matrix(geno_hibrido %*% matriz_efeitos)
    
    colnames(gv_testador) <- colnames(matriz_efeitos)
    
    lista_gv[[i]] <- gv_testador
  }
  
  gv_medio <- Reduce("+", lista_gv) / length(lista_gv)
  
  colnames(gv_medio) <- colnames(matriz_efeitos)
  
  return(gv_medio)
}

avaliar_estagio <- function(
    genotipos_DH,
    ids_DH,
    cross_ids,
    testador_ids,
    matriz_pool,
    matriz_efeitos,
    cenario,
    estagio,
    n_sel
) {
  
  gv <- calcular_gv_testcross_media(
    genotipos_DH = genotipos_DH,
    testador_ids = testador_ids,
    matriz_pool = matriz_pool,
    matriz_efeitos = matriz_efeitos
  )
  
  h2_estagio <- h2_traits * fator_h2_estagio[estagio]
  
  h2_estagio[h2_estagio >= 0.99] <- 0.99
  
  pheno <- adicionar_erro(
    gv = gv,
    h2_estagio = h2_estagio
  )
  
  df <- as.data.frame(pheno)
  
  df$ID <- ids_DH
  df$Cross_ID <- cross_ids
  
  df <- df %>%
    dplyr::select(
      ID,
      Cross_ID,
      AP,
      AE,
      FF,
      FM
    )
  
  df$Score_selecao <- score_por_cenario(
    gv = as.matrix(df[, c("AP", "AE", "FF", "FM")]),
    cenario = cenario
  )
  
  df_sel <- selecionar_por_cenario(
    df = df,
    n_sel = n_sel,
    cenario = cenario
  )
  
  return(list(
    selecionados = df_sel,
    gv_verdadeiro = gv,
    pheno = pheno
  ))
}

selecionar_cruzamentos_culling <- function(candidatos, n_cruzamentos) {
  
  n1 <- max(n_cruzamentos, ceiling(nrow(candidatos) * 0.70))
  n2 <- max(n_cruzamentos, ceiling(n1 * 0.70))
  n3 <- max(n_cruzamentos, ceiling(n2 * 0.70))
  
  candidatos_sel <- candidatos %>%
    dplyr::arrange(dplyr::desc(UC_AP)) %>%
    dplyr::slice_head(n = n1) %>%
    dplyr::arrange(dplyr::desc(UC_AE)) %>%
    dplyr::slice_head(n = n2) %>%
    dplyr::arrange(dplyr::desc(UC_FF)) %>%
    dplyr::slice_head(n = n3) %>%
    dplyr::arrange(dplyr::desc(UC_FM)) %>%
    dplyr::slice_head(n = n_cruzamentos)
  
  return(candidatos_sel)
}

calcular_uc_cruzamentos <- function(
    parentais_ids,
    matriz_pool,
    matriz_efeitos,
    cenario,
    n_cruzamentos,
    max_candidatos
) {
  
  combinacoes <- utils::combn(parentais_ids, 2)
  
  candidatos <- data.frame(
    Parent1 = combinacoes[1, ],
    Parent2 = combinacoes[2, ],
    stringsAsFactors = FALSE
  )
  
  if (nrow(candidatos) > max_candidatos) {
    candidatos <- candidatos[
      sample(seq_len(nrow(candidatos)), max_candidatos),
      ,
      drop = FALSE
    ]
  }
  
  geno_p1 <- matriz_pool[candidatos$Parent1, , drop = FALSE]
  geno_p2 <- matriz_pool[candidatos$Parent2, , drop = FALSE]
  
  gv_pool <- as.matrix(
    matriz_pool[parentais_ids, , drop = FALSE] %*% matriz_efeitos
  )
  
  colnames(gv_pool) <- colnames(matriz_efeitos)
  
  gv_p1 <- gv_pool[candidatos$Parent1, , drop = FALSE]
  gv_p2 <- gv_pool[candidatos$Parent2, , drop = FALSE]
  
  media_cruzamento <- (gv_p1 + gv_p2) / 2
  
  efeito2 <- matriz_efeitos^2
  
  sd_cruzamento <- matrix(
    NA_real_,
    nrow = nrow(candidatos),
    ncol = ncol(matriz_efeitos)
  )
  
  colnames(sd_cruzamento) <- colnames(matriz_efeitos)
  
  for (i in seq_len(nrow(candidatos))) {
    
    segregantes <- geno_p1[i, ] != geno_p2[i, ]
    
    if (sum(segregantes) == 0) {
      sd_cruzamento[i, ] <- 0
    } else {
      var_i <- colSums(efeito2[segregantes, , drop = FALSE])
      sd_cruzamento[i, ] <- sqrt(var_i)
    }
  }
  
  uc_trait <- media_cruzamento + i_uc * sd_cruzamento
  
  colnames(uc_trait) <- colnames(matriz_efeitos)
  
  candidatos$UC_AP <- uc_trait[, "AP"]
  candidatos$UC_AE <- uc_trait[, "AE"]
  candidatos$UC_FF <- uc_trait[, "FF"]
  candidatos$UC_FM <- uc_trait[, "FM"]
  candidatos$UC_MT <- calc_indice_MT(uc_trait)
  
  if (cenario %in% c("AP", "AE", "FF", "FM")) {
    
    candidatos$Y <- candidatos[[paste0("UC_", cenario)]]
    
    candidatos <- candidatos %>%
      dplyr::arrange(dplyr::desc(Y)) %>%
      dplyr::slice_head(n = n_cruzamentos)
  }
  
  if (cenario == "MT") {
    
    candidatos$Y <- candidatos$UC_MT
    
    candidatos <- candidatos %>%
      dplyr::arrange(dplyr::desc(Y)) %>%
      dplyr::slice_head(n = n_cruzamentos)
  }
  
  if (cenario == "Culling") {
    
    candidatos <- selecionar_cruzamentos_culling(
      candidatos = candidatos,
      n_cruzamentos = n_cruzamentos
    )
    
    candidatos$Y <- candidatos$UC_MT
  }
  
  candidatos$K <- NA_real_
  
  candidatos <- candidatos %>%
    dplyr::mutate(
      cruzamento_id = paste(Parent1, Parent2, sep = "_")
    ) %>%
    dplyr::select(
      Parent1,
      Parent2,
      Y,
      K,
      cruzamento_id,
      UC_AP,
      UC_AE,
      UC_FF,
      UC_FM,
      UC_MT
    )
  
  return(candidatos)
}

juntar_populacoes <- function(pop1, pop2) {
  
  pop_junta <- mergePops(
    list(
      pop1,
      pop2
    )
  )
  
  return(pop_junta)
}

# ============================================================
# 8. Objetos para guardar resultados
# ============================================================

resultados_ciclos <- list()
resultados_elites <- list()
resultados_cruzamentos <- list()
resultados_pool <- list()
resultados_genotipos_elites <- list()

contador_ciclos <- 1
contador_elites <- 1
contador_cruzamentos <- 1
contador_pool <- 1
contador_genotipos_elites <- 1

# ============================================================
# 9. Loop principal da simulação
# ============================================================

cat("\nIniciando simulação recorrente v2...\n")

for (cenario in cenarios) {
  
  cat("\n============================================================\n")
  cat("CENÁRIO:", cenario, "\n")
  cat("============================================================\n")
  
  cruzamentos_iniciais_cenario <- cruzamentos_iniciais_por_cenario[[cenario]]
  
  if (is.null(cruzamentos_iniciais_cenario)) {
    stop(paste("Não há cruzamentos iniciais para o cenário", cenario))
  }
  
  for (rep_i in seq_len(n_rep)) {
    
    cat("\n--- Repetição", rep_i, "de", n_rep, "| Cenário", cenario, "---\n")
    
    set.seed(1000 + rep_i + which(cenarios == cenario) * 100)
    
    pop_trabalho <- pop_base
    matriz_pool <- matriz_marcadores
    
    parental_pool_ids <- parentais_grupo_B
    
    for (ciclo in seq_len(n_cycles)) {
      
      cat("\nCiclo", ciclo, "| Cenário", cenario, "| Repetição", rep_i, "\n")
      
      # --------------------------------------------------------
      # 9.1 Definir cruzamentos do ciclo
      # --------------------------------------------------------
      
      if (ciclo == 1) {
        
        cruzamentos_ciclo <- cruzamentos_iniciais_cenario %>%
          dplyr::arrange(Ordem_no_cenario) %>%
          dplyr::slice_head(n = n_cruzamentos_por_ciclo) %>%
          dplyr::select(
            Parent1,
            Parent2,
            Y,
            K,
            cruzamento_id,
            dplyr::any_of(c("UC_AP", "UC_AE", "UC_FF", "UC_FM", "UC_MT"))
          )
        
        metodo_cruzamento <- "SimpleMating_inicial"
        
      } else {
        
        cruzamentos_ciclo <- calcular_uc_cruzamentos(
          parentais_ids = parental_pool_ids,
          matriz_pool = matriz_pool,
          matriz_efeitos = matriz_efeitos,
          cenario = cenario,
          n_cruzamentos = n_cruzamentos_por_ciclo,
          max_candidatos = max_candidatos_uc
        )
        
        metodo_cruzamento <- "UC_analogico"
      }
      
      cruzamentos_ciclo$Cenario <- cenario
      cruzamentos_ciclo$Rep <- rep_i
      cruzamentos_ciclo$Ciclo <- ciclo
      cruzamentos_ciclo$Metodo <- metodo_cruzamento
      
      resultados_cruzamentos[[contador_cruzamentos]] <- cruzamentos_ciclo
      contador_cruzamentos <- contador_cruzamentos + 1
      
      # --------------------------------------------------------
      # 9.2 Fazer cruzamentos no AlphaSimR
      # --------------------------------------------------------
      
      ids_pop <- pop_trabalho@id
      
      indice_p1 <- match(cruzamentos_ciclo$Parent1, ids_pop)
      indice_p2 <- match(cruzamentos_ciclo$Parent2, ids_pop)
      
      if (any(is.na(indice_p1)) || any(is.na(indice_p2))) {
        stop("Algum parental do ciclo não foi encontrado em pop_trabalho.")
      }
      
      plano_indices <- cbind(
        indice_p1,
        indice_p2
      )
      
      pop_F1 <- makeCross(
        pop = pop_trabalho,
        crossPlan = plano_indices,
        simParam = parametros_simulacao
      )
      
      pop_DH <- makeDH(
        pop = pop_F1,
        nDH = n_DH_por_cruzamento,
        simParam = parametros_simulacao
      )
      
      ids_DH <- paste0(
        cenario,
        "_R",
        rep_i,
        "_C",
        ciclo,
        "_DH",
        seq_len(pop_DH@nInd)
      )
      
      pop_DH@id <- ids_DH
      
      genotipos_DH <- pullSnpGeno(
        pop = pop_DH,
        simParam = parametros_simulacao
      )
      
      rownames(genotipos_DH) <- ids_DH
      
      cross_ids_DH <- rep(
        cruzamentos_ciclo$cruzamento_id,
        each = n_DH_por_cruzamento
      )
      
      # --------------------------------------------------------
      # 9.3 Avaliação TC1
      # --------------------------------------------------------
      
      avaliacao_TC1 <- avaliar_estagio(
        genotipos_DH = genotipos_DH,
        ids_DH = ids_DH,
        cross_ids = cross_ids_DH,
        testador_ids = testador_ids,
        matriz_pool = matriz_pool,
        matriz_efeitos = matriz_efeitos,
        cenario = cenario,
        estagio = "TC1",
        n_sel = n_TC1
      )
      
      selecionados_TC1 <- avaliacao_TC1$selecionados
      
      genotipos_TC1 <- genotipos_DH[
        selecionados_TC1$ID,
        ,
        drop = FALSE
      ]
      
      cross_ids_TC1 <- selecionados_TC1$Cross_ID
      
      # --------------------------------------------------------
      # 9.4 Avaliação TC2
      # --------------------------------------------------------
      
      avaliacao_TC2 <- avaliar_estagio(
        genotipos_DH = genotipos_TC1,
        ids_DH = rownames(genotipos_TC1),
        cross_ids = cross_ids_TC1,
        testador_ids = testador_ids,
        matriz_pool = matriz_pool,
        matriz_efeitos = matriz_efeitos,
        cenario = cenario,
        estagio = "TC2",
        n_sel = n_TC2
      )
      
      selecionados_TC2 <- avaliacao_TC2$selecionados
      
      genotipos_TC2 <- genotipos_TC1[
        selecionados_TC2$ID,
        ,
        drop = FALSE
      ]
      
      cross_ids_TC2 <- selecionados_TC2$Cross_ID
      
      # --------------------------------------------------------
      # 9.5 Avaliação TC3
      # --------------------------------------------------------
      
      avaliacao_TC3 <- avaliar_estagio(
        genotipos_DH = genotipos_TC2,
        ids_DH = rownames(genotipos_TC2),
        cross_ids = cross_ids_TC2,
        testador_ids = testador_ids,
        matriz_pool = matriz_pool,
        matriz_efeitos = matriz_efeitos,
        cenario = cenario,
        estagio = "TC3",
        n_sel = n_TC3
      )
      
      selecionados_TC3 <- avaliacao_TC3$selecionados
      
      # --------------------------------------------------------
      # 9.6 Selecionar elites
      # --------------------------------------------------------
      
      elites_ciclo <- selecionar_por_cenario(
        df = selecionados_TC3,
        n_sel = n_elite,
        cenario = cenario
      )
      
      elite_ids_antigos <- elites_ciclo$ID
      
      elite_ids_novos <- paste0(
        cenario,
        "_R",
        rep_i,
        "_C",
        ciclo,
        "_E",
        seq_len(n_elite)
      )
      
      genotipos_elite <- genotipos_DH[
        elite_ids_antigos,
        ,
        drop = FALSE
      ]
      
      rownames(genotipos_elite) <- elite_ids_novos
      
      matriz_pool <- rbind(
        matriz_pool,
        genotipos_elite
      )
      
      resultados_genotipos_elites[[contador_genotipos_elites]] <- genotipos_elite
      contador_genotipos_elites <- contador_genotipos_elites + 1
      
      indices_elite_pop <- match(
        elite_ids_antigos,
        pop_DH@id
      )
      
      pop_elite <- pop_DH[indices_elite_pop]
      pop_elite@id <- elite_ids_novos
      
      pop_trabalho <- juntar_populacoes(
        pop1 = pop_trabalho,
        pop2 = pop_elite
      )
      
      elites_ciclo$Elite_ID <- elite_ids_novos
      elites_ciclo$ID_original_DH <- elite_ids_antigos
      elites_ciclo$Cenario <- cenario
      elites_ciclo$Rep <- rep_i
      elites_ciclo$Ciclo <- ciclo
      
      resultados_elites[[contador_elites]] <- elites_ciclo
      contador_elites <- contador_elites + 1
      
      # --------------------------------------------------------
      # 9.7 Substituir piores parentais do pool
      # --------------------------------------------------------
      
      gv_parentais <- as.matrix(
        matriz_pool[parental_pool_ids, , drop = FALSE] %*% matriz_efeitos
      )
      
      colnames(gv_parentais) <- colnames(matriz_efeitos)
      
      score_parentais <- score_por_cenario(
        gv = gv_parentais,
        cenario = cenario
      )
      
      tabela_parentais <- data.frame(
        Linhagem = parental_pool_ids,
        Score = score_parentais,
        stringsAsFactors = FALSE
      )
      
      parentais_removidos <- tabela_parentais %>%
        dplyr::arrange(Score) %>%
        dplyr::slice_head(n = n_elite) %>%
        dplyr::pull(Linhagem)
      
      parental_pool_ids <- setdiff(
        parental_pool_ids,
        parentais_removidos
      )
      
      parental_pool_ids <- c(
        parental_pool_ids,
        elite_ids_novos
      )
      
      if (length(parental_pool_ids) != length(parentais_grupo_B)) {
        stop("O tamanho do pool parental mudou. Verificar substituição.")
      }
      
      # --------------------------------------------------------
      # 9.8 Guardar resumo do ciclo
      # --------------------------------------------------------
      
      gv_elites <- as.matrix(
        matriz_pool[elite_ids_novos, , drop = FALSE] %*% matriz_efeitos
      )
      
      colnames(gv_elites) <- colnames(matriz_efeitos)
      
      score_elites <- score_por_cenario(
        gv = gv_elites,
        cenario = cenario
      )
      
      gv_pool_atual <- as.matrix(
        matriz_pool[parental_pool_ids, , drop = FALSE] %*% matriz_efeitos
      )
      
      colnames(gv_pool_atual) <- colnames(matriz_efeitos)
      
      score_pool_atual <- score_por_cenario(
        gv = gv_pool_atual,
        cenario = cenario
      )
      
      resumo_ciclo <- data.frame(
        Cenario = cenario,
        Rep = rep_i,
        Ciclo = ciclo,
        Metodo_cruzamento = metodo_cruzamento,
        n_parentais_pool = length(parental_pool_ids),
        n_cruzamentos = nrow(cruzamentos_ciclo),
        n_DH = nrow(genotipos_DH),
        n_TC1 = nrow(selecionados_TC1),
        n_TC2 = nrow(selecionados_TC2),
        n_TC3 = nrow(selecionados_TC3),
        n_elite = nrow(elites_ciclo),
        Media_score_elites = mean(score_elites),
        Melhor_score_elite = max(score_elites),
        Media_score_pool = mean(score_pool_atual),
        Melhor_score_pool = max(score_pool_atual),
        Parentais_removidos = paste(parentais_removidos, collapse = ";"),
        Elites_adicionados = paste(elite_ids_novos, collapse = ";"),
        stringsAsFactors = FALSE
      )
      
      resultados_ciclos[[contador_ciclos]] <- resumo_ciclo
      contador_ciclos <- contador_ciclos + 1
      
      resumo_pool <- data.frame(
        Cenario = cenario,
        Rep = rep_i,
        Ciclo = ciclo,
        Linhagem = parental_pool_ids,
        stringsAsFactors = FALSE
      )
      
      resultados_pool[[contador_pool]] <- resumo_pool
      contador_pool <- contador_pool + 1
      
      cat(
        "Resumo ciclo:",
        "média elites =", round(mean(score_elites), 4),
        "| melhor elite =", round(max(score_elites), 4),
        "| média pool =", round(mean(score_pool_atual), 4),
        "\n"
      )
    }
  }
}

# ============================================================
# 10. Consolidar resultados
# ============================================================

cat("\nConsolidando resultados...\n")

tabela_ciclos <- dplyr::bind_rows(resultados_ciclos)
tabela_elites <- dplyr::bind_rows(resultados_elites)
tabela_cruzamentos <- dplyr::bind_rows(resultados_cruzamentos)
tabela_pool <- dplyr::bind_rows(resultados_pool)

matriz_genotipos_elites <- do.call(
  rbind,
  resultados_genotipos_elites
)

# ============================================================
# 11. Conferências finais
# ============================================================

cat("\nDimensões finais:\n")
cat("tabela_ciclos:", dim(tabela_ciclos), "\n")
cat("tabela_elites:", dim(tabela_elites), "\n")
cat("tabela_cruzamentos:", dim(tabela_cruzamentos), "\n")
cat("tabela_pool:", dim(tabela_pool), "\n")
cat("matriz_genotipos_elites:", dim(matriz_genotipos_elites), "\n")

cat("\nCiclos por cenário:\n")
print(table(tabela_ciclos$Cenario, tabela_ciclos$Ciclo))

cat("\nRepetições por cenário:\n")
print(table(tabela_ciclos$Cenario, tabela_ciclos$Rep))

valores_elites <- sort(unique(as.vector(matriz_genotipos_elites)))

cat("\nValores encontrados nos genótipos dos elites:\n")
print(valores_elites)

if (!all(valores_elites %in% c(0, 2))) {
  stop("Os genótipos dos elites possuem valores diferentes de 0 e 2.")
}

if (nrow(tabela_ciclos) != length(cenarios) * n_rep * n_cycles) {
  stop("Número inesperado de linhas em tabela_ciclos.")
}

if (nrow(tabela_elites) != length(cenarios) * n_rep * n_cycles * n_elite) {
  stop("Número inesperado de linhas em tabela_elites.")
}

# ============================================================
# 12. Salvar objeto final
# ============================================================

resultados_simulacao_v2 <- list(
  
  parametros = list(
    n_rep = n_rep,
    n_cycles = n_cycles,
    cenarios = cenarios,
    n_cruzamentos_por_ciclo = n_cruzamentos_por_ciclo,
    n_DH_por_cruzamento = n_DH_por_cruzamento,
    n_TC1 = n_TC1,
    n_TC2 = n_TC2,
    n_TC3 = n_TC3,
    n_elite = n_elite,
    prop_sel_uc = prop_sel_uc,
    i_uc = i_uc,
    max_candidatos_uc = max_candidatos_uc,
    h2_traits = h2_traits,
    fator_h2_estagio = fator_h2_estagio,
    media_base_indice = media_base_indice,
    sd_base_indice = sd_base_indice,
    base_indice_MT = "Grupo B inicial"
  ),
  
  tabela_ciclos = tabela_ciclos,
  
  tabela_elites = tabela_elites,
  
  tabela_cruzamentos = tabela_cruzamentos,
  
  tabela_pool = tabela_pool,
  
  matriz_genotipos_elites = matriz_genotipos_elites,
  
  testadores = testadores,
  
  parentais_grupo_B_inicial = parentais_grupo_B,
  
  observacao = "Simulacao recorrente v2 com dados reais. MT e Culling usam indice com padronizacao fixa baseada no Grupo B inicial. Genotipos dos elites foram salvos para analise de diversidade. Maior score = melhor; efeitos ja orientados para reducao."
)

saveRDS(
  resultados_simulacao_v2,
  "resultados_simulacao_recorrente_alphasimr_v2.rds"
)

# ============================================================
# 13. Mensagem final
# ============================================================

cat("\n============================================================\n")
cat("SCRIPT 04 v2 CONCLUÍDO COM SUCESSO\n")
cat("============================================================\n")

cat("\nArquivo salvo:\n")
cat("- resultados_simulacao_recorrente_alphasimr_v2.rds\n")

cat("\nResumo:\n")
cat("Cenários:", paste(cenarios, collapse = ", "), "\n")
cat("Repetições:", n_rep, "\n")
cat("Ciclos:", n_cycles, "\n")
cat("Elites salvos:", nrow(matriz_genotipos_elites), "\n")

cat("\nEsse RDS está pronto para os scripts de resultado:\n")
