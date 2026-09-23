library(dplyr)

set.seed(123)

# ============================================================
# 1. Carregar entradas
# ============================================================

entrada <- readRDS("sim_inputs_real.rds")

planos_todos <- entrada$planos_todos

grupos_linhagens <- readRDS(
  "grupos_geneticos_operacionais_PCA_kmeans.rds"
)

testadores <- readRDS(
  "testadores_selecionados_via_PCA.rds"
)

cat("\nDimensão dos planos de cruzamento:\n")
print(dim(planos_todos))

cat("\nColunas dos planos de cruzamento:\n")
print(colnames(planos_todos))

cat("\nDimensão da tabela de grupos:\n")
print(dim(grupos_linhagens))

cat("\nTestadores selecionados:\n")
print(testadores$Linhagem)

# ============================================================
# 2. Padronizar nomes das colunas dos planos
# ============================================================

possiveis_colunas_cenario <- c(
  "Cenario",
  "Cenário",
  "cenario",
  "cenário",
  "Scenario",
  "scenario",
  "Trait",
  "trait",
  "Caracter",
  "caracter",
  "Variavel",
  "variavel"
)

coluna_cenario_encontrada <- intersect(
  possiveis_colunas_cenario,
  colnames(planos_todos)
)

if (!"Cenario" %in% colnames(planos_todos)) {
  
  if (length(coluna_cenario_encontrada) == 0) {
    cat("\nColunas disponíveis em planos_todos:\n")
    print(colnames(planos_todos))
    stop("Não encontrei coluna de cenário em planos_todos.")
  }
  
  names(planos_todos)[
    names(planos_todos) == coluna_cenario_encontrada[1]
  ] <- "Cenario"
}

if (!"Parent1" %in% colnames(planos_todos) && "parent1" %in% colnames(planos_todos)) {
  names(planos_todos)[names(planos_todos) == "parent1"] <- "Parent1"
}

if (!"Parent2" %in% colnames(planos_todos) && "parent2" %in% colnames(planos_todos)) {
  names(planos_todos)[names(planos_todos) == "parent2"] <- "Parent2"
}

cat("\nColunas padronizadas dos planos:\n")
print(colnames(planos_todos))

# ============================================================
# 3. Conferir colunas obrigatórias
# ============================================================

colunas_planos <- c(
  "Cenario",
  "Parent1",
  "Parent2",
  "Y",
  "K"
)

colunas_grupos <- c(
  "Linhagem",
  "Grupo",
  "Grupo_operacional"
)

colunas_faltantes_planos <- setdiff(
  colunas_planos,
  colnames(planos_todos)
)

colunas_faltantes_grupos <- setdiff(
  colunas_grupos,
  colnames(grupos_linhagens)
)

if (length(colunas_faltantes_planos) > 0) {
  cat("\nColunas faltantes em planos_todos:\n")
  print(colunas_faltantes_planos)
  stop("A tabela planos_todos não possui todas as colunas obrigatórias.")
}

if (length(colunas_faltantes_grupos) > 0) {
  cat("\nColunas faltantes em grupos_linhagens:\n")
  print(colunas_faltantes_grupos)
  stop("A tabela de grupos não possui todas as colunas obrigatórias.")
}

# ============================================================
# 4. Preparar tabela de grupos
# ============================================================

if (!"Testador_selecionado" %in% colnames(grupos_linhagens)) {
  grupos_linhagens$Testador_selecionado <- ifelse(
    grupos_linhagens$Linhagem %in% testadores$Linhagem,
    "Sim",
    "Nao"
  )
}

grupos_resumo <- grupos_linhagens %>%
  dplyr::select(
    Linhagem,
    Grupo,
    Grupo_operacional,
    Testador_selecionado
  )

cat("\nResumo dos grupos A/B:\n")
print(table(grupos_resumo$Grupo))

cat("\nResumo dos grupos operacionais:\n")
print(table(grupos_resumo$Grupo_operacional))

# Conferir se os testadores estão no Grupo A
checagem_testadores <- grupos_resumo %>%
  dplyr::filter(Linhagem %in% testadores$Linhagem)

cat("\nGrupo dos testadores selecionados:\n")
print(checagem_testadores)

if (nrow(checagem_testadores) != nrow(testadores)) {
  stop("Nem todos os testadores foram encontrados na tabela de grupos.")
}

if (!all(checagem_testadores$Grupo == "A")) {
  stop("Nem todos os testadores selecionados estão no Grupo A.")
}

# ============================================================
# 5. Classificar cruzamentos por grupo dos parentais
# ============================================================

cat("\nClassificando cruzamentos por grupo genético...\n")

grupos_parental_1 <- grupos_resumo %>%
  dplyr::select(
    Parent1 = Linhagem,
    Grupo_P1 = Grupo,
    Grupo_operacional_P1 = Grupo_operacional
  )

grupos_parental_2 <- grupos_resumo %>%
  dplyr::select(
    Parent2 = Linhagem,
    Grupo_P2 = Grupo,
    Grupo_operacional_P2 = Grupo_operacional
  )

cruzamentos_classificados <- planos_todos %>%
  dplyr::left_join(
    grupos_parental_1,
    by = "Parent1"
  ) %>%
  dplyr::left_join(
    grupos_parental_2,
    by = "Parent2"
  ) %>%
  dplyr::mutate(
    Tipo_cruzamento = dplyr::case_when(
      Grupo_P1 == "A" & Grupo_P2 == "A" ~ "A_x_A",
      Grupo_P1 == "B" & Grupo_P2 == "B" ~ "B_x_B",
      Grupo_P1 != Grupo_P2 ~ "A_x_B",
      TRUE ~ "Nao_classificado"
    )
  )

# ============================================================
# 6. Conferir parentais
# ============================================================

parentais_nao_encontrados <- cruzamentos_classificados %>%
  dplyr::filter(
    is.na(Grupo_P1) | is.na(Grupo_P2)
  )

if (nrow(parentais_nao_encontrados) > 0) {
  
  cat("\nCruzamentos com parentais não encontrados:\n")
  print(parentais_nao_encontrados)
  
  stop("Existem parentais dos planos SimpleMating que não foram encontrados na tabela de grupos.")
}

cat("\nTodos os parentais foram encontrados na tabela de grupos.\n")

cat("\nResumo geral dos tipos de cruzamento:\n")
print(table(cruzamentos_classificados$Tipo_cruzamento))

# ============================================================
# 7. Filtrar cruzamentos B x B
# ============================================================

cruzamentos_BxB <- cruzamentos_classificados %>%
  dplyr::filter(Tipo_cruzamento == "B_x_B")

cat("\nNúmero de cruzamentos B x B por cenário:\n")
print(table(cruzamentos_BxB$Cenario))

# ============================================================
# 8. Selecionar 50 melhores cruzamentos B x B por cenário
# ============================================================

n_cruzamentos_iniciais <- 50

checagem_n_BxB <- cruzamentos_BxB %>%
  dplyr::group_by(Cenario) %>%
  dplyr::summarise(
    n_BxB = dplyr::n(),
    .groups = "drop"
  )

cat("\nNúmero de cruzamentos B x B disponíveis por cenário:\n")
print(checagem_n_BxB)

if (any(checagem_n_BxB$n_BxB < n_cruzamentos_iniciais)) {
  stop("Algum cenário possui menos de 50 cruzamentos B x B disponíveis.")
}

cruzamentos_iniciais_top50 <- cruzamentos_BxB %>%
  dplyr::group_by(Cenario) %>%
  dplyr::arrange(
    dplyr::desc(Y),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(
    n = n_cruzamentos_iniciais
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Cenario) %>%
  dplyr::mutate(
    Ordem_no_cenario = dplyr::row_number()
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Cenario, Ordem_no_cenario)

cat("\nCruzamentos iniciais selecionados por cenário:\n")
print(table(cruzamentos_iniciais_top50$Cenario))

cat("\nPrimeiros cruzamentos selecionados:\n")
print(head(cruzamentos_iniciais_top50, 20))

# ============================================================
# 9. Conferência final dos cruzamentos selecionados
# ============================================================

if (!all(cruzamentos_iniciais_top50$Tipo_cruzamento == "B_x_B")) {
  stop("Nem todos os cruzamentos selecionados são B x B.")
}

if (!all(cruzamentos_iniciais_top50$Grupo_P1 == "B")) {
  stop("Existem Parent1 fora do Grupo B.")
}

if (!all(cruzamentos_iniciais_top50$Grupo_P2 == "B")) {
  stop("Existem Parent2 fora do Grupo B.")
}

cat("\nTodos os cruzamentos iniciais selecionados são B x B.\n")

# ============================================================
# 10. Criar lista por cenário para uso na simulação
# ============================================================

cruzamentos_iniciais_por_cenario <- split(
  cruzamentos_iniciais_top50,
  cruzamentos_iniciais_top50$Cenario
)

cat("\nCenários preparados:\n")
print(names(cruzamentos_iniciais_por_cenario))

# ============================================================
# 11. Salvar somente o objeto necessário para o próximo script
# ============================================================

entradas_cruzamentos_iniciais <- list(
  
  cruzamentos_iniciais_top50 = cruzamentos_iniciais_top50,
  
  cruzamentos_iniciais_por_cenario = cruzamentos_iniciais_por_cenario,
  
  testadores = testadores,
  
  grupos_resumo = grupos_resumo,
  
  checagem_n_BxB = checagem_n_BxB,
  
  n_cruzamentos_iniciais = n_cruzamentos_iniciais,
  
  observacao = "Cruzamentos iniciais B x B selecionados dentro do grupo de selecao recorrente. Grupo A mantido como grupo oposto/testador."
)

saveRDS(
  entradas_cruzamentos_iniciais,
  "entradas_cruzamentos_iniciais_simulacao.rds"
)

# ============================================================
# 12. Mensagem final
# ============================================================

cat("\n============================================================\n")
cat("SCRIPT 03 CONCLUÍDO COM SUCESSO\n")
cat("============================================================\n")

cat("\nArquivo salvo:\n")
cat("- entradas_cruzamentos_iniciais_simulacao.rds\n")

cat("\nResumo dos grupos:\n")
print(table(grupos_resumo$Grupo))

cat("\nResumo dos tipos de cruzamento nos planos SimpleMating:\n")
print(table(cruzamentos_classificados$Tipo_cruzamento))

cat("\nCruzamentos B x B disponíveis por cenário:\n")
print(checagem_n_BxB)

cat("\nCruzamentos iniciais selecionados por cenário:\n")
print(table(cruzamentos_iniciais_top50$Cenario))

cat("\nTestadores que serão usados depois na simulação:\n")
print(testadores$Linhagem)

cat("\nObservação metodológica:\n")
cat("Os cruzamentos iniciais foram selecionados dentro do Grupo B, definido como grupo de seleção recorrente.\n")
cat("O Grupo A foi mantido como grupo oposto/testador, com três testadores representativos selecionados via PCA.\n")
cat("Como os efeitos já estão orientados para redução, maior valor de Y foi considerado melhor.\n")
