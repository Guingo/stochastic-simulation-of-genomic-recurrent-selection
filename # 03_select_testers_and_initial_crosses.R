# ============================================================
# 03_select_testers_and_initial_crosses.R

# Select three testers from genetic Group A and select the top
# 50 B x B crosses per scenario to initialize the recurrent
# selection simulation.
#
# Tester selection from Group A:
# - one central line in the PCA space;
# - one line with low PC2 score;
# - one line with high PC2 score.
#
# Initial crosses:
# - top 50 B x B crosses per scenario, ranked by Y.

# Outputs:
# - testadores_grupo_A.rds
# - testadores_grupo_A.csv
# - cruzamentos_programa_B_top50_por_cenario.rds
# - cruzamentos_programa_B_top50_por_cenario.csv
# - resumo_cruzamentos_B_top50_por_cenario.csv
# - PCA_testadores_grupo_A.png
# - PCA_testadores_grupo_A.pdf
# ============================================================

rm(list = ls())

library(dplyr)
library(ggplot2)

set.seed(123)

# ------------------------------------------------------------
# 1. Load inputs
# ------------------------------------------------------------

grupos_linhagens <- readRDS("grupos_linhagens_A_B.rds")
pca_df <- readRDS("pca_grupos_linhagens_A_B.rds")
cruzamentos_programa_B <- readRDS("cruzamentos_programa_B_BxB.rds")

# ------------------------------------------------------------
# 2. Initial checks
# ------------------------------------------------------------

cat("\nNumber of lines per genetic group:\n")
print(table(grupos_linhagens$Grupo))

cat("\nGroups available in PCA data:\n")
print(table(pca_df$Grupo))

cat("\nNumber of B x B crosses per scenario available:\n")
print(table(cruzamentos_programa_B$scenario))

required_cols_pca <- c("Linhagem", "Grupo", "PC1", "PC2")

missing_cols_pca <- required_cols_pca[
  !(required_cols_pca %in% colnames(pca_df))
]

if (length(missing_cols_pca) > 0) {
  stop(
    paste(
      "Missing required columns in pca_df:",
      paste(missing_cols_pca, collapse = ", ")
    )
  )
}

required_cols_crosses <- c("Parent1", "Parent2", "scenario", "Y")

missing_cols_crosses <- required_cols_crosses[
  !(required_cols_crosses %in% colnames(cruzamentos_programa_B))
]

if (length(missing_cols_crosses) > 0) {
  stop(
    paste(
      "Missing required columns in cruzamentos_programa_B:",
      paste(missing_cols_crosses, collapse = ", ")
    )
  )
}

if (!("cross_id" %in% colnames(cruzamentos_programa_B))) {
  
  cruzamentos_programa_B <- cruzamentos_programa_B %>%
    dplyr::mutate(
      cross_id = paste(Parent1, Parent2, sep = "_")
    )
}

# ------------------------------------------------------------
# 3. Separate Group A in the PCA space
# ------------------------------------------------------------

grupo_A_pca <- pca_df %>%
  dplyr::filter(Grupo == "A")

cat("\nNumber of lines in Group A:\n")
print(nrow(grupo_A_pca))

if (nrow(grupo_A_pca) < 3) {
  stop("Group A has fewer than three lines. It is not possible to select three testers.")
}

# ------------------------------------------------------------
# 4. Calculate distance to the Group A centroid
# ------------------------------------------------------------

centro_PC1 <- mean(grupo_A_pca$PC1, na.rm = TRUE)
centro_PC2 <- mean(grupo_A_pca$PC2, na.rm = TRUE)

grupo_A_pca <- grupo_A_pca %>%
  dplyr::mutate(
    dist_centro = sqrt(
      (PC1 - centro_PC1)^2 +
        (PC2 - centro_PC2)^2
    )
  )

# ------------------------------------------------------------
# 5. Select three testers from Group A
# ------------------------------------------------------------
# Selection criteria:
# 1. Central line in the PCA space;
# 2. Extreme line with low PC2;
# 3. Extreme line with high PC2.
#
# If duplicated lines are selected, additional lines far from
# the centroid are used to complete three unique testers.

tester_central <- grupo_A_pca %>%
  dplyr::arrange(dist_centro) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(Tipo_testador = "Central")

tester_extremo_baixo <- grupo_A_pca %>%
  dplyr::arrange(PC2) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(Tipo_testador = "Extremo_PC2_baixo")

tester_extremo_alto <- grupo_A_pca %>%
  dplyr::arrange(dplyr::desc(PC2)) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(Tipo_testador = "Extremo_PC2_alto")

testadores_A <- dplyr::bind_rows(
  tester_central,
  tester_extremo_baixo,
  tester_extremo_alto
) %>%
  dplyr::distinct(Linhagem, .keep_all = TRUE)

if (nrow(testadores_A) < 3) {
  
  extras <- grupo_A_pca %>%
    dplyr::filter(!(Linhagem %in% testadores_A$Linhagem)) %>%
    dplyr::arrange(dplyr::desc(dist_centro)) %>%
    dplyr::slice_head(n = 3 - nrow(testadores_A)) %>%
    dplyr::mutate(
      Tipo_testador = paste0("Extra_", dplyr::row_number())
    )
  
  testadores_A <- dplyr::bind_rows(
    testadores_A,
    extras
  )
}

testadores_A <- testadores_A %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::select(
    Tipo_testador,
    Linhagem,
    Grupo,
    PC1,
    PC2,
    dist_centro
  )

cat("\nSelected testers from Group A:\n")
print(testadores_A)

if (nrow(testadores_A) != 3) {
  stop("The final tester set does not contain exactly three testers.")
}

if (anyDuplicated(testadores_A$Linhagem) > 0) {
  stop("Duplicated testers were found.")
}

# ------------------------------------------------------------
# 6. Select top 50 B x B crosses per scenario
# ------------------------------------------------------------

n_crosses_scenario <- 50

cruzamentos_B_top50 <- cruzamentos_programa_B %>%
  dplyr::group_by(scenario) %>%
  dplyr::arrange(dplyr::desc(Y), .by_group = TRUE) %>%
  dplyr::slice_head(n = n_crosses_scenario) %>%
  dplyr::ungroup()

cat("\nNumber of selected B x B crosses per scenario:\n")
print(table(cruzamentos_B_top50$scenario))

# Check if any scenario has fewer than 50 selected crosses
resumo_top50 <- cruzamentos_B_top50 %>%
  dplyr::group_by(scenario) %>%
  dplyr::summarise(
    n_selected_crosses = dplyr::n(),
    max_Y = max(Y, na.rm = TRUE),
    min_Y = min(Y, na.rm = TRUE),
    mean_Y = mean(Y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(scenario)

cat("\nSummary of selected B x B crosses:\n")
print(resumo_top50)

if (any(resumo_top50$n_selected_crosses < n_crosses_scenario)) {
  
  warning(
    "At least one scenario has fewer than 50 B x B crosses available."
  )
}

# ------------------------------------------------------------
# 7. PCA plot highlighting selected testers
# ------------------------------------------------------------

plot_testadores <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Grupo)
) +
  geom_point(size = 2.5, alpha = 0.65) +
  geom_point(
    data = testadores_A,
    aes(x = PC1, y = PC2),
    inherit.aes = FALSE,
    shape = 8,
    size = 5,
    stroke = 1.3
  ) +
  geom_text(
    data = testadores_A,
    aes(x = PC1, y = PC2, label = Linhagem),
    inherit.aes = FALSE,
    vjust = -1.0,
    size = 4
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Genetic group"
  )

print(plot_testadores)

# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------

saveRDS(
  testadores_A,
  "testadores_grupo_A.rds"
)


saveRDS(
  cruzamentos_B_top50,
  "cruzamentos_programa_B_top50_por_cenario.rds"
)
ggsave(
  filename = "PCA_testadores_grupo_A.png",
  plot = plot_testadores,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = "PCA_testadores_grupo_A.pdf",
  plot = plot_testadores,
  width = 7.5,
  height = 5.5
)

cat("\n============================================================\n")
cat("TESTER AND INITIAL CROSS SELECTION COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nFiles saved:\n")
cat("- testadores_grupo_A.rds\n")
cat("- testadores_grupo_A.csv\n")
cat("- cruzamentos_programa_B_top50_por_cenario.rds\n")
cat("- cruzamentos_programa_B_top50_por_cenario.csv\n")
cat("- resumo_cruzamentos_B_top50_por_cenario.csv\n")
cat("- PCA_testadores_grupo_A.png\n")
cat("- PCA_testadores_grupo_A.pdf\n")
