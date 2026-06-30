 ============================================================
# 02_classify_crosses_by_genetic_group.R
#
# Purpose:
# Classify candidate crosses according to the genetic groups

# Crosses are classified as:
# - A_x_A
# - B_x_B
# - A_x_B
#
# The B_x_B crosses are retained as the initial recurrent
# selection crossing set.
#
# Inputs:
# - sim_inputs_real.rds
# - grupos_linhagens_A_B.rds
#
# Outputs:
# - cruzamentos_classificados_por_grupo.rds
# - cruzamentos_classificados_por_grupo.csv
# - cruzamentos_programa_B_BxB.rds
# - cruzamentos_programa_B_BxB.csv
# - resumo_cruzamentos_por_grupo.csv
# ============================================================

rm(list = ls())
library(dplyr)
set.seed(123)

# ------------------------------------------------------------
# 1. Load inputs
# ------------------------------------------------------------

sim_inputs_real <- readRDS("sim_inputs_real.rds")
grupos_linhagens <- readRDS("grupos_linhagens_A_B.rds")
plans_by_scenario <- sim_inputs_real$plans_by_scenario


# ------------------------------------------------------------
# 2. Combine crossing plans across scenarios
# ------------------------------------------------------------

cruzamentos_all <- dplyr::bind_rows(
  lapply(
    names(plans_by_scenario),
    function(sc) {
      
      df <- plans_by_scenario[[sc]]
      
      df$scenario <- sc
      
      return(df)
    }
  )
)

# ------------------------------------------------------------
# 3. Standardize required columns
# ------------------------------------------------------------

required_cols <- c("Parent1", "Parent2", "scenario")

missing_cols <- required_cols[
  !(required_cols %in% colnames(cruzamentos_all))
]

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns in crossing plans:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

if (!("cross_id" %in% colnames(cruzamentos_all))) {
  
  cruzamentos_all <- cruzamentos_all %>%
    dplyr::mutate(
      cross_id = paste(Parent1, Parent2, sep = "_")
    )
}

# ------------------------------------------------------------
# 4. Add genetic group of each parent
# ------------------------------------------------------------

grupos_p1 <- grupos_linhagens %>%
  dplyr::rename(
    Parent1 = Linhagem,
    Grupo_P1 = Grupo
  )

grupos_p2 <- grupos_linhagens %>%
  dplyr::rename(
    Parent2 = Linhagem,
    Grupo_P2 = Grupo
  )

cruzamentos_classificados <- cruzamentos_all %>%
  dplyr::left_join(
    grupos_p1,
    by = "Parent1"
  ) %>%
  dplyr::left_join(
    grupos_p2,
    by = "Parent2"
  )

# ------------------------------------------------------------
# 5. Check missing parental group information
# ------------------------------------------------------------

missing_group <- cruzamentos_classificados %>%
  dplyr::filter(
    is.na(Grupo_P1) | is.na(Grupo_P2)
  )

if (nrow(missing_group) > 0) {
  
  cat("\nWARNING: Some crosses have missing parental group information.\n")
  cat("Number of crosses with missing group:\n")
  print(nrow(missing_group))
  
  print(
    head(
      missing_group %>%
        dplyr::select(scenario, Parent1, Parent2, Grupo_P1, Grupo_P2),
      20
    )
  )
  
  stop("Fix missing parental group information before continuing.")
}

# ------------------------------------------------------------
# 6. Classify cross type
# ------------------------------------------------------------

cruzamentos_classificados <- cruzamentos_classificados %>%
  dplyr::mutate(
    Tipo_cruzamento = dplyr::case_when(
      Grupo_P1 == "A" & Grupo_P2 == "A" ~ "A_x_A",
      Grupo_P1 == "B" & Grupo_P2 == "B" ~ "B_x_B",
      Grupo_P1 != Grupo_P2 ~ "A_x_B",
      TRUE ~ NA_character_
    )
  )

cat("\nNumber of crosses by cross type:\n")
print(table(cruzamentos_classificados$Tipo_cruzamento, useNA = "ifany"))

cat("\nNumber of crosses by scenario and cross type:\n")
print(table(
  cruzamentos_classificados$scenario,
  cruzamentos_classificados$Tipo_cruzamento
))

# ------------------------------------------------------------
# 7. Retain B x B crosses for the recurrent selection program
# ------------------------------------------------------------

cruzamentos_programa_B <- cruzamentos_classificados %>%
  dplyr::filter(Tipo_cruzamento == "B_x_B")

cat("\nNumber of B x B crosses per scenario:\n")
print(table(cruzamentos_programa_B$scenario))

# ------------------------------------------------------------
# 8. Create summary table
# ------------------------------------------------------------

resumo_cruzamentos <- cruzamentos_classificados %>%
  dplyr::group_by(scenario, Tipo_cruzamento) %>%
  dplyr::summarise(
    n_crosses = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(scenario, Tipo_cruzamento)

cat("\nCross classification summary:\n")
print(resumo_cruzamentos)

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

saveRDS(cruzamentos_classificados, "cruzamentos_classificados_por_grupo.rds")
saveRDS(cruzamentos_programa_B,"cruzamentos_programa_B_BxB.rds")
