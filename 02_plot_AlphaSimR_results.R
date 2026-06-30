# ============================================================
# 02_plot_AlphaSimR_results.R
#
# Outputs:
# - Figures and summary tables saved in graficos_AlphaSimR/
# ============================================================

rm(list = ls())

library(tidyr)
library(ggplot2)
library(dplyr)

dir.create("graficos_AlphaSimR", showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load results
# ------------------------------------------------------------

cycle_summary <- read.csv("AlphaSimR_cycle_summary_recorrente.csv")
elite_summary <- read.csv("AlphaSimR_elites_recorrente.csv")
TC_selected_all <- read.csv("AlphaSimR_TC_selected_recorrente.csv")
crosses_used_all <- read.csv("AlphaSimR_crosses_used_recorrente.csv")

scenario_order <- c("AP", "AE", "FF", "FM", "MT", "Culling")

cycle_summary$scenario <- factor(cycle_summary$scenario, levels = scenario_order)
elite_summary$scenario <- factor(elite_summary$scenario, levels = scenario_order)
TC_selected_all$scenario <- factor(TC_selected_all$scenario, levels = scenario_order)
crosses_used_all$scenario <- factor(crosses_used_all$scenario, levels = scenario_order)

# ------------------------------------------------------------
# 2. Save function
# ------------------------------------------------------------

salvar_grafico <- function(plot_obj, nome, width = 8.5, height = 5.2) {
  
  ggsave(
    filename = paste0("graficos_AlphaSimR/", nome, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
  
  ggsave(
    filename = paste0("graficos_AlphaSimR/", nome, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height
  )
}

# ------------------------------------------------------------
# 3. Accumulated gain relative to cycle 1
# ------------------------------------------------------------

ganho_relativo <- cycle_summary %>%
  dplyr::group_by(scenario) %>%
  dplyr::arrange(cycle, .by_group = TRUE) %>%
  dplyr::mutate(
    base_ciclo1 = mean_elite[cycle == min(cycle)][1],
    ganho_acumulado = mean_elite - base_ciclo1
  ) %>%
  dplyr::ungroup()

# ------------------------------------------------------------
# Figure 1: accumulated gain
# ------------------------------------------------------------

fig1 <- ggplot(
  ganho_relativo,
  aes(x = cycle, y = ganho_acumulado, color = scenario, group = scenario)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 1.25) +
  geom_point(size = 2.8) +
  theme_classic(base_size = 15) +
  labs(
    x = "Selection cycle",
    y = "Accumulated gain of selected elites",
    color = "Scenario"
  )

print(fig1)

salvar_grafico(
  fig1,
  "01_AlphaSimR_accumulated_gain_six_scenarios",
  width = 9,
  height = 5.5
)

# ------------------------------------------------------------
# Figure 2: mean elite per cycle
# ------------------------------------------------------------

fig2 <- ggplot(
  cycle_summary,
  aes(x = cycle, y = mean_elite, color = scenario, group = scenario)
) +
  geom_line(linewidth = 1.25) +
  geom_point(size = 2.8) +
  theme_classic(base_size = 15) +
  labs(
    x = "Selection cycle",
    y = "Mean simulated phenotype of selected elites",
    color = "Scenario"
  )

print(fig2)

salvar_grafico(
  fig2,
  "02_AlphaSimR_mean_elite_per_cycle",
  width = 9,
  height = 5.5
)

# ------------------------------------------------------------
# Figure 3: testcross selection funnel
# ------------------------------------------------------------

cycle_long <- cycle_summary %>%
  dplyr::select(
    scenario,
    cycle,
    mean_TC1,
    mean_TC2,
    mean_TC3,
    mean_elite
  ) %>%
  tidyr::pivot_longer(
    cols = c(mean_TC1, mean_TC2, mean_TC3, mean_elite),
    names_to = "Stage",
    values_to = "Mean"
  )

cycle_long$Stage <- factor(
  cycle_long$Stage,
  levels = c("mean_TC1", "mean_TC2", "mean_TC3", "mean_elite"),
  labels = c("TC1", "TC2", "TC3", "Elite")
)

fig3 <- ggplot(
  cycle_long,
  aes(x = cycle, y = Mean, color = Stage, group = Stage)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  facet_wrap(~ scenario, scales = "free_y") +
  theme_classic(base_size = 13) +
  labs(
    x = "Selection cycle",
    y = "Mean simulated phenotype",
    color = "Stage"
  )

print(fig3)

salvar_grafico(
  fig3,
  "03_AlphaSimR_testcross_funnel_by_scenario",
  width = 11,
  height = 7
)

# ------------------------------------------------------------
# Figure 4: final accumulated gain
# ------------------------------------------------------------

ganho_final <- ganho_relativo %>%
  dplyr::filter(cycle == max(cycle)) %>%
  dplyr::arrange(dplyr::desc(ganho_acumulado))

ganho_final$scenario <- factor(
  ganho_final$scenario,
  levels = ganho_final$scenario
)

fig4 <- ggplot(
  ganho_final,
  aes(x = scenario, y = ganho_acumulado)
) +
  geom_col(width = 0.7) +
  theme_classic(base_size = 15) +
  labs(
    x = "Scenario",
    y = "Final accumulated gain"
  )

print(fig4)

salvar_grafico(
  fig4,
  "04_AlphaSimR_final_accumulated_gain",
  width = 8,
  height = 5
)

# ------------------------------------------------------------
# Figure 5: selected elites boxplot
# ------------------------------------------------------------

fig5 <- ggplot(
  elite_summary,
  aes(x = scenario, y = phenotype_score)
) +
  geom_boxplot(width = 0.65, outlier.alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.8) +
  theme_classic(base_size = 15) +
  labs(
    x = "Scenario",
    y = "Simulated phenotype of selected elites"
  )

print(fig5)

salvar_grafico(
  fig5,
  "05_AlphaSimR_boxplot_selected_elites",
  width = 8,
  height = 5
)

# ------------------------------------------------------------
# Figure 6: individual trait scenarios vs MT and Culling
# ------------------------------------------------------------

plot_trait_vs_multi <- function(trait_name) {
  
  dados_plot <- ganho_relativo %>%
    dplyr::filter(scenario %in% c(trait_name, "MT", "Culling"))
  
  dados_plot$scenario <- factor(
    dados_plot$scenario,
    levels = c(trait_name, "MT", "Culling")
  )
  
  g <- ggplot(
    dados_plot,
    aes(x = cycle, y = ganho_acumulado, color = scenario, group = scenario)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 1.25) +
    geom_point(size = 2.8) +
    theme_classic(base_size = 15) +
    labs(
      x = "Selection cycle",
      y = "Accumulated gain of selected elites",
      color = "Scenario"
    )
  
  print(g)
  
  salvar_grafico(
    g,
    paste0("06_AlphaSimR_accumulated_gain_", trait_name, "_vs_MT_Culling"),
    width = 8.5,
    height = 5.2
  )
  
  return(g)
}

fig6_AP <- plot_trait_vs_multi("AP")
fig6_AE <- plot_trait_vs_multi("AE")
fig6_FF <- plot_trait_vs_multi("FF")
fig6_FM <- plot_trait_vs_multi("FM")

# ------------------------------------------------------------
# Figure 7: heatmap of mean elite performance
# ------------------------------------------------------------

fig7 <- ggplot(
  cycle_summary,
  aes(x = cycle, y = scenario, fill = mean_elite)
) +
  geom_tile(color = "white") +
  theme_classic(base_size = 15) +
  labs(
    x = "Selection cycle",
    y = "Scenario",
    fill = "Mean elite"
  )

print(fig7)

salvar_grafico(
  fig7,
  "07_AlphaSimR_heatmap_mean_elite",
  width = 8.5,
  height = 5
)

# ------------------------------------------------------------
# Figure 8: number of crosses used
# ------------------------------------------------------------

cross_summary <- crosses_used_all %>%
  dplyr::group_by(scenario, cycle) %>%
  dplyr::summarise(
    n_crosses_used = dplyr::n(),
    .groups = "drop"
  )

fig8 <- ggplot(
  cross_summary,
  aes(x = cycle, y = n_crosses_used, color = scenario, group = scenario)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.6) +
  theme_classic(base_size = 15) +
  labs(
    x = "Selection cycle",
    y = "Number of crosses used",
    color = "Scenario"
  )

print(fig8)

salvar_grafico(
  fig8,
  "08_AlphaSimR_number_crosses_used",
  width = 8.5,
  height = 5.2
)

# ------------------------------------------------------------
# Figure 9: accumulated recycled elites
# ------------------------------------------------------------

fig9 <- ggplot(
  cycle_summary,
  aes(x = cycle, y = n_elites_acumuladas, color = scenario, group = scenario)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.6) +
  theme_classic(base_size = 15) +
  labs(
    x = "Selection cycle",
    y = "Accumulated number of recycled elites",
    color = "Scenario"
  )

print(fig9)

salvar_grafico(
  fig9,
  "09_AlphaSimR_accumulated_recycled_elites",
  width = 8.5,
  height = 5.2
)

# ------------------------------------------------------------
# Save summary tables
# ------------------------------------------------------------

write.csv(
  ganho_relativo,
  "graficos_AlphaSimR/table_AlphaSimR_accumulated_gain_by_cycle.csv",
  row.names = FALSE
)

write.csv(
  ganho_final,
  "graficos_AlphaSimR/table_AlphaSimR_final_accumulated_gain.csv",
  row.names = FALSE
)

write.csv(
  cycle_long,
  "graficos_AlphaSimR/table_AlphaSimR_testcross_funnel.csv",
  row.names = FALSE
)

write.csv(
  cross_summary,
  "graficos_AlphaSimR/table_AlphaSimR_crosses_used_by_cycle.csv",
  row.names = FALSE
)

cat("\nPlot script completed successfully.\n")
cat("Figures saved in: graficos_AlphaSimR\n")
