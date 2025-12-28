#!/usr/bin/env Rscript
# =============================================================================
# 06_network_visualization.R
# Interactive network graph of antibiotic-species associations
# Uses visNetwork for interactive HTML output
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(visNetwork)
library(igraph)
library(htmlwidgets)

# =============================================================================
# Configuration
# =============================================================================

project_dir <- "/home/david/projects/abx_metagenomic_correlations"
results_dir <- file.path(project_dir, "results/new_analysis_legacy_data")

# Output directory for network visualizations
network_dir <- file.path(results_dir, "network_graphs")
dir.create(network_dir, showWarnings = FALSE, recursive = TRUE)

# Network parameters
config <- list(
  # Significance threshold
  pval_threshold = 0.1,


  # Minimum methods for robust associations (if using concordance data)
  min_methods = 2,

  # Maximum species per antibiotic to display (prevents overcrowding)
  max_species_per_abx = 15,

  # Minimum absolute effect size to include
  min_effect_size = 0.1,

  # Use robust associations if available, otherwise fall back to ALDEx2

  prefer_robust = TRUE
)

# =============================================================================
# 1. Load Data
# =============================================================================

cat("=== Loading Data ===\n\n")

# Try to load robust associations first (from 3-method concordance analysis)
robust_file <- file.path(results_dir, "individual_antibiotics_v2_with_covariates",
                          "robust_associations.csv")
concordance_file <- file.path(results_dir, "individual_antibiotics_v2_with_covariates",
                               "method_concordance.csv")

# Fall back to ALDEx2-only results
aldex_file <- file.path(results_dir, "individual_antibiotics",
                         "all_antibiotics_combined.csv")

data_source <- NULL
associations <- NULL

if (config$prefer_robust && file.exists(robust_file)) {
  cat("Loading robust associations (2+ methods concordant)...\n")
  associations <- read_csv(robust_file, show_col_types = FALSE)
  data_source <- "robust"
  cat("  Found", nrow(associations), "robust associations\n")

} else if (config$prefer_robust && file.exists(concordance_file)) {
  cat("Loading method concordance data...\n")
  concordance <- read_csv(concordance_file, show_col_types = FALSE)
  associations <- concordance %>%
    filter(n_methods >= config$min_methods)
  data_source <- "concordance"
  cat("  Found", nrow(associations), "associations with", config$min_methods, "+ methods\n")

} else if (file.exists(aldex_file)) {
  cat("Loading ALDEx2 results (robust associations not yet available)...\n")
  cat("  NOTE: This data does NOT adjust for co-administered antibiotics\n")
  aldex_data <- read_csv(aldex_file, show_col_types = FALSE)

  # Filter to significant associations
  associations <- aldex_data %>%
    filter(we.eBH < config$pval_threshold) %>%
    mutate(
      mean_effect = effect,
      direction = ifelse(effect > 0, "increased", "decreased"),
      n_methods = 1,
      methods = "ALDEx2"
    ) %>%
    select(antibiotic, species, mean_effect, direction, n_methods, methods)

  data_source <- "aldex2_no_covariates"
  cat("  Found", nrow(associations), "significant ALDEx2 associations\n")

} else {
  stop("No association data found. Please run the analysis scripts first.")
}

cat("\nData source:", data_source, "\n")

# =============================================================================
# 2. Filter and Prepare Network Data
# =============================================================================

cat("\n=== Preparing Network Data ===\n\n")

# Filter by effect size
associations_filtered <- associations %>%
  filter(abs(mean_effect) >= config$min_effect_size)

cat("After effect size filter (>=", config$min_effect_size, "):",
    nrow(associations_filtered), "associations\n")

# Limit species per antibiotic to prevent overcrowding
associations_top <- associations_filtered %>%
  group_by(antibiotic) %>%
  arrange(desc(abs(mean_effect))) %>%
  slice_head(n = config$max_species_per_abx) %>%
  ungroup()

cat("After limiting to top", config$max_species_per_abx,
    "per antibiotic:", nrow(associations_top), "associations\n")

# Clean up species names for display
associations_top <- associations_top %>%
  mutate(
    species_clean = str_replace_all(species, "\\.", " "),
    species_short = ifelse(
      str_length(species_clean) > 30,
      paste0(str_sub(species_clean, 1, 27), "..."),
      species_clean
    )
  )

# =============================================================================
# 3. Create Network Nodes
# =============================================================================

cat("\n=== Building Network ===\n\n")

# Antibiotic nodes
antibiotics <- unique(associations_top$antibiotic)
abx_nodes <- data.frame(
  id = antibiotics,
  label = antibiotics,
  group = "antibiotic",
  shape = "box",
  color.background = "#E74C3C",
  color.border = "#C0392B",
  color.highlight.background = "#FF6B6B",
  color.highlight.border = "#E74C3C",
  font.color = "white",
  font.size = 16,
  font.bold = TRUE,
  size = 30,
  title = paste0("<b>", antibiotics, "</b><br>",
                 "Antibiotic<br>",
                 sapply(antibiotics, function(a) {
                   n <- sum(associations_top$antibiotic == a)
                   paste0(n, " associated species")
                 }))
)

# Species nodes
species_list <- unique(associations_top$species)
species_data <- associations_top %>%
  group_by(species, species_clean, species_short) %>%
  summarise(
    n_antibiotics = n(),
    antibiotics_list = paste(antibiotic, collapse = ", "),
    avg_effect = mean(mean_effect),
    max_effect = max(abs(mean_effect)),
    .groups = "drop"
  )

# Color species by predominant direction
species_nodes <- data.frame(
  id = species_data$species,
  label = species_data$species_short,
  group = "species",
  shape = "dot",
  size = 10 + species_data$n_antibiotics * 5,  # Larger if connected to more antibiotics
  color.background = ifelse(species_data$avg_effect > 0, "#27AE60", "#3498DB"),
  color.border = ifelse(species_data$avg_effect > 0, "#1E8449", "#2980B9"),
  color.highlight.background = ifelse(species_data$avg_effect > 0, "#2ECC71", "#5DADE2"),
  color.highlight.border = ifelse(species_data$avg_effect > 0, "#27AE60", "#3498DB"),
  font.size = 12,
  title = paste0(
    "<b>", species_data$species_clean, "</b><br>",
    "Connected to ", species_data$n_antibiotics, " antibiotics<br>",
    "Antibiotics: ", species_data$antibiotics_list, "<br>",
    "Avg effect: ", round(species_data$avg_effect, 3)
  )
)

# Combine nodes
nodes <- bind_rows(abx_nodes, species_nodes)

cat("Nodes: ", nrow(nodes), " (", length(antibiotics), " antibiotics, ",
    length(species_list), " species)\n", sep = "")

# =============================================================================
# 4. Create Network Edges
# =============================================================================

# Edge colors based on direction
edges <- data.frame(
  from = associations_top$antibiotic,
  to = associations_top$species,
  value = abs(associations_top$mean_effect) * 3,  # Edge width
  color = ifelse(associations_top$direction == "increased",
                 "#27AE60",  # Green for increased
                 "#E74C3C"), # Red for decreased
  arrows = "to",
  smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.2),
  title = paste0(
    "<b>", associations_top$antibiotic, " → ", associations_top$species_clean, "</b><br>",
    "Direction: ", associations_top$direction, "<br>",
    "Effect size: ", round(associations_top$mean_effect, 3), "<br>",
    ifelse(data_source != "aldex2",
           paste0("Methods: ", associations_top$methods, " (", associations_top$n_methods, ")"),
           "Method: ALDEx2")
  ),
  dashes = ifelse(associations_top$n_methods >= 3, FALSE,
                  ifelse(associations_top$n_methods == 2, FALSE, TRUE))
)

cat("Edges:", nrow(edges), "\n")

# =============================================================================
# 5. Create Interactive Network
# =============================================================================

cat("\n=== Generating Interactive Visualization ===\n\n")

# Dynamic title based on data source
if (data_source == "aldex2_no_covariates") {
  main_title <- "Antibiotic-Microbiome Interaction Network (PRELIMINARY)"
  sub_title <- paste0("Top ", config$max_species_per_abx,
                      " species per antibiotic | ",
                      "ALDEx2 only - NO covariate adjustment for co-administered antibiotics | ",
                      "Effect threshold: ", config$min_effect_size)
  title_color <- "#E67E22"  # Orange for preliminary
} else {

  main_title <- "Antibiotic-Microbiome Interaction Network"
  sub_title <- paste0("Top ", config$max_species_per_abx,
                      " species associations per antibiotic | ",
                      "Data: ", data_source, " | ",
                      "Effect threshold: ", config$min_effect_size)
  title_color <- "#2C3E50"
}

network <- visNetwork(nodes, edges,
                       main = list(
                         text = main_title,
                         style = paste0("font-family: Arial; font-size: 24px; font-weight: bold; color: ", title_color, ";")
                       ),
                       submain = list(
                         text = sub_title,
                         style = "font-family: Arial; font-size: 14px; color: #7F8C8D;"
                       ),
                       footer = list(
                         text = "Green edges = increased abundance | Red edges = decreased abundance | Node size = # connections",
                         style = "font-family: Arial; font-size: 12px; color: #95A5A6;"
                       ),
                       width = "100%",
                       height = "900px") %>%

  # Physics for layout
  visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(
      gravitationalConstant = -100,
      centralGravity = 0.01,
      springLength = 150,
      springConstant = 0.08,
      damping = 0.4,
      avoidOverlap = 0.5
    ),
    stabilization = list(
      enabled = TRUE,
      iterations = 1000,
      updateInterval = 25
    )
  ) %>%

  # Interaction options
  visInteraction(
    hover = TRUE,
    hoverConnectedEdges = TRUE,
    selectConnectedEdges = TRUE,
    multiselect = TRUE,
    navigationButtons = TRUE,
    keyboard = TRUE,
    tooltipDelay = 100,
    zoomView = TRUE,
    dragView = TRUE
  ) %>%

  # Options
  visOptions(
    highlightNearest = list(
      enabled = TRUE,
      degree = 1,
      hover = TRUE,
      hideColor = "rgba(200,200,200,0.3)"
    ),
    nodesIdSelection = list(
      enabled = TRUE,
      main = "Select Antibiotic or Species"
    ),
    selectedBy = list(
      variable = "group",
      main = "Filter by Type"
    ),
    collapse = FALSE
  ) %>%

  # Layout
  visLayout(randomSeed = 42) %>%

  # Legend
  visLegend(
    addNodes = list(
      list(label = "Antibiotic", shape = "box", color = "#E74C3C", font = list(color = "white")),
      list(label = "Species (↑)", shape = "dot", color = "#27AE60"),
      list(label = "Species (↓)", shape = "dot", color = "#3498DB")
    ),
    addEdges = list(
      list(label = "Increased", color = "#27AE60", arrows = "to"),
      list(label = "Decreased", color = "#E74C3C", arrows = "to")
    ),
    useGroups = FALSE,
    position = "right",
    width = 0.15,
    stepY = 75
  ) %>%

  # Export button
  visExport(
    type = "png",
    name = "antibiotic_species_network",
    float = "left",
    label = "Export as PNG",
    style = "margin: 10px;"
  )

# =============================================================================
# 6. Save Outputs
# =============================================================================

# Check if pandoc is available for self-contained HTML
# Look in common locations if not in PATH
pandoc_paths <- c(
  Sys.which("pandoc"),
  "/home/david/miniforge3/bin/pandoc",
  "/usr/bin/pandoc",
  "/usr/local/bin/pandoc"
)
pandoc_found <- pandoc_paths[file.exists(pandoc_paths) & pandoc_paths != ""]
if (length(pandoc_found) > 0) {
  Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc_found[1]))
  pandoc_available <- TRUE
  cat("Using pandoc at:", pandoc_found[1], "\n")
} else {
  pandoc_available <- FALSE
}

# Save interactive HTML
html_file <- file.path(network_dir, "antibiotic_species_network.html")
saveWidget(network, html_file, selfcontained = pandoc_available)
cat("Saved interactive network:", html_file, "\n")
if (!pandoc_available) {
  cat("  (Note: pandoc not found - HTML requires lib/ folder to work)\n")
}

# =============================================================================
# 7. Create Summary Statistics
# =============================================================================

cat("\n=== Network Statistics ===\n\n")

# Create igraph object for statistics
g <- graph_from_data_frame(
  d = edges %>% select(from, to),
  directed = TRUE,
  vertices = nodes %>% select(id, group)
)

stats <- list(
  data_source = data_source,
  n_antibiotics = length(antibiotics),
  n_species = length(species_list),
  n_edges = nrow(edges),
  n_increased = sum(associations_top$direction == "increased"),
  n_decreased = sum(associations_top$direction == "decreased"),
  avg_degree_species = mean(degree(g, v = V(g)[V(g)$group == "species"])),
  avg_degree_antibiotics = mean(degree(g, v = V(g)[V(g)$group == "antibiotic"])),
  density = edge_density(g),
  config = config
)

# Print stats
cat("Data source:", stats$data_source, "\n")
cat("Antibiotics:", stats$n_antibiotics, "\n")
cat("Species:", stats$n_species, "\n")
cat("Edges:", stats$n_edges, "\n")
cat("  - Increased abundance:", stats$n_increased, "\n")
cat("  - Decreased abundance:", stats$n_decreased, "\n")
cat("Avg connections per antibiotic:", round(stats$avg_degree_antibiotics, 1), "\n")
cat("Avg connections per species:", round(stats$avg_degree_species, 2), "\n")
cat("Network density:", round(stats$density, 4), "\n")

# Save stats
saveRDS(stats, file.path(network_dir, "network_statistics.rds"))

# =============================================================================
# 8. Create Additional Views
# =============================================================================

cat("\n=== Creating Additional Network Views ===\n\n")

# --- View 2: Hierarchical layout (antibiotics on left, species on right) ---

hier_title <- ifelse(data_source == "aldex2_no_covariates",
                     "Antibiotic-Microbiome Network - Hierarchical View (PRELIMINARY - No Covariate Adjustment)",
                     "Antibiotic-Microbiome Network (Hierarchical View)")

network_hierarchical <- visNetwork(nodes, edges,
                                    main = hier_title,
                                    width = "100%", height = "900px") %>%
  visHierarchicalLayout(
    direction = "LR",
    sortMethod = "directed",
    levelSeparation = 300,
    nodeSpacing = 50
  ) %>%
  visInteraction(hover = TRUE, tooltipDelay = 100) %>%
  visOptions(highlightNearest = TRUE) %>%
  visExport(type = "png", name = "network_hierarchical")

html_hier <- file.path(network_dir, "antibiotic_species_network_hierarchical.html")
saveWidget(network_hierarchical, html_hier, selfcontained = pandoc_available)
cat("Saved hierarchical view:", html_hier, "\n")

# --- View 3: Per-antibiotic subnetworks ---

cat("\nGenerating per-antibiotic networks...\n")

for (abx in antibiotics) {
  abx_edges <- edges %>% filter(from == abx)
  abx_species <- unique(abx_edges$to)
  abx_nodes <- nodes %>% filter(id == abx | id %in% abx_species)

  if (nrow(abx_edges) == 0) next

  abx_title <- ifelse(data_source == "aldex2_no_covariates",
                      paste0(abx, " - Associated Species (PRELIMINARY)"),
                      paste0(abx, " - Associated Species"))

  abx_network <- visNetwork(abx_nodes, abx_edges,
                             main = abx_title,
                             width = "100%", height = "700px") %>%
    visPhysics(solver = "forceAtlas2Based",
               forceAtlas2Based = list(gravitationalConstant = -50)) %>%
    visInteraction(hover = TRUE) %>%
    visOptions(highlightNearest = TRUE) %>%
    visExport(type = "png", name = paste0("network_", gsub("/", "_", abx)))

  abx_file <- file.path(network_dir, paste0("network_", gsub("/", "_", abx), ".html"))
  saveWidget(abx_network, abx_file, selfcontained = pandoc_available)
}

cat("  Saved", length(antibiotics), "individual antibiotic networks\n")

# =============================================================================
# 9. Summary Table Export
# =============================================================================

# Export the filtered associations for reference
write_csv(associations_top, file.path(network_dir, "network_associations.csv"))

# Summary per antibiotic
abx_summary <- associations_top %>%
  group_by(antibiotic) %>%
  summarise(
    n_species = n(),
    n_increased = sum(direction == "increased"),
    n_decreased = sum(direction == "decreased"),
    top_increased = paste(head(species[direction == "increased"], 3), collapse = "; "),
    top_decreased = paste(head(species[direction == "decreased"], 3), collapse = "; "),
    mean_effect_magnitude = mean(abs(mean_effect)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_species))

write_csv(abx_summary, file.path(network_dir, "network_summary_by_antibiotic.csv"))

cat("\n=== Complete ===\n")
cat("\nOutput files in:", network_dir, "\n")
cat("  - antibiotic_species_network.html (main interactive network)\n")
cat("  - antibiotic_species_network_hierarchical.html (left-to-right layout)\n")
cat("  - network_<antibiotic>.html (individual antibiotic views)\n")
cat("  - network_associations.csv (data used for network)\n")
cat("  - network_summary_by_antibiotic.csv (summary table)\n")
cat("  - network_statistics.rds (network metrics)\n")
cat("\nOpen the HTML files in a browser to interact with the networks!\n")
