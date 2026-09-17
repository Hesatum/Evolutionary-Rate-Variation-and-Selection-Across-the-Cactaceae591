# ===============================================================================
#   ANALISE PHYLOP PARA INTRONS DE CEREUS
#
#   Pipeline completo de analise filogenomica usando PhyloP
#   Detecta selecao positiva/negativa em diferentes linhagens
#   (adaptado para rodar em Linux, caminhos locais desta maquina)
# ===============================================================================

library(ape)
library(rphast)
library(ggplot2)
library(dplyr)
library(purrr)
library(tools)
library(tidyr)
library(viridis)

# Configuracao dos diretorios
base_dir   <- "/home/user/Desktop/projetos/phast_rerun"
tree_file  <- file.path(base_dir, "filtered_tree.newick")
intron_dir <- file.path(base_dir, "introns/")
plot_dir   <- file.path(base_dir, "plots_filtered/")
csv_dir    <- file.path(base_dir, "csv_filtered/")

plot_individual_dir <- file.path(plot_dir, "perfis_individuais")
plot_comparison_dir <- file.path(plot_dir, "comparacoes_metodos")
csv_detailed_dir <- file.path(csv_dir, "resultados_detalhados")

dir.create(plot_individual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_comparison_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_detailed_dir, recursive = TRUE, showWarnings = FALSE)

# Parametros de controle de qualidade
MIN_INFORMATIVE_SITES <- 50
MIN_ALIGNMENT_LENGTH <- 100
DIVERSITY_THRESHOLD <- 0.05

cat("=== INICIANDO ANALISE PHYLOP ===\n")
cat("Parametros:\n")
cat("- Minimo de sitios informativos:", MIN_INFORMATIVE_SITES, "\n")
cat("- Comprimento minimo:", MIN_ALIGNMENT_LENGTH, "\n")
cat("- Threshold de diversidade:", DIVERSITY_THRESHOLD, "\n\n")

# Carrega arvore e lista de arquivos
tree <- read.tree(file = tree_file)
all_intron_files <- list.files(intron_dir, pattern = "\\.fasta$", full.names = TRUE)
cat("Arquivos encontrados:", length(all_intron_files), "\n")

# Le todos os alinhamentos
all_msas <- purrr::map(all_intron_files, read.msa)
names(all_msas) <- basename(all_intron_files)

# Funcao de controle de qualidade dos alinhamentos
quality_check <- function(msa_list) {
  map_dfr(seq_along(msa_list), function(i) {
    msa <- msa_list[[i]]
    filename <- names(msa_list)[i]

    seqs <- as.character(msa$seqs)
    mat <- do.call(rbind, strsplit(seqs, ""))

    gap_prop <- mean(mat == "-")
    missing_prop <- mean(mat == "N")

    informative_sites <- sum(apply(mat, 2, function(col) {
      non_gap <- col[col != "-" & col != "N"]
      length(unique(non_gap)) >= 2 & length(non_gap) >= 2
    }))

    diversity <- mean(apply(mat, 2, function(col) {
      non_gap <- col[col != "-" & col != "N"]
      length(unique(non_gap)) > 1
    }), na.rm = TRUE)

    data.frame(
      file = filename,
      length = ncol.msa(msa),
      n_species = nrow.msa(msa),
      gap_proportion = gap_prop,
      missing_proportion = missing_prop,
      informative_sites = informative_sites,
      diversity = diversity,
      quality_passed = informative_sites >= MIN_INFORMATIVE_SITES &
        ncol.msa(msa) >= MIN_ALIGNMENT_LENGTH &
        diversity >= DIVERSITY_THRESHOLD
    )
  })
}

quality_report <- quality_check(all_msas)
write.csv(quality_report, file.path(csv_dir, "01_quality_control.csv"), row.names = FALSE)

good_indices <- which(quality_report$quality_passed)
cat("Alinhamentos aprovados:", length(good_indices), "de", nrow(quality_report), "\n")

if (length(good_indices) == 0) {
  stop("Nenhum alinhamento passou no controle de qualidade.")
}

selected_files <- all_intron_files[good_indices]
selected_msas <- all_msas[good_indices]

# Identifica especies comuns
tree_taxa <- tree$tip.label
msa_taxa_list <- purrr::map(selected_msas, ~ .x$names)
common_taxa <- Reduce(intersect, msa_taxa_list)
common_taxa <- intersect(common_taxa, tree_taxa)

cat("Especies comuns a todos os dados:", length(common_taxa), "\n")

tree <- keep.tip(tree, common_taxa)

filter_msa_robust <- function(msa_obj, keep_taxa) {
  keep_indices <- match(keep_taxa, msa_obj$names)
  keep_indices <- keep_indices[!is.na(keep_indices)]
  if (length(keep_indices) == 0) return(NULL)
  filtered_msa <- msa_obj
  filtered_msa$seqs <- msa_obj$seqs[keep_indices]
  filtered_msa$names <- msa_obj$names[keep_indices]
  filtered_msa$nseq <- length(keep_indices)
  return(filtered_msa)
}

selected_msas_filtered <- purrr::map(selected_msas, ~ filter_msa_robust(.x, common_taxa))
selected_msas_filtered <- purrr::compact(selected_msas_filtered)

# Define os clados biologicos
clades <- list(
  A1 = c('Cereus_trigonodendron_S184A4', 'Cereus_bicolor_S102A10',
         'Cereus_jamacaru_S180V2', 'Cereushexagonus_S155A2'),
  A2 = c('Cereussp.Nov_S169A2', 'Cereusfernambucensis_S80F9',
         'Cereusfernambucensissericifer_S88F2'),
  B  = c('Cereus_spegazzini_S180A14', 'Cereus_vargasianus_S184A5',
         'Cereusspegazzinii_S77A31'),
  C  = c('Cereussaddianus_S103D4', 'Cereusphatnospermus_S149V10'),
  D  = c('Cereus_mirabella_S162A26', 'Cereusalbicaulis_S143F5'),
  E  = c('Cereus_mortensenii_S184A3', 'Cereusrepandus_S162VA22'),
  Outgroup = c('Cipocereuslaniflorus_S174A6', 'Cipocereusminensis_S153A1')
)

clade_colors <- c(
  A1 = "#E69F00",
  A2 = "#56B4E9",
  B = "#009E73",
  C = "#F0E442",
  D = "#0072B2",
  E = "#D55E00",
  Outgroup = "#999999"
)

clades_validated <- map(clades, ~ .x[.x %in% common_taxa])
clades_validated <- keep(clades_validated, ~ length(.x) >= 2)

cat("Clados validos:", paste(names(clades_validated), collapse = ", "), "\n")

tree$node.label <- character(tree$Nnode)
branch_mapping <- list()

for (clade_name in names(clades_validated)) {
  species_in_clade <- clades_validated[[clade_name]]
  mrca_node <- getMRCA(tree, species_in_clade)
  if (!is.null(mrca_node) && mrca_node > Ntip(tree)) {
    node_index <- mrca_node - Ntip(tree)
    node_label <- paste0(clade_name, "_ancestor")
    tree$node.label[node_index] <- node_label
    branch_mapping[[clade_name]] <- node_label
  }
}

clade_info <- data.frame(
  clade = names(clades_validated),
  n_species = sapply(clades_validated, length),
  ancestor_node = sapply(names(clades_validated), function(x) branch_mapping[[x]]),
  color = clade_colors[names(clades_validated)]
)
write.csv(clade_info, file.path(csv_dir, "02_clade_info.csv"), row.names = FALSE)

# Seleciona MSAs para o modelo neutro
diversities <- quality_report$diversity[good_indices]
lengths <- quality_report$length[good_indices]

high_div_indices <- order(diversities, decreasing = TRUE)[1:min(8, length(diversities))]
varied_length_indices <- order(lengths, decreasing = TRUE)[1:min(5, length(lengths))]
neutral_indices <- unique(c(high_div_indices, varied_length_indices))

neutral_msas <- selected_msas_filtered[neutral_indices]
concatenated_msa <- concat.msa(neutral_msas)

cat("MSAs para modelo neutro:", length(neutral_msas), "\n")
cat("Comprimento concatenado:", ncol.msa(concatenated_msa), "bp\n")

neutral_info <- data.frame(
  intron = basename(selected_files[neutral_indices]),
  length = sapply(neutral_msas, ncol.msa),
  diversity = diversities[neutral_indices]
)
write.csv(neutral_info, file.path(csv_dir, "03_neutral_model_data.csv"), row.names = FALSE)

# Ajusta modelo neutro com multiplas tentativas
fit_robust_model <- function(msa, tree_obj, attempt_models = c("HKY85", "REV", "JC69")) {
  tree_string <- write.tree(tree_obj)

  configs <- list(
    list(precision = "HIGH", min_informative = 50),
    list(precision = "MED", min_informative = 25),
    list(precision = "LOW", min_informative = 10)
  )

  for (model in attempt_models) {
    for (config in configs) {
      cat("Testando", model, "com precisao", config$precision, "...\n")

      result <- tryCatch({
        phyloFit(
          msa = msa,
          tree = tree_string,
          subst.mod = model,
          precision = config$precision,
          ninf.sites = config$min_informative,
          init.backgd.from.data = TRUE,
          quiet = TRUE
        )
      }, error = function(e) NULL)

      if (!is.null(result) && !is.null(result$likelihood) && is.finite(result$likelihood)) {
        cat("  Sucesso com", model, "\n")
        return(list(model = result, method = model, config = config))
      }
    }
  }
  return(NULL)
}

cat("\n=== AJUSTANDO MODELO NEUTRO ===\n")
neutral_result <- fit_robust_model(concatenated_msa, tree)

if (is.null(neutral_result)) {
  stop("Erro: nao foi possivel ajustar o modelo neutro")
}

neutral_model <- neutral_result$model
cat("Modelo selecionado:", neutral_result$method, "\n")
cat("Log-likelihood:", neutral_model$likelihood, "\n")

model_params <- data.frame(
  parameter = c("model", "log_likelihood", "precision", "min_informative"),
  value = c(neutral_result$method, neutral_model$likelihood,
            neutral_result$config$precision, neutral_result$config$min_informative)
)
write.csv(model_params, file.path(csv_dir, "04_neutral_model_parameters.csv"), row.names = FALSE)

# Analise PhyloP com multiplos metodos
analyze_intron_robust <- function(model, msa, filename, branch_mapping, methods = c("LRT", "SCORE")) {
  all_results <- list()

  for (method in methods) {
    for (clade_name in names(branch_mapping)) {
      node_name <- branch_mapping[[clade_name]]

      result <- tryCatch({
        phyloP(
          mod = model,
          msa = msa,
          method = method,
          mode = "CONACC",
          branches = node_name
        )
      }, error = function(e) NULL)

      if (!is.null(result)) {
        df <- as.data.frame(result)
        if (nrow(df) > 0) {
          df$clade <- clade_name
          df$gene <- filename
          df$method <- method
          df$mode <- "CONACC"
          df$node_name <- node_name
          result_key <- paste(filename, clade_name, method, sep = "_")
          all_results[[result_key]] <- df
        }
      }
    }
  }

  return(all_results)
}

# Analise principal
cat("\n=== EXECUTANDO ANALISE PHYLOP ===\n")
all_results <- list()
failed_analyses <- character(0)

for (i in seq_along(selected_files)) {
  filename <- tools::file_path_sans_ext(basename(selected_files[i]))
  msa <- selected_msas_filtered[[i]]

  if (i %% 5 == 0) {
    cat(sprintf("Progresso: %d/%d (%.1f%%)\n", i, length(selected_files),
                100 * i / length(selected_files)))
  }

  intron_results <- analyze_intron_robust(neutral_model, msa, filename, branch_mapping)

  if (length(intron_results) == 0) {
    failed_analyses <- c(failed_analyses, filename)
    next
  }

  # Gera grafico individual (metodo LRT)
  lrt_results <- intron_results[grepl("_LRT$", names(intron_results))]

  if (length(lrt_results) > 0) {
    combined_df <- dplyr::bind_rows(lrt_results)

    n_total <- nrow(combined_df)
    n_conserved <- sum(combined_df$score > 1, na.rm = TRUE)
    n_accelerated <- sum(combined_df$score < -1, na.rm = TRUE)

    p <- ggplot2::ggplot(combined_df, aes(x = coord, y = score, color = clade)) +
      ggplot2::geom_point(size = 1.5, alpha = 0.6, shape = 16) +
      ggplot2::scale_color_manual(values = clade_colors, name = "Clade") +
      ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "gray30", linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.4) +
      ggplot2::labs(
        title = paste("Evolutionary Profile:", filename),
        subtitle = sprintf("Sites: %d | Conserved: %d (%.1f%%) | Accelerated: %d (%.1f%%)",
                            n_total, n_conserved, 100*n_conserved/n_total,
                            n_accelerated, 100*n_accelerated/n_total),
        x = "Genomic Position (bp)",
        y = "PhyloP Score (LRT)",
        caption = "Dashed lines: +-1 threshold | Positive scores: conservation | Negative scores: acceleration"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 14),
        plot.subtitle = ggplot2::element_text(hjust = 0, size = 10, color = "gray20"),
        plot.caption = ggplot2::element_text(size = 8, color = "gray40", hjust = 0),
        legend.title = ggplot2::element_text(face = "bold", size = 11),
        legend.text = ggplot2::element_text(size = 10),
        legend.position = "right",
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
        panel.border = ggplot2::element_rect(fill = NA, color = "gray70", linewidth = 0.5),
        axis.title = ggplot2::element_text(face = "bold", size = 11),
        axis.text = ggplot2::element_text(size = 10)
      )

    ggplot2::ggsave(
      file.path(plot_individual_dir, paste0(filename, "_profile_LRT.png")),
      p, width = 12, height = 6, dpi = 300, bg = "white"
    )

    write.csv(combined_df,
              file.path(csv_detailed_dir, paste0(filename, "_LRT_scores.csv")),
              row.names = FALSE)
  }

  for (result_name in names(intron_results)) {
    all_results[[result_name]] <- intron_results[[result_name]]
  }
}

cat("\nAnalise concluida!\n")
cat("Sucessos:", length(selected_files) - length(failed_analyses), "\n")
cat("Falhas:", length(failed_analyses), "\n")

if (length(failed_analyses) > 0) {
  writeLines(failed_analyses, file.path(csv_dir, "05_failed_analyses.txt"))
}

# Consolida todos os resultados
if (length(all_results) > 0) {
  final_df <- dplyr::bind_rows(all_results, .id = "result_id")
  write.csv(final_df, file.path(csv_dir, "06_all_phylop_results.csv"), row.names = FALSE)

  method_stats <- final_df %>%
    dplyr::group_by(method, clade) %>%
    dplyr::summarise(
      n_sites = n(),
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      sd_score = sd(score, na.rm = TRUE),
      conserved_sites = sum(score > 1, na.rm = TRUE),
      accelerated_sites = sum(score < -1, na.rm = TRUE),
      neutral_sites = sum(abs(score) <= 1, na.rm = TRUE),
      pct_conserved = 100 * conserved_sites / n_sites,
      pct_accelerated = 100 * accelerated_sites / n_sites,
      .groups = 'drop'
    )

  write.csv(method_stats, file.path(csv_dir, "07_statistics_by_method_clade.csv"), row.names = FALSE)

  gene_stats <- final_df %>%
    dplyr::filter(method == "LRT") %>%
    dplyr::group_by(gene, clade) %>%
    dplyr::summarise(
      n_sites = n(),
      mean_score = mean(score, na.rm = TRUE),
      conserved_sites = sum(score > 1, na.rm = TRUE),
      accelerated_sites = sum(score < -1, na.rm = TRUE),
      pct_conserved = 100 * conserved_sites / n_sites,
      pct_accelerated = 100 * accelerated_sites / n_sites,
      .groups = 'drop'
    )

  write.csv(gene_stats, file.path(csv_dir, "08_statistics_by_gene.csv"), row.names = FALSE)

  cat("\n=== ESTATISTICAS FINAIS ===\n")
  cat("Total de sitios analisados:", nrow(final_df), "\n")
  cat("Metodos:", paste(unique(final_df$method), collapse = ", "), "\n")
  cat("Clados:", paste(unique(final_df$clade), collapse = ", "), "\n")

  significant_sites <- sum(abs(final_df$score) > 1, na.rm = TRUE)
  conserved_sites <- sum(final_df$score > 1, na.rm = TRUE)
  accelerated_sites <- sum(final_df$score < -1, na.rm = TRUE)

  cat(sprintf("\nSitios significativos: %d (%.2f%%)\n",
              significant_sites, 100 * significant_sites / nrow(final_df)))
  cat(sprintf("  - Conservados: %d (%.2f%%)\n",
              conserved_sites, 100 * conserved_sites / nrow(final_df)))
  cat(sprintf("  - Acelerados: %d (%.2f%%)\n",
              accelerated_sites, 100 * accelerated_sites / nrow(final_df)))

  summary_df <- data.frame(
    metric = c("total_sites", "significant_sites", "conserved_sites", "accelerated_sites",
               "neutral_sites", "pct_significant", "pct_conserved", "pct_accelerated"),
    value = c(nrow(final_df), significant_sites, conserved_sites, accelerated_sites,
              nrow(final_df) - significant_sites,
              100 * significant_sites / nrow(final_df),
              100 * conserved_sites / nrow(final_df),
              100 * accelerated_sites / nrow(final_df))
  )
  write.csv(summary_df, file.path(csv_dir, "09_summary_statistics.csv"), row.names = FALSE)
}

cat("\n=== ANALISE COMPLETA ===\n")
cat("CSVs salvos em:", csv_dir, "\n")
cat("Graficos salvos em:", plot_dir, "\n")
