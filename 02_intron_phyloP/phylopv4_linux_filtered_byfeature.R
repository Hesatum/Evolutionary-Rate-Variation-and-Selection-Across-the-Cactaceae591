# ===============================================================================
#   PHYLOP POR FEATURE (elemento completo) - metodo nativo/recomendado do rphast
#
#   Em vez de rodar phyloP() sitio-a-sitio (basewise, o que gera 1 p-valor por
#   coluna do alinhamento e obriga a lidar com autocorrelacao espacial na hora
#   de corrigir por multiplos testes), usamos o argumento "features" do proprio
#   phyloP() (ver rphast-manual.pdf, pag. 80-81) para testar o intron INTEIRO
#   como um unico elemento por clado. Isso e o que a documentacao do pacote
#   descreve como o uso pretendido para testar elementos (genes, exons, UTRs,
#   ou aqui, introns completos) em vez de posicoes individuais.
#
#   Reusa exatamente a mesma calibracao do script principal (mesmo QC, mesma
#   arvore, mesmos clados, mesmo modelo neutro ajustado uma vez sobre o MSA
#   concatenado) - a UNICA diferenca e a chamada final ao phyloP().
# ===============================================================================

library(ape)
library(rphast)
library(dplyr)
library(purrr)
library(tools)

base_dir   <- "/home/user/Desktop/projetos/phast_rerun"
tree_file  <- file.path(base_dir, "filtered_tree.newick")
intron_dir <- file.path(base_dir, "introns/")
csv_dir    <- file.path(base_dir, "csv_filtered/")

MIN_INFORMATIVE_SITES <- 50
MIN_ALIGNMENT_LENGTH <- 100
DIVERSITY_THRESHOLD <- 0.05

set.seed(1)

cat("=== PHYLOP POR FEATURE (elemento completo, metodo nativo) ===\n")

tree <- read.tree(file = tree_file)
all_intron_files <- list.files(intron_dir, pattern = "\\.fasta$", full.names = TRUE)
cat("Arquivos encontrados:", length(all_intron_files), "\n")

all_msas <- purrr::map(all_intron_files, read.msa)
names(all_msas) <- basename(all_intron_files)

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
      file = filename, length = ncol.msa(msa), n_species = nrow.msa(msa),
      gap_proportion = gap_prop, missing_proportion = missing_prop,
      informative_sites = informative_sites, diversity = diversity,
      quality_passed = informative_sites >= MIN_INFORMATIVE_SITES &
        ncol.msa(msa) >= MIN_ALIGNMENT_LENGTH & diversity >= DIVERSITY_THRESHOLD
    )
  })
}

quality_report <- quality_check(all_msas)
good_indices <- which(quality_report$quality_passed)
cat("Alinhamentos aprovados:", length(good_indices), "de", nrow(quality_report), "\n")

selected_files <- all_intron_files[good_indices]
selected_msas <- all_msas[good_indices]

tree_taxa <- tree$tip.label
msa_taxa_list <- purrr::map(selected_msas, ~ .x$names)
common_taxa <- Reduce(intersect, msa_taxa_list)
common_taxa <- intersect(common_taxa, tree_taxa)
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

# ---- MESMA selecao de MSAs para o modelo neutro que o script principal ----
diversities <- quality_report$diversity[good_indices]
lengths <- quality_report$length[good_indices]
high_div_indices <- order(diversities, decreasing = TRUE)[1:min(8, length(diversities))]
varied_length_indices <- order(lengths, decreasing = TRUE)[1:min(5, length(lengths))]
neutral_indices <- unique(c(high_div_indices, varied_length_indices))
neutral_msas <- selected_msas_filtered[neutral_indices]
concatenated_msa <- concat.msa(neutral_msas)
cat("MSAs para modelo neutro:", length(neutral_msas), "| comprimento:", ncol.msa(concatenated_msa), "bp\n")

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
        phyloFit(msa = msa, tree = tree_string, subst.mod = model,
                 precision = config$precision, ninf.sites = config$min_informative,
                 init.backgd.from.data = TRUE, quiet = TRUE)
      }, error = function(e) NULL)
      if (!is.null(result) && !is.null(result$likelihood) && is.finite(result$likelihood)) {
        cat("  Sucesso com", model, "\n")
        return(list(model = result, method = model, config = config))
      }
    }
  }
  return(NULL)
}

cat("\n=== AJUSTANDO MODELO NEUTRO (identico ao script principal) ===\n")
neutral_result <- fit_robust_model(concatenated_msa, tree)
if (is.null(neutral_result)) stop("Erro: nao foi possivel ajustar o modelo neutro")
neutral_model <- neutral_result$model
cat("Modelo selecionado:", neutral_result$method, "| logLik:", neutral_model$likelihood, "\n")

# ===============================================================================
# ANALISE POR FEATURE: 1 p-valor por (intron, clado), via features= do phyloP()
# ===============================================================================
cat("\n=== EXECUTANDO PHYLOP POR FEATURE (intron inteiro como 1 elemento) ===\n")

all_feature_results <- list()
failed_analyses <- character(0)

for (i in seq_along(selected_files)) {
  filename <- tools::file_path_sans_ext(basename(selected_files[i]))
  msa <- selected_msas_filtered[[i]]

  if (i %% 10 == 0) cat(sprintf("Progresso: %d/%d\n", i, length(selected_files)))

  whole_feat <- feat(seqname = "MSA", start = 1, end = ncol.msa(msa),
                      feature = filename)

  for (clade_name in names(branch_mapping)) {
    node_name <- branch_mapping[[clade_name]]
    res <- tryCatch({
      phyloP(mod = neutral_model, msa = msa, method = "LRT", mode = "CONACC",
             branches = node_name, features = whole_feat)
    }, error = function(e) NULL)

    if (!is.null(res) && nrow(res) > 0) {
      res$og <- filename
      res$clade <- clade_name
      res$node_name <- node_name
      res$n_sites <- ncol.msa(msa)
      all_feature_results[[paste(filename, clade_name, sep = "_")]] <- res
    } else {
      failed_analyses <- c(failed_analyses, paste(filename, clade_name, sep = "_"))
    }
  }
}

feature_df <- dplyr::bind_rows(all_feature_results)
cat("\nLinhas obtidas:", nrow(feature_df), "| falhas:", length(failed_analyses), "\n")
cat("Colunas:", paste(names(feature_df), collapse=", "), "\n")

write.csv(feature_df, file.path(csv_dir, "10_phylop_por_feature_intron_completo.csv"), row.names = FALSE)
if (length(failed_analyses) > 0) {
  writeLines(failed_analyses, file.path(csv_dir, "10_failed_feature_analyses.txt"))
}

cat("\n=== FEITO ===\n")
cat("Salvo em:", file.path(csv_dir, "10_phylop_por_feature_intron_completo.csv"), "\n")
