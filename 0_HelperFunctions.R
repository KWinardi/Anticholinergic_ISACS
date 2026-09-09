# Helper Functions
# 25/09/2025 
# Kevin Winardi

compute_metrics_from_map <- function(graph, drug_targets, pathway_proteins,
                                     closeness_all    = NULL,
                                     include_centre   = FALSE,
                                     verbose          = TRUE) {
  library(pbapply)
  
  drug_targets     <- intersect(drug_targets,     merged.PDI.df2$Protein_A_SYMBOL)
  pathway_proteins <- intersect(pathway_proteins, V(graph)$name)
  
  if (include_centre && is.null(closeness_all)) {
    if (verbose) message("Computing closeness centrality...")
    closeness_all        <- closeness(graph, normalized = TRUE)
    names(closeness_all) <- V(graph)$name
  }
  
  if (verbose) message("Preparing drug-target map...")
  drug_target_map <- lapply(drug_targets, function(drug) {
    if (drug %in% merged.PDI.df2$Protein_A_SYMBOL)
      merged.PDI.df2 %>% dplyr::filter(Protein_A_SYMBOL == drug) %>% pull(Protein_B_SYMBOL)
    else character(0)
  })
  names(drug_target_map) <- drug_targets
  
  if (verbose) message("Computing observed distance metrics...")
  
  results <- pbapply::pblapply(seq_along(drug_target_map), function(i) {
    drug <- names(drug_target_map)[i]
    A    <- intersect(drug_target_map[[i]], V(graph)$name)
    B    <- intersect(pathway_proteins,     V(graph)$name)
    
    if (length(A) == 0 || length(B) == 0) return(NULL)
    
    dAB      <- distances(graph, v = A, to = B)
    finiteAB <- is.finite(dAB)
    
    min_a  <- apply(dAB, 1, function(row) if (any(is.finite(row))) min(row[is.finite(row)]) else NA_real_)
    min_b  <- apply(dAB, 2, function(col) if (any(is.finite(col))) min(col[is.finite(col)]) else NA_real_)
    d_c_b  <- mean(min_b, na.rm = TRUE)
    d_c_ab <- mean(c(min_a, min_b), na.rm = TRUE)
    d_s    <- mean(dAB[finiteAB], na.rm = TRUE)
    
    d_k_b_sum <- 0
    for (j in seq_len(ncol(dAB))) {
      dj_fin <- dAB[, j][is.finite(dAB[, j])]
      if (length(dj_fin) > 0)
        d_k_b_sum <- d_k_b_sum + log(mean(exp(-(dj_fin + 1))))
    }
    d_k_b <- -d_k_b_sum / length(B)
    
    dBA       <- distances(graph, v = B, to = A)
    d_k_a_sum <- 0
    for (j in seq_len(ncol(dBA))) {
      dj_fin <- dBA[, j][is.finite(dBA[, j])]
      if (length(dj_fin) > 0)
        d_k_a_sum <- d_k_a_sum + log(mean(exp(-(dj_fin + 1))))
    }
    d_k_ab <- -(d_k_b_sum + d_k_a_sum) / (length(A) + length(B))
    
    # Centre proximity - only computed when include_centre = TRUE
    d_cc <- if (include_centre) {
      sub_cc       <- closeness_all[A]
      max_cc       <- max(sub_cc, na.rm = TRUE)
      center_nodes <- names(sub_cc)[sub_cc == max_cc]
      cdists       <- numeric()
      for (cnode in center_nodes) {
        if (cnode %in% V(graph)$name) {
          ds     <- distances(graph, v = cnode, to = B)
          cdists <- c(cdists, ds[is.finite(ds)])
        }
      }
      if (length(cdists) > 0) mean(cdists) else NA_real_
    } else {
      NULL
    }
    
    out <- data.frame(
      Drug              = drug,
      DrugTarget_degree = length(A),
      Closest_b         = d_c_b,
      Closest_ab        = d_c_ab,
      Shortest          = d_s,
      Kernel_b          = d_k_b,
      Kernel_ab         = d_k_ab,
      stringsAsFactors  = FALSE
    )
    if (include_centre) out$Centre <- d_cc
    out
  })
  
  do.call(rbind, results)
}


compute_distance_significance <- function(graph,
                                          drug_targets,
                                          pathway_proteins,
                                          distance_param  = NULL,
                                          distance_value  = NULL,
                                          closeness_all   = NULL,
                                          include_centre  = FALSE,
                                          n_perm          = 1000,
                                          randomize_drugs = TRUE,
                                          verbose         = TRUE,
                                          seed            = NULL,
                                          n_cores         = parallel::detectCores() - 2) {
  start_time <- Sys.time()
  
  if (!is.null(seed)) {
    set.seed(seed)
    if (verbose) message("Random seed set to: ", seed)
  }
  
  sample_degree_preserving <- function(graph, target_nodes) {
    degs      <- degree(graph, v = target_nodes)
    all_nodes <- V(graph)$name
    all_degs  <- degree(graph, v = all_nodes)
    sampled   <- character(length(degs))
    for (i in seq_along(degs)) {
      d          <- degs[i]
      candidates <- setdiff(all_nodes[all_degs == d],          sampled[seq_len(i - 1)])
      if (length(candidates) == 0)
        candidates <- setdiff(all_nodes[abs(all_degs - d) <= 1], sampled[seq_len(i - 1)])
      sampled[i] <- if (length(candidates) > 0) sample(candidates, 1) else NA_character_
    }
    sampled[!is.na(sampled)]
  }
  
  compute_metrics_from_map.v2 <- function(graph, drug_target_map, pathway_proteins,
                                          include_centre = FALSE,
                                          verbose        = FALSE) {
    library(pbapply)
    
    results <- pbapply::pblapply(seq_along(drug_target_map), function(i) {
      drug <- names(drug_target_map)[i]
      A    <- intersect(drug_target_map[[i]], V(graph)$name)
      B    <- intersect(pathway_proteins,     V(graph)$name)
      
      if (length(A) == 0 || length(B) == 0) return(NULL)
      
      dAB      <- distances(graph, v = A, to = B)
      finiteAB <- is.finite(dAB)
      
      min_a  <- apply(dAB, 1, function(row) if (any(is.finite(row))) min(row[is.finite(row)]) else NA_real_)
      min_b  <- apply(dAB, 2, function(col) if (any(is.finite(col))) min(col[is.finite(col)]) else NA_real_)
      d_c_b  <- mean(min_b, na.rm = TRUE)
      d_c_ab <- mean(c(min_a, min_b), na.rm = TRUE)
      d_s    <- mean(dAB[finiteAB], na.rm = TRUE)
      
      d_k_b_sum <- 0
      for (j in seq_len(ncol(dAB))) {
        dj_fin <- dAB[, j][is.finite(dAB[, j])]
        if (length(dj_fin) > 0)
          d_k_b_sum <- d_k_b_sum + log(mean(exp(-(dj_fin + 1))))
      }
      d_k_b <- -d_k_b_sum / length(B)
      
      dBA       <- distances(graph, v = B, to = A)
      d_k_a_sum <- 0
      for (j in seq_len(ncol(dBA))) {
        dj_fin <- dBA[, j][is.finite(dBA[, j])]
        if (length(dj_fin) > 0)
          d_k_a_sum <- d_k_a_sum + log(mean(exp(-(dj_fin + 1))))
      }
      d_k_ab <- -(d_k_b_sum + d_k_a_sum) / (length(A) + length(B))
      
      # Centre proximity - only if requested (expensive; default FALSE in perms)
      d_cc <- if (include_centre) {
        sub_cc       <- closeness_all[A]
        max_cc       <- max(sub_cc, na.rm = TRUE)
        center_nodes <- names(sub_cc)[sub_cc == max_cc]
        cdists       <- numeric()
        for (cnode in center_nodes) {
          if (cnode %in% V(graph)$name) {
            ds     <- distances(graph, v = cnode, to = B)
            cdists <- c(cdists, ds[is.finite(ds)])
          }
        }
        if (length(cdists) > 0) mean(cdists) else NA_real_
      } else {
        NULL
      }
      
      out <- data.frame(
        Drug              = drug,
        DrugTarget_degree = length(A),
        Closest_b         = d_c_b,
        Closest_ab        = d_c_ab,
        Shortest          = d_s,
        Kernel_b          = d_k_b,
        Kernel_ab         = d_k_ab,
        stringsAsFactors  = FALSE
      )
      if (include_centre) out$Centre <- d_cc
      out
    })
    
    do.call(rbind, results)
  }
  
  if (include_centre && is.null(closeness_all)) {
    if (verbose) message("Computing closeness centrality...")
    closeness_all        <- closeness(graph, normalized = TRUE)
    names(closeness_all) <- V(graph)$name
  }
  
  drug_targets     <- intersect(drug_targets,     merged.PDI.df2$Protein_A_SYMBOL)
  if (verbose) message("Assessing ", length(drug_targets), " drugs...")
  pathway_proteins <- intersect(pathway_proteins, V(graph)$name)
  
  if (verbose) message("Preparing drug-target map...")
  drug_target_map <- lapply(drug_targets, function(drug) {
    if (drug %in% merged.PDI.df2$Protein_A_SYMBOL)
      merged.PDI.df2 %>% dplyr::filter(Protein_A_SYMBOL == drug) %>% pull(Protein_B_SYMBOL)
    else character(0)
  })
  names(drug_target_map) <- drug_targets
  
  if (is.null(distance_value)) {
    if (verbose) message("Computing observed distance metrics...")
    observed <- compute_metrics_from_map(graph, drug_targets, pathway_proteins,
                                         closeness_all  = closeness_all,
                                         include_centre = include_centre,
                                         verbose        = verbose)
  } else {
    if (verbose) message("Using pre-calculated distance metrics...")
    observed <- distance_value
  }
  
  if (is.null(distance_param)) {
    metric_candidates <- setdiff(names(observed), c("Drug", "DrugTarget_degree"))
    metric_classes    <- sapply(observed[metric_candidates], class)
    distance_param    <- metric_candidates[metric_classes %in% c("numeric", "integer")]
  }
  if (verbose) message("Distance parameters assessed: ", paste(distance_param, collapse = ", "))
  
  permute_once <- function(i) {
    rand_pathway <- sample_degree_preserving(graph, pathway_proteins)
    
    rand_drug_target_map <- if (randomize_drugs) {
      nm <- lapply(drug_target_map, function(targets) sample_degree_preserving(graph, targets))
      names(nm) <- names(drug_target_map)
      nm
    } else {
      drug_target_map
    }
    
    rand_result <- compute_metrics_from_map.v2(graph, rand_drug_target_map, rand_pathway,
                                               include_centre = include_centre)
    
    res <- lapply(distance_param, function(metric) {
      sapply(observed$Drug, function(drug) {
        row_idx <- which(rand_result$Drug == drug)
        if (length(row_idx) == 1) rand_result[[metric]][row_idx] else NA_real_
      })
    })
    names(res) <- distance_param
    res
  }
  
  # Parallel permutations
  if (verbose) message("Running ", n_perm, " permutations using ", n_cores, " cores...")
  
  export_vars <- c("graph", "drug_target_map", "pathway_proteins",
                   "randomize_drugs", "distance_param",
                   "observed", "sample_degree_preserving",
                   "compute_metrics_from_map.v2", "include_centre",
                   "merged.PDI.df2")
  if (include_centre) export_vars <- c(export_vars, "closeness_all")
  
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  clusterExport(cl, varlist = export_vars, envir = environment())
  clusterEvalQ(cl, library(igraph))
  
  perm_results <- pbapply::pblapply(1:n_perm, permute_once, cl = cl)
  
  #  Null distributions
  if (verbose) message("Calculating z-scores and p-values...")
  null_metrics <- lapply(distance_param, function(metric) {
    mat             <- do.call(cbind, lapply(perm_results, `[[`, metric))
    rownames(mat)   <- observed$Drug
    mat
  })
  names(null_metrics) <- distance_param
  
  # z-scores and p-values
  z_to_p <- function(z, one.sided = NULL) {
    if (is.null(one.sided)) {
      pval        <- pnorm(-abs(z)) * 2
      pval[pval > 1] <- 1
    } else if (one.sided == "-") {
      pval <- pnorm(z)
    } else {
      pval <- pnorm(-z)
    }
    pval
  }
  
  for (metric in distance_param) {
    obs_vals  <- observed[[metric]]
    null_mat  <- null_metrics[[metric]]
    zscores   <- mapply(function(obs, null_dist) {
      mu  <- mean(null_dist, na.rm = TRUE)
      sdv <- sd(null_dist,   na.rm = TRUE)
      if (is.na(sdv) || sdv == 0) NA_real_ else (obs - mu) / sdv
    }, obs_vals, lapply(seq_len(nrow(null_mat)), function(i) null_mat[i, ]))
    
    observed[[paste0(metric, "_pval")]]   <- z_to_p(zscores)
    observed[[paste0(metric, "_zscore")]] <- zscores
  }
  base_cols   <- c("Drug", "DrugTarget_degree")
  metric_cols <- unlist(lapply(distance_param, function(m)
    c(m, paste0(m, "_pval"), paste0(m, "_zscore"))))
  observed    <- observed[, c(base_cols, metric_cols)]
  end_time   <- Sys.time()
  total_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
  if (verbose) {
    message("Total run time: ",
            floor(total_secs / 3600), " hr ",
            floor((total_secs %% 3600) / 60), " min ",
            round(total_secs %% 60, 1), " sec")
  }
  
  observed
}

compute_shortest_paths <- function(graph, drug_targets, pathway_proteins) {
  drug_targets <- intersect(drug_targets, merged.PDI.df2$Protein_A_SYMBOL)
  pathway_proteins <- intersect(pathway_proteins, V(graph)$name)
  n_cores <- max(1, parallel::detectCores() - 2)
  cl <- makeCluster(n_cores)
  # Export necessary variables and libraries to cluster
  clusterExport(cl, varlist = c("graph", "drug_targets", "pathway_proteins", "merged.PDI.df2"), envir = environment())
  clusterEvalQ(cl, {
    library(igraph)
    library(dplyr)
  })
  
  process_drug <- function(d) {
    result_list <- list()
    
    d.targets <- merged.PDI.df2 %>%
      filter(Protein_A_SYMBOL == d) %>%
      pull(Protein_B_SYMBOL) %>%
      unique()
    d.targets <- intersect(d.targets, V(graph)$name)
    
    if (length(d.targets) == 0) return(NULL)
    
    for (to_node in pathway_proteins) {
      for (from_node in d.targets) {
        all_paths <- all_shortest_paths(graph, from = from_node, to = to_node)$res
        
        if (length(all_paths) > 0) {
          for (p in all_paths) {
            path_nodes <- V(graph)[p]$name
            path_str <- paste(path_nodes, collapse = " --> ")
            path_length <- length(path_nodes) - 1
            
            result_list[[length(result_list) + 1]] <- data.frame(
              drug = d,
              from = from_node,
              to = to_node,
              path = path_str,
              distance = path_length,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    if (length(result_list) > 0) {
      return(do.call(rbind, result_list))
    } else {
      return(NULL)
    }
  }
  
  results_list <- pbapply::pblapply(drug_targets, process_drug, cl = cl)
  
  # Stop the cluster
  stopCluster(cl)
  
  # Combine results
  final_results <- do.call(rbind, results_list)
  
  return(final_results)
}
