

#' @title Indexes of PGs embedded in each others
#' @description Returns indexes of PGs that are embedded in others
#'
#' @param adj An adjacency matrix between precursors/peptides and PGs
#'
#' @return A vector of indices
#' @export
#'
#' @examples
#' data(subbouyssie)
#' get_indexes_embedded_prots(subbouyssie$adj)
#'
get_indexes_embedded_prots <- function(adj) {
  similarity.mat = t(adj) %*% adj
  mat.size = ncol(similarity.mat)
  idx.prot = 1
  idx.prot.rm = c()
  for (idx.prot in 1:mat.size) {
    if (any(similarity.mat[idx.prot, idx.prot] ==
            similarity.mat[idx.prot, -c(idx.prot, idx.prot.rm)])) {
      idx.prot.rm = c(idx.prot.rm, idx.prot)
    }
  }
  return(idx.prot.rm)
}


#' @title Remove PGs by index and merge
#' @description Remove PG by index and merge transcripts (if transcriptomic 
#' information is available) of PG included in one another (under condition 
#' that they have peptide). Then it removes transcripts without PG. Do not 
#' remove peptides that are left without PG.
#'
#' @param l_pep_rna A list representing dataset, formatted as in 
#' pipeline_llkimpute function
#' @param pg_idx Vector of indices
#'
#' @return A list representing dataset.
#' @export
#'
#' @examples
#' NULL
#'
rm_pg_from_idx_merge_pg <- function(l_pep_rna, pg_idx) {
  if (!(length(pg_idx) == 0) & !is.null(pg_idx)) {
    adj2keep = l_pep_rna$adj[, -pg_idx, drop = FALSE]
    adj2rm = l_pep_rna$adj[, pg_idx, drop = FALSE]
    l_pep_rna$adj = adj2keep
    if (!is.null(l_pep_rna$adj_rna_pg)) {
      n_pep_per_pg2rm = colSums(adj2rm)
      pg_sim_mat = ((t(adj2rm) %*% adj2keep) == n_pep_per_pg2rm) &
        (n_pep_per_pg2rm != 0)
      idx.pg2merge = which(pg_sim_mat, arr.ind = TRUE)
      adj2keep_rna = l_pep_rna$adj_rna_pg[, -pg_idx, drop = FALSE]
      adj2rm_rna = l_pep_rna$adj_rna_pg[, pg_idx, drop = FALSE] 
      if (length(idx.pg2merge) != 0) {
        for (i in 1:nrow(idx.pg2merge)) {
          adj2keep_rna[, idx.pg2merge[i, 2]] <-
            adj2keep_rna[, idx.pg2merge[i, 2]] | 
            adj2rm_rna[, idx.pg2merge[i, 1]]
        }
      }
      l_pep_rna$adj_rna_pg = adj2keep_rna
      i_rna_rm = which(rowSums(adj2keep_rna) == 0)
      if (length(i_rna_rm) != 0) {
        l_pep_rna$rnas_ab = l_pep_rna$rnas_ab[,-i_rna_rm]
        l_pep_rna$adj_rna_pg = l_pep_rna$adj_rna_pg[-i_rna_rm, ]
      }
    } else {
      l_pep_rna$adj = adj2keep
      if (!is.null(l_pep_rna$mask_prot_diff)) {
        l_pep_rna$mask_prot_diff = l_pep_rna$mask_prot_diff[-pg_idx]
      }
    }
    
  }
  return(l_pep_rna)
}




# TODO: Rename this function in code. 
# TODO: Impfunc shoud not be passed as parameter, by we should directly call 
# TODO: "estimate params & impute.
# TODO: Some parameters are not used anymore in pipeline_ll_imp file, need to 
# TODO: remove them
#' @title Impute each PG.
#' @description Imputes each PG separately and return the results for each PG. 
#'
#' @param data.pep.rna.crop A list representing dataset
# #' @param impfunc Imputation function
#' @param psi Inverse scale parameter for IW prior of peptides abundances
#' @param pep_ab_or In case we impute a dataset with pseudo-MVS, we can provide
#'  the ground truth abundance table, 
#' such that imputation will by done only for pseudo-MVs. This will accelerate 
#' imputation algorithm.
#' @param df Estimate degree of freedom of the IG distribution fitted on 
#' observed variance.
#' @param nu_factor Multiplication factor on degree of freedom. 2 by default.
#' @param max_pg_size Maximum PGs size authorized for imputation. PG size is 
#' plitted if its size is above this threshold.
#' @param min.pg.size2imp Minimum PG size to impute after splitting. PGs for 
#' which size is greater are not imputed. Should be lower than max_pg_size to 
#' have effect. 
#' @param verbose A boolean (FALSE as default) which indicates whether to 
#' display more details ont the process
#' @param ... xxx
#'
#' @return A list containing imputation results for each PG, the execution time,
#'  and adjacency matrix between peptides and PGs corresponding to the 
#'  imputed PGs.
#' @export
#'
#' @examples
#' NULL
#'
impute_block_llk_reset = function(data.pep.rna.crop,
                                  psi,
                                  pep_ab_or = NULL,
                                  df = 1,
                                  nu_factor = 2,
                                  max_pg_size = NULL,
                                  min.pg.size2imp = 1,
                                  verbose = FALSE,
                                  ...) {
  
  adj = data.pep.rna.crop$adj
  if (!is.null(max_pg_size)) {
    adj = split_large_pg(adj, max_pg_size)
  }
  prot.idxs = 1:ncol(adj)
  
  nsamples = nrow(data.pep.rna.crop$peptides_ab)
  logs = list()
  npseudoNA = 0
  begtime = Sys.time()
  if (min.pg.size2imp > 1) { 
    # If need to remove too small PGs, change adj and PG indexes
    pg.idxs = which(colSums(adj) >= min.pg.size2imp)
    adj = adj[, pg.idxs, drop = FALSE]
    prot.idxs = 1:ncol(adj)
  }
  n_params = sum(colSums(adj)^2)
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent 
                                   [Elapsed time: :elapsedfull || Estimated 
                                   time remaining: :eta]",
                         total = n_params,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  for (i in prot.idxs) {
    #message("##### peptide group=", i, "#####")
    if(verbose)
      message("Peptide_group ", i," of ", ncol(adj))
    idx_cur_pep = which(adj[,i] == 1)
    pb$tick(length(idx_cur_pep)^2)
    cur_ab = matrix(data.pep.rna.crop$peptides_ab[,idx_cur_pep], 
                    nrow = nsamples)
    colnames(cur_ab) = colnames(data.pep.rna.crop$peptides_ab)[idx_cur_pep]
    if (is.null(pep_ab_or)) {
      X_gt = NULL
      subpp_ab = as.matrix(cur_ab, nrow = nsamples)
    } else {
      cur_ab_gt = matrix(pep_ab_or[,idx_cur_pep], nrow = nsamples)
      subpp_ab = as.matrix(cur_ab, nrow = nsamples)
      X_gt = as.matrix(cur_ab_gt, nrow = nsamples)
    }
    n_pep_cur = ncol(subpp_ab)
    if (sum(is.na(subpp_ab)) == 0 || # No missing values
        (all(is.na(X_gt) == is.na(subpp_ab)) & !is.null(X_gt))) {# No pseudo-MVs
      logs[[i]] = list()
    } else {
      if (all(is.na(X_gt) == is.na(subpp_ab))) {
        X_gt = NULL
      }
      
      K = (nu_factor*df + n_pep_cur - 1) + n_pep_cur + 1
      psimat = psi*diag(n_pep_cur)
      #res_imp = impfunc(subpp_ab, true_X = X_gt, K = K, psi = psimat, ...) # + max(colSums(is.na(subpp_ab)))
      
      res_imp = py$estimate_params_and_impute(
        subpp_ab, 
        true_X = X_gt, 
        K = K, 
        psi = psimat,
        ...) # + max(colSums(is.na(subpp_ab)))
      
      
      
      ermsg = res_imp$error_msg
      stopifnot(ermsg == "success")
      if ((!is.null(X_gt)) & (ermsg == "success")) {
        npseudos = sum((!is.na(X_gt)) & is.na(subpp_ab))
        npseudoNA = npseudoNA + npseudos
      }
      logs[[i]] = res_imp
    }
  }
  endtime = Sys.time()
  logs[["time"]] = endtime - begtime
  if (!is.null(max_pg_size)) {
    logs[["new_adj"]] = adj
  }
  return(logs)
}



# TODO: Rename this function in code. 
# TODO: Impfunc shoud not be passed as parameter, by we should directly call
# "estimate params & impute.
# TODO: Some parameters are not used anymore in pipeline_ll_imp file, need to 
# remove them
#' @title Impute each PG.
#' @description Imputes each PG separately accounting for transcriptomic 
#' dataset and returns the results for each PG. 
#'
#' @param data.pep.rna.crop A list representing dataset, with mRNA normalized 
#' counts and mRNA/PGs adjacecy table.
# #' @param impfunc Imputation function
#' @param psi Inverse scale parameter for IW prior of peptides abundances
#' @param psi_rna Inverse scale parameter for IW prior of mRNA abundances
#' @param rna.cond.mask Vector of size equal to the number of samples in mRNA 
#' abundance table, 
#' containing indices of conditions of each sample.
#' @param pep.cond.mask Vector of size equal to the number of samples in 
#' peptide abundance table, 
#' containing indices of conditions of each sample.
#' @param pep_ab_or In case we impute a dataset with pseudo-MVS, we can provide 
#' the ground truth abundance table, such that imputation will by done only for 
#' pseudo-MVs. This will accelerate imputation algorithm.
#' @param df Estimate degree of freedom of the IG distribution fitted on 
#' observed variance.
#' @param nu_factor Multiplication factor on degree of freedom. 2 by default.
#' @param max_pg_size Maximum PGs size authorized for imputation. PG size is 
#' plitted if its size is above this threshold.
#' @param max.pg.size2imp Maximum PG size to impute after splitting. PGs for 
#' which size is greater are not imputed. Should be lower than max_pg_size to 
#' have effect. 
#' @param verbose A boolean (FALSE as default) which indicates whether to 
#' display more details ont the process
#' @param ... xxx
#'
#' @return A list containing imputation results for each PG, the execution time,
#'  and adjacency matrix between peptides and PGs corresponding to the 
#'  imputed PGs.
#' @export
#' @import reticulate
#' 
#' @examples
#' NULL
#' 
#'
impute_block_llk_reset_PG = function(data.pep.rna.crop,
                                     psi,
                                     psi_rna,
                                     rna.cond.mask,
                                     pep.cond.mask,
                                     pep_ab_or = NULL,
                                     df = 2,
                                     nu_factor = 1,
                                     max_pg_size = NULL,
                                     max.pg.size2imp = 1,
                                     verbose = FALSE,
                                     ...) {

  if (!is.null(max_pg_size)) {
    adjs = split_large_pg_PG(data.pep.rna.crop$adj, max_pg_size,
                             data.pep.rna.crop$adj_rna_pg)
    adj = adjs$adj
    adj_rna_pg = adjs$adj_rna_pg
  }
  prot.idxs = 1:ncol(adj)
  niter = length(prot.idxs)
  nsamples = nrow(data.pep.rna.crop$peptides_ab)
  logs = list()
  begtime = Sys.time()
  n_pg = length(prot.idxs)
  n_cond = length(unique(rna.cond.mask))
  rnas_ab = matrix(NA, length(pep.cond.mask), nrow(adj_rna_pg))
  colnames(rnas_ab) = colnames(data.pep.rna.crop$rnas_ab)
  
  for (i in unique(rna.cond.mask)) {
    nrep_rna = sum(rna.cond.mask == i)
    nrep_pep = sum(pep.cond.mask == i)
    rnas_means = colMeans(matrix(
      data.pep.rna.crop$rnas_ab[rna.cond.mask == i, ,drop = FALSE], nrep_rna))
    # rnas_sds = apply(matrix(
    #   data.pep.rna.crop$rnas_ab[rna.cond.mask == i, ], nrep_rna), 2, sd, na.rm = TRUE)
    # rnas_sds[is.na(rnas_sds)] = 0
    # rnas_ab[pep.cond.mask == i, ] = matrix(rnorm(
    #   nrow(adj_rna_pg) * nrep_pep, rnas_means, rnas_sds), nrep_pep, byrow = TRUE) # This line enables to sample from statistics of each condition instead of setting the mean.
    rnas_ab[pep.cond.mask == i, ] = matrix(rep(rnas_means, nrep_pep), 
                                           nrep_pep, 
                                           byrow =TRUE)
    
  }
  if (!is.null(max.pg.size2imp)) {
    pg.idxs = which(colSums(adj) <= max.pg.size2imp)
    adj = adj[, pg.idxs, drop = FALSE]
    adj_rna_pg = adj_rna_pg[, pg.idxs, drop = FALSE]
    prot.idxs = 1:ncol(adj)
  }
  n_params = sum((colSums(adj[, prot.idxs, drop = FALSE]) + 1)^2)
  pb <- progress::progress_bar$new(
    format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || 
    Estimated time remaining: :eta]",
                         total = n_params,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  for (i in prot.idxs) {
    if(verbose)
      message("Peptide_group ", i," of ", ncol(adj))
    
    # cat("\n##### PROT ", i, "/ ", n_pg, "#####\n")
    idx_cur_pep = which(adj[,i] == 1)
    pb$tick((length(idx_cur_pep) + 1)^2)
    idx_cur_rna = which(adj_rna_pg[,i] == 1)
    cur_ab = matrix(data.pep.rna.crop$peptides_ab[,idx_cur_pep], 
                    nrow = nsamples)
    colnames(cur_ab) = colnames(data.pep.rna.crop$peptides_ab)[idx_cur_pep]
    cur_ab_rna = matrix(rnas_ab[,idx_cur_rna], nrow = nsamples)
    colnames(cur_ab_rna) = colnames(rnas_ab)[idx_cur_rna]
    subpp_ab = cbind(cur_ab, cur_ab_rna)
    if (is.null(pep_ab_or)) {
      X_gt = NULL
    } else {
      cur_ab_gt = matrix(pep_ab_or[,idx_cur_pep], nrow = nsamples)
      X_gt = cbind(cur_ab_gt, cur_ab_rna)
    }
    if (sum(is.na(subpp_ab)) == 0 |
        (all(is.na(X_gt) == is.na(subpp_ab)) & !is.null(X_gt))) {
      logs[[i]] = list()
    } else {
      if (all(is.na(X_gt) == is.na(subpp_ab))) {
        X_gt = NULL
      }
      n_pep_cur = ncol(subpp_ab)
      K = (nu_factor*df + n_pep_cur - 1) + n_pep_cur + 1
      psimat = c(rep(psi, ncol(cur_ab)), rep(psi_rna, ncol(cur_ab_rna))) * 
        diag(n_pep_cur)
      #res_imp = impfunc(subpp_ab, true_X = NULL, K = K, psi = psimat, ...) # + max(colSums(is.na(subpp_ab)))
      
      res_imp = py$estimate_params_and_impute(subpp_ab, 
                                              true_X = NULL, 
                                              K = K, 
                                              psi = psimat, 
                                              ...) # + max(colSums(is.na(subpp_ab)))
      
      
      
      
      res_imp$Xhat = res_imp$Xhat[, 1:ncol(cur_ab)]
      # res_imp = list(Xhat=matrix(10, nsamples, ncol(cur_ab)), error_msg="success")
      ermsg = res_imp$error_msg
      #print(ermsg)
      stopifnot(ermsg == "success")
      logs[[i]] = res_imp
    }
  }
  endtime = Sys.time()
  logs[["time"]] = endtime - begtime
  if (!is.null(max_pg_size)) {
    logs[["new_adj"]] = adj
  }
  return(logs)
}

# TODO: Change function name
#' @title Impute abundance table from PGs results
#' @description From imputation results in each PG and the associate adjacency 
#' peptide/PG matrix,imputes the original abundance table.  .
#'
#' @param logs.blocks List of PGs imputation results, that also contains 
#' related peptide/PGs adjacency matrix.
#' @param data.pep.rna List representing the dataset not yet imputed
#' @param idx_blocks Indices of PGs for which imputation results should be 
#' integrated
#'
#' @return The original peptide abundance table with imputed values.
#' @export
#'
#' @examples
#' NULL
#' 
#'
impute_from_blocks = function(logs.blocks,
                              data.pep.rna,
                              idx_blocks = NULL) {
  if (!is.null(logs.blocks$new_adj)) {
    adj = logs.blocks$new_adj
  } else {
    adj = data.pep.rna$adj
  }
  if (is.null(idx_blocks)) {
    idx_blocks = 1:ncol(adj)
  }
  npeps = ncol(data.pep.rna$peptides_ab)
  nsamples = nrow(data.pep.rna$peptides_ab)
  pep.imputed = matrix(0, nsamples, npeps)
  colnames(pep.imputed) = colnames(data.pep.rna$peptides_ab)
  rownames(pep.imputed) = rownames(data.pep.rna$peptides_ab)
  n_imputations = matrix(0, nsamples, npeps)
  for (iblock in idx_blocks) {
    cur.block = logs.blocks[[iblock]]
    if (is.list(cur.block)) {
      idxpeps = which(adj[,iblock] == 1)
      if (!is.null(cur.block$Xhat)) {
        pep.imputed[, idxpeps] = pep.imputed[, idxpeps] + cur.block$Xhat
        n_imputations[, idxpeps] = n_imputations[, idxpeps] + 1
      }
    }
  }
  pep.imputed = pep.imputed / pmax(n_imputations, 1)
  pep.imputed[pep.imputed == 0] = data.pep.rna$peptides_ab[pep.imputed == 0]
  return(pep.imputed)
}

