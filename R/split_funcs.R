

#' @title Split too large PGs
#' 
#' @description Split PGs with too many peptides/precursors, while keeping other PG untouched. It creates 
#' new PGs with size equal to size max. Hence, some peptides can be
#'  duplicated in the new PGscreated. The seed is fixed such that the split generated remains constant for a given
#'  adjacency matrix. 
#'  
#'
#' @param adj Adjacency matrix between peptides and PGs.
#' @param size_max Maximum PG size desired.
#'
#' @return New adjacency matrix between peptides and PGs.
#' @export
#'
#' @examples
#' data(subbouyssie)
#' split_large_pg(subbouyssie$adj, 3)
#'
split_large_pg <- function(adj, 
  size_max) {
  
  set.seed(1234)
  
  idx_pg_too_large = which(colSums(adj) >= size_max)
  if (length(idx_pg_too_large) > 0) {
    new_adj = adj
    npeps = nrow(adj)
    l_adjs2bind = list()
    idx_adj2bind = 1
    for (i in 1:length(idx_pg_too_large)) {
      cur_pg_idx = idx_pg_too_large[i]
      cur_pg = adj[, cur_pg_idx, drop = FALSE]
      idx_pg_pep = which(cur_pg == 1)
      n_groups = ceiling(sum(cur_pg) / size_max)
      groups_idx_pg_pep = suppressWarnings(
        split(sample(idx_pg_pep, replace = FALSE), 1:n_groups)
      )
      for (j in 1:length(groups_idx_pg_pep)) {
        new_pg = rep(FALSE, npeps)
        idx2add = setdiff(idx_pg_pep, groups_idx_pg_pep[[j]])
        n2samples = size_max - length(groups_idx_pg_pep[[j]])
        idx2add = sample(idx2add, n2samples)
        new_pg[c(groups_idx_pg_pep[[j]], idx2add)] = TRUE
        l_adjs2bind[[idx_adj2bind]] = new_pg
        idx_adj2bind = idx_adj2bind + 1
      }
    }
    return(cbind(new_adj, do.call(cbind, l_adjs2bind))[, -idx_pg_too_large])
  } else {
    return(adj)
  }
}





# TODO: Rename the function
#' @title Splits too large PGs in proteogenomics context
#' 
#' @description Split PGs with too many peptides/precursors while keeping other PGs untouched,
#' and adapts adjacency matrix between mRNA and PGs accordingly. It creates new PGs 
#' with size equal to size max (including peptides and mRNAs). Hence, some 
#' peptides and mRNA can be duplicated in the new PGs.
#' The seed is fixed such that the split generated remains constant for a given
#' adjacency matrix. 
#'
#' @param adj Adjacency matrix between peptides and PGs.
#' @param size_max Maximum PG size desired.
#' @param adj_rna_pg Adjacency matrix between mRNA and PGs.
#'
#' @return List containing new adjacency matrix between peptides and PGs, 
#' and new adjacency matrix between mRNA and PGs. 
#' @export
#'
#' @examples
#' NULL
#'
split_large_pg_PG <- function(adj, 
  size_max, 
  adj_rna_pg) {
  
  set.seed(1234)
  
  idx_pg_too_large = which((colSums(adj) + colSums(adj_rna_pg)) >=
      size_max)
  if (length(idx_pg_too_large) > 0) {
    new_adj = adj
    new_adj_rna_pg = adj_rna_pg
    nrnas = nrow(adj_rna_pg)
    npeps = nrow(adj)
    adj2bind = matrix(NA, nrow = npeps, ncol = 0)
    adj_rna_pg_2bind = matrix(NA, nrow = nrnas, ncol = 0)
    for (i in 1:length(idx_pg_too_large)) {
      cur_pg_idx = idx_pg_too_large[i]
      cur_pg = adj[, cur_pg_idx, drop = FALSE]
      idx_pg_pep = which(cur_pg == 1)
      cur_rna_pg = matrix(adj_rna_pg[, cur_pg_idx, drop = FALSE], ncol = 1)
      colnames(cur_rna_pg) = colnames(adj_rna_pg)[cur_pg_idx]
      n_rna_cur_pg = sum(cur_rna_pg)
      stopifnot("Size max of pg too low when including RNA expression" =
          n_rna_cur_pg < size_max)
      n_groups = ceiling(sum(cur_pg) / (size_max - n_rna_cur_pg) )
      groups_idx_pg_pep = suppressWarnings(
        split(sample(idx_pg_pep, replace = FALSE), 1:n_groups)
      )
      for (j in 1:length(groups_idx_pg_pep)) {
        new_pg = rep(0, npeps)
        idx2add = setdiff(idx_pg_pep, groups_idx_pg_pep[[j]])
        n2samples = size_max - n_rna_cur_pg - length(groups_idx_pg_pep[[j]])
        idx2add = sample(idx2add, n2samples)
        new_pg[c(groups_idx_pg_pep[[j]], idx2add)] = 1
        adj2bind = cbind(adj2bind, new_pg)
        adj_rna_pg_2bind = cbind(adj_rna_pg_2bind, cur_rna_pg)
      }
    }
    return(list(adj = cbind(new_adj, adj2bind)[, -idx_pg_too_large],
      adj_rna_pg = cbind(new_adj_rna_pg,
        adj_rna_pg_2bind)[, -idx_pg_too_large]))
  } else {
    return(list(adj = adj, adj_rna_pg = adj_rna_pg))
  }
}