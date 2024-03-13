#' @title COnvert Pirat dataset to SummarizedExperiment
#' @description This function converts the original dataset structure into 
#' a SummarizedExperiment .
#' 
#' @param peptides_ab xxx
#' @param adj xxx
#' @param mask_prot_diff xxx
#' @param mask_pep_diff xxx
#' 
#' @return An instance of the class `SummarizedExperiment`
#' 
#' @import SummarizedExperiment
#' 
#' @export
#' 
#' @examples
#' data(subbouyssie)
#' peptides_ab <- subbouyssie$peptides_ab
#' adj <- subbouyssie$adj
#' mask_prot_diff <- subbouyssie$mask_prot_diff
#' mask_pep_diff <- subbouyssie$mask_pep_diff
#' obj <- pirat2SE(peptides_ab, adj, mask_prot_diff, mask_pep_diff )
#' obj
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment
#' 
pirat2SE <- function(peptides_ab, 
                     adj, 
                     mask_prot_diff, 
                     mask_pep_diff){

 obj <- SummarizedExperiment::SummarizedExperiment(
   assays = as.matrix(t(peptides_ab), row.names = colnames(t(peptides_ab))), 
   colData = data.frame(Condition = colnames(t(peptides_ab)),
                        Sample = colnames(t(peptides_ab)),
                        row.names = colnames(t(peptides_ab)))
   )
 
 metadata(obj)$adj <- adj
 metadata(obj)$mask_prot_diff <- mask_prot_diff
 metadata(obj)$mask_pep_diff <- mask_pep_diff
  
  obj
}




#' @title xxx
#' @description xxx
#' 
#' @param se An instance of the class SummarizedExperiment
#' @param ... Additional arguments to pass to `pipeline_llkimpute()`
#' 
#' 
#' @export
#' 
#' @examples
#' data(subbouyssie)
#' obj <- pirat2SE(subbouyssie$peptides_ab, subbouyssie$adj, subbouyssie$mask_prot_diff, 
#' subbouyssie$mask_pep_diff )
#' res <- wrapper_pipeline_llkimpute(obj)
#' 
#' @return See pipeline_llkimpute() function
#' 
#' 
wrapper_pipeline_llkimpute <- function(se, ...){
  
  stopifnot(inherits(obj, 'SummarizedExperiment'))
  
  obj <- list(
    peptides_ab = t(assay(se)),
    adj = metadata(se)$adj,
    mask_prot_diff = metadata(se)$mask_prot_diff,
    mask_pep_diff = metadata(se)$mask_pep_diff
  )

  my_pipeline_llkimpute(obj, ...)
}


