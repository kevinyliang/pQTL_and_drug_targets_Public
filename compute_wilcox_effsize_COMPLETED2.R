#' for computing effsize suitable for wilcox rank sum test
#' 
#' Notes:
#'  the rank biserial correlation computed here based on common language effect size (f-(1-f)) is the same as the glass rank biserial correlation computed in rcompanion::wilcoxonRG()

.compare_x_geq_y_with_ties <- function(x,y){
  greater <- x>y
  equals <- x==y
  return(greater + (equals/2))
}

#' compute_ref_geq_effsizes
#' compute the common language effect size as described in 
#'  Statistical analysis of ordinal categorical status after therapies
#'  https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Effect_sizes
#'
#' Notes: assumes only 2 groups.
#' 
#' @params metrics: vector metric of all data
#' @params groups: vector of classification of all data
#'  
compute_ref_geq_effsizes <- function(metrics,groups,refgroup){
  assertthat::assert_that(
    length(unique(groups))==2 &
    length(metrics) == length(groups) & 
    sum(is.na(metrics)) == 0 &
    sum(is.na(groups)) == 0 &
    refgroup %in% groups
  )
  
  g2 <- setdiff(unique(groups),refgroup)
  g1_data <- metrics[groups==refgroup]
  g2_data <- metrics[groups==g2]
  assertthat::assert_that(
    length(g1_data) + length(g2_data) == length(metrics)
  )
  #' compute CLES of Pr(G1 > G2)
  #' # This would be equivalent to AUROC
  g1_cles <- sum(outer(g1_data,g2_data,.compare_x_geq_y_with_ties))/(length(g1_data) * length(g2_data))
  g2_cles <- sum(outer(g2_data,g1_data,.compare_x_geq_y_with_ties))/(length(g1_data) * length(g2_data))
  assertthat::assert_that(g1_cles + g2_cles == 1)

  #' compute rbs of Pr(G1 > G2)
  #' this simplified to 2(AUROC) - 1
  g1_rbs <- g1_cles - (1-g1_cles)

  return(data.frame(
    refgroup_cles_geq=g1_cles,
    refgroup_rbs_geq=g1_rbs
  ))
}


#' compute_ustat
#' computes the rank sum as
#' 
#' U1 = n1n2 + ((n1 * (n1+1))/2) - R1
#' 
#' Interpretation as proportion of G1 < G2
.compute_ustat <- function(n1,n2,r1){
  a <-  n1*n2
  b <-  (n1 * (n1+1))/2
  u <- a + b - r1
  return(u)
}
#' compute U stat and ranksum
#' 
#' obtain the U statistic that is reported by wilcox.test and the rank sum
#' 
compute_U_stat_ranksum_ref_smaller <- function(metrics,groups,refgroup){
  assertthat::assert_that(
    length(unique(groups))==2 &
    length(metrics) == length(groups) & 
    sum(is.na(metrics)) == 0 &
    sum(is.na(groups)) == 0 &
    refgroup %in% groups
  )
  g2 <- setdiff(unique(groups),refgroup)
  
  data_ranks <- rank(metrics)
  g1_ranksum <- sum(data_ranks[groups==refgroup])
  g2_ranksum <- sum(data_ranks[groups==g2])

  n1 <- sum(groups == refgroup)
  n2 <- sum(groups == g2)
  g1_u <- .compute_ustat(n1,n2,g1_ranksum)
  g2_u <- .compute_ustat(n2,n1,g2_ranksum)
  return(
    data.frame(
      refgroup_ranksum = g1_ranksum,
      refgroup_u = g1_u,
      refgroup_n = n1,
      g2_ranksum = g2_ranksum,
      g2_u = g2_u,
      g2_n = n2,
      min_u = min(g1_u,g2_u)
    )
  )
}