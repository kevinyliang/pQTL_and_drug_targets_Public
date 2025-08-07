gather_data <- function(metrics,groups,refgroup,method,metric_name,target_src){
  
  assertthat::assert_that(
    length(unique(groups))==2 &
    length(metrics) == length(groups) & 
    sum(is.na(metrics)) == 0 &
    sum(is.na(groups)) == 0 &
    refgroup %in% groups
  )
  g2 <- setdiff(unique(groups),refgroup)
  #' wilcox test
  wilcox_res <- wilcox.test(
    metrics[groups == refgroup],
    metrics[groups == g2],
    conf.int=T
  )
  #' CLES = AUROC and own Rbs
  wilcox_effsize <- compute_ref_geq_effsizes(
    metrics,groups,refgroup
  )
  #' AUROC
  auroc <- roc.curve(
    scores.class0 = metrics,
    weights.class0 = as.numeric(groups==refgroup),
    curve=T
  )
  assertthat::assert_that(
    all.equal(wilcox_effsize$refgroup_cles_geq,auroc$auc)
  )
  wilcox_effsize_boot <- boot(
    data = data.frame(metrics=metrics,groups=groups,refgroup = refgroup),
    statistic = function(df,indicies){
      res <- compute_ref_geq_effsizes(
        metrics = df$metrics[indicies],
        groups = df$groups[indicies],
        refgroup = unique(df$refgroup)
      )
      return(c("auroc" = res$refgroup_cles_geq,'rbs' = res$refgroup_rbs_geq))
    },
    parallel = 'multicore',
    ncpus=10,
    R=n_bootstrap
  )
  assertthat::assert_that(
    wilcox_effsize$refgroup_cles_geq == wilcox_effsize_boot$t0['auroc'] &
    wilcox_effsize$refgroup_rbs_geq == wilcox_effsize_boot$t0['rbs']
  )
  auroc_ci <- quantile(wilcox_effsize_boot$t[,1],c(0.025,0.975))
  rbs_ci <- quantile(wilcox_effsize_boot$t[,2],c(0.025,0.975))
  #' own wilcox stat
  wilcox_stat <- compute_U_stat_ranksum_ref_smaller(
    metrics,groups,refgroup
  )
  #' Glass Rbs
  glass_rbs <- rcompanion::wilcoxonRG(
    metrics,
    groups != refgroup
  )

  AUPRC <- pr.curve(
    scores.class0 = metrics,
    weights.class0 = as.numeric(groups==refgroup),
    curve=T
  )
  auprc_boot <- boot(
    data = data.frame(metrics = metrics,groups=as.numeric(groups==refgroup)),
    statistic = function(df,idx){
      res <- pr.curve(
        scores.class0 = df$metrics[idx],
        weights.class0 = df$groups[idx],curve=T
      )
      return(res$auc.integral)
    },
    parallel='multicore',
    R=n_bootstrap
  )
  assertthat::assert_that(
    auprc_boot$t0 == AUPRC$auc.integral
  )
  auprc_ci <- quantile(auprc_boot$t,c(0.025,0.975))
  df <- data.frame(
      #' sample size
      causal_gene_counts = wilcox_stat$refgroup_n,
      non_causal_gene_counts = wilcox_stat$g2_n,
      #' wilcox stat
      wilcox_stat = wilcox_res$statistic,
      own_stat_refgroup = wilcox_stat$refgroup_u,
      own_stat_altgroup = wilcox_stat$g2_u,
      wilcox_pval = wilcox_res$p.value,
      #' CLES
      AUROC = wilcox_effsize$refgroup_cles_geq,
      auroc_lower_ci =auroc_ci[1],
      auroc_upper_ci =auroc_ci[2],
      #' Rbs
      Rbs = wilcox_effsize$refgroup_rbs_geq,
      rbs_lower_ci = rbs_ci[1],
      rbs_upper_ci = rbs_ci[2],
      glass_rbg =  glass_rbs,
      #' AUPRC
      AUPRC = AUPRC$auc.integral,
      auprc_lower_ci = auprc_ci[1],
      auprc_upper_ci = auprc_ci[2],
      auprc_null = sum(as.numeric(groups==refgroup))/length(groups)
    ) %>% 
    dplyr::mutate(
      method = method,
      metric_name = metric_name,
      target_src = target_src
    ) %>% 
    dplyr::select(
      method,metric_name,target_src,
      causal_gene_counts,
      non_causal_gene_counts,
      wilcox_pval,
      AUROC,auroc_lower_ci,auroc_upper_ci,
      AUPRC,auprc_lower_ci,auprc_upper_ci,
      Rbs,rbs_lower_ci,rbs_upper_ci,
      wilcox_stat,own_stat_refgroup,own_stat_altgroup,
      glass_rbg,auprc_null
    )
  return(
    df 
  )
}