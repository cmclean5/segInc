#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# RStudio Workbench is strictly for use by Public Health Scotland staff and     
# authorised users only, and is governed by an <Acceptable Usage Policy>.
#
# This is a shared resource and is hosted on a pay-as-you-go cloud computing
# platform.  Your usage will incur direct financial cost to Public Health
# Scotland.  As such, please ensure
#
#   1. that this session is appropriately sized with the minimum number of CPUs
#      and memory required for the size and scale of your analysis;
#   2. the code you write in this script is optimal and only writes out the
#      data required, nothing more.
#   3. you close this session when not in use; idle sessions still cost PHS
#      money!
#
# For further guidance, please see <insert link>.
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Set up some formatting for flextable to be applied to most tables
my_ft_format <- function(ft) {
  ft %>%
    bold(part = "header") %>%
    bg(bg = "#43358B", part = "header") %>%
    color(color = "white", part = "header") %>%
    align(align = "left", part = "header") %>%
    valign(valign = "center", part = "header") %>%
    valign(valign = "top", part = "body") %>%
    colformat_num(big.mark = "") %>%
    fontsize(size = 12, part = "all") %>%
    border(border = fp_border_default(color = "#000000", width = 0.5),
           part = "all")
}


get_segi_doll_pop_model <- function(){
  setNames( c(
    12000, 10000, 9000, 9000,
    8000, 8000, 6000, 6000,
    6000, 6000, 5000, 4000,
    4000, 3000, 2000, 1000,
    500, 500) , c(
      "0-4", "5-9", "10-14", "15-19",
      "20-24", "25-29", "30-34", "35-39",
      "40-44", "45-49", "50-54", "55-59",
      "60-64", "65-69", "70-74", "75-79",
      "80-84", "85+"
    ))
}

get_emp2013_pop_model <- function(){
  setNames( c(
    5000, 5500, 5500, 5500,
    6000, 6000, 6500, 7000,
    7000, 7000, 7000, 6500,
    6000, 5500, 5000, 4000,
    2500, 2500), c(
      "0-4", "5-9", "10-14", "15-19",
      "20-24", "25-29", "30-34", "35-39",
      "40-44", "45-49", "50-54", "55-59",
      "60-64", "65-69", "70-74", "75-79",
      "80-84", "85+"
    ))
}


set_pop_model <- function(model="segi"){
  switch(model,
         "segi" = get_segi_doll_pop_model(),
         "epm13"= get_emp2013_pop_model())
}


compute_age_group_asr <- function(x,  
                                  early_age_bands     = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39"),
                                  early_sub_age_bands = c("40-44", "45-49"),
                                  late_age_bands      = c("50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+"),
                                  pop_model = "segi", 
                                  conf_level = 0.95){
  
  z = qnorm(1-(1-conf_level)/2)
  
  ## Choose age standardisation model: c("segi", "epm13")
  pop_model = set_pop_model(model=pop_model)
  
  x %<>%
    mutate( pop_weight      = recode(age_group, !!!pop_model)) %>%
    mutate( crude_rate      = (num/dem)*1e5,
            weighted_rate   = crude_rate * pop_weight) %>%
    group_by(year, group) %>%
    summarise(
      sum_weighted_rate   = sum(weighted_rate),
      total_cases         = sum(num),
      total_pop           = sum(dem),
      ## Calculate variance components for CI
      variance_sum        = sum( (pop_weight^2 * num) / (dem^2) ), .groups="drop") %>%
    mutate( early_tot     = sum(pop_model[names(pop_model) %in% early_age_bands]),
            sub_early_tot = sum(pop_model[names(pop_model) %in% early_sub_age_bands]),
            late_tot      = sum(pop_model[names(pop_model) %in% late_age_bands])) %>%
    ## Total standard population for this age group 
    mutate( total_std_pop  = case_when( group=="early"     ~ early_tot, 
                                        group=="sub_early" ~ sub_early_tot, 
                                        group=="late"      ~ late_tot)
    ) %>% 
    ## Calculate ASR and it's 95% CI
    mutate(
      
      ## Calculate ASR
      asr = sum_weighted_rate / total_std_pop,
      
      ## Dobson method for calculating 95% CI
      ## Calculate variance of ASR
      var_asr  = (variance_sum / total_std_pop^2) * 1e5^2, 
      
      ## Calculate standard error
      se_asr   = sqrt(var_asr), 
      
      ## observed lower confidence limit
      obs_lw = total_cases * ( 1 - 1/(9*total_cases) - z/(3*sqrt(total_cases)))^3,
      
      ## observed upper confidence limit
      obs_up = (total_cases + 1) * ( 1 - 1/(9*(total_cases+1)) + z/(3*sqrt((total_cases+1))))^3,
      
      ## asr_segi lower
      asr_lw = asr + sqrt(var_asr/total_cases) * ( obs_lw - total_cases),
      
      ## asr_segi upper
      asr_up = asr + sqrt(var_asr/total_cases) * ( obs_up - total_cases),
      
      ## Calculate 95% confidence interval
      asr_ci_lw_dobson = asr - z * se_asr,
      asr_ci_up_dobson = asr + z * se_asr,
      
      ## Gamma-based confidence-intervals - should be more accurate for small case counts?
      asr_ci_lw_gamma = if_else(total_cases > 0, 
                                asr * qgamma(0.025, shape=total_cases,   rate=1)/total_cases,
                                0),
      
      asr_ci_up_gamma = if_else(total_cases > 0, 
                                asr * qgamma(0.975, shape=total_cases+1, rate=1)/total_cases,
                                asr * qgamma(0.975, shape=1, rate=1)),
      
     crude_rate = (total_cases / total_pop ) * 1e5) %>%
     ## For plotting
     mutate( diag_year = year ) %>%
     dplyr::select(-c(early_tot, sub_early_tot, late_tot)) %>%
     mutate(group = fct_relevel(group, "late", "sub_early", "early"))
  
   return(x)
  
}

incidence_rates_ci5 <- function(x, cancer_id= NULL, country_id="82602099", model_type="segi",
                                early_age_bands     = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39"),
                                early_sub_age_bands = c("40-44", "45-49"),
                                late_age_bands      = c("50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+"),
                                conf_level = 0.95){
  
  Ncodes = length(country_id)

  if( !is.null(cancer_id) ){
    x %<>% filter(cancer_code %in% cancer_id)
  }
    
  ## We should only include years that include data entries for each nations (sub)-registry
  valid_years <- x %>%
    filter( id_code %in% country_id ) %>%
    dplyr::select(id_code, year) %>%
    distinct() %>%
    group_by(year) %>%
    count() %>%
    mutate( use_year = ifelse( n==Ncodes, TRUE, FALSE))
  
  
  x %>%
    filter( id_code %in% country_id ) %>%
    filter( year %in% valid_years$year[valid_years$use_year]) %>%
    ## sum over cancer groups
    group_by(age, year, sex) %>% 
    summarise(num = sum(cases), dem=sum(py), .groups="drop") %>%
    mutate( age_group = recode(as.character(age), !!!age_band_map)) %>%
    dplyr::select(-age) %>%
    mutate( group = case_when( age_group %in% early_age_bands     ~ "early",
                               age_group %in% early_sub_age_bands ~ "sub_early",
                               age_group %in% late_age_bands      ~ "late")) %>%
    filter(!is.na(group)) -> ci5_crude_rate_by_age_band

  
 ## Male
 ci5_rate_male <- ci5_crude_rate_by_age_band %>%
                    filter(sex==1) %>%
                    compute_age_group_asr(x=.,
                                          pop_model           = model_type,
                                          early_age_bands     = early_age_bands,
                                          early_sub_age_bands = early_sub_age_bands,
                                          late_age_bands      = late_age_bands,
                                          conf_level          = conf_level)
 
 
 ## Female
 ci5_rate_female <- ci5_crude_rate_by_age_band %>%
   filter(sex==2) %>%
   compute_age_group_asr(x=.,
                         pop_model           = model_type,
                         early_age_bands     = early_age_bands,
                         early_sub_age_bands = early_sub_age_bands,
                         late_age_bands      = late_age_bands,
                         conf_level          = conf_level)
 
   
  
  ## Both Sexes
  ci5_final_rate_est <- compute_age_group_asr(x=ci5_crude_rate_by_age_band,
                                              pop_model           = model_type,
                                              early_age_bands     = early_age_bands,
                                              early_sub_age_bands = early_sub_age_bands,
                                              late_age_bands      = late_age_bands,
                                              conf_level          = conf_level)

  ci5_final_rate_est %<>%
    mutate( diag_year      = year)
  
 
  return(list(age_bands    = ci5_crude_rate_by_age_band, 
              rates        = ci5_final_rate_est,
              rates_male   = ci5_rate_male,
              rates_female = ci5_rate_female))
  
}


ci5_data_fit <- function(x, est_var="asr"){
  
  
  ## linear regression models for each flag
  x %>%
    filter(!!sym(est_var) > 0) %>%
    mutate( est    = !!sym(est_var)) %>%
    mutate(log_est = log(est)) %>% 
    mutate( asr_ci_lw = asr_ci_lw_gamma,
            asr_ci_up = asr_ci_up_gamma) %>%
    dplyr::select(diag_year, group, log_est, total_cases, total_pop, asr_ci_lw, asr_ci_up) -> df
  
  dat = list()
  
  jp_early_models     = NULL
  jp_sub_early_models = NULL
  jp_late_models      = NULL
  
  k=1
  
  ## parameter set
  param_set = c("diag_year", "group", "log_est", "total_cases", "total_pop")
  
  if( "early" %in% levels(x$group) ){
  
  ## early onset
  dat[[k]] = df %>% filter(group=="early") %>% mutate(group ="early") %>% dplyr::select(param_set) # diag_year, group, log_est)
  
  jp_early_models = joint_point_model(x=dat[[k]])
  
  k=k+1
  
  } 
  
  if( "sub_early" %in% levels(x$group) ){
  
  ## sub_early onset
  dat[[k]] = df %>% filter(group=="sub_early") %>% mutate(group ="sub_early") %>% dplyr::select(param_set)
    #dplyr::select(diag_year, group, log_est)
  
  jp_sub_early_models = joint_point_model(x=dat[[k]])
  
  k=k+1
  
  }
  
  if( "late" %in% levels(x$group) ){
  
  ## late onset
  dat[[k]] = df %>% filter(group=="late") %>%  mutate(group ="late") %>% dplyr::select(param_set) #dplyr::select(diag_year, group, log_est)
  
  jp_late_models = joint_point_model(x=dat[[k]])
  
  }
  
  return(list(df        = df, 
              early     = jp_early_models,
              sub_early = jp_sub_early_models,
              late      = jp_late_models))

}


ci5_data_fit_plot <- function(x,
                              model_fit = list("late"=list("fit"="jp_none", "labs"="late, > 50"),
                                            "sub_early"=list("fit"="jp_none", "labs"="early, 40-49"),
                                            "early"=list("fit"="jp_none", "labs"="early, 20-39")),
                              y_axis_breaks   = c(0.1,1,2,3,4,5,7,10,20,30,50,70,100,200),
                              show_error_band = FALSE){
  
  gp = NULL
  
  Nm = length(names(model_fit))
  
  if(!any(names(model_fit) %in% "sub_early")){
  
    model_ci_indx = rep(NA,Nm)
    
    for( i in 1:Nm ){
      model_ci_indx[i] = switch(model_fit[[i]]$fit,
                                 "jp_none" = 1,
                                 "jp_one"  = 2,
                                 "jp_two"  = 3)
    }
    
    ## Plot
    df_fit = x$df %>% left_join(dplyr::bind_rows(x$early$dat[[1]],
                                                 x$late$dat[[1]]), 
                                by=c("diag_year", "group", "log_est")) %>%
      mutate(group = fct_relevel(group, "late", "early")) 
    
    df_fit %>%
      rename( est     = log_est) %>%
      mutate( est     = exp(est),
              jp_none = exp(jp_none),
              jp_one  = exp(jp_one),
              jp_two  = exp(jp_two)) %>% 
      mutate( fit = case_when( group == "early"     ~ !!sym(model_fit$early$fit),    
                               group == "late"      ~ !!sym(model_fit$late$fit))) %>% 
      ggplot(., aes(x=diag_year, y=est, color=group))+
      geom_point()+ 
      geom_line(aes(y=fit))+ 
      {if(show_error_band) geom_ribbon(aes(ymin=asr_ci_lw, ymax=asr_ci_up), alpha=0.2, color=NA)}+
        scale_color_manual(values=c("early"     = "#FF0000", 
                                    "late"      = "blue"),
                           labels=c("early"     = model_fit$early$lab, 
                                    "late"      = model_fit$late$lab))+
      scale_y_log10(
        limits = range(y_axis_breaks), 
        breaks = y_axis_breaks,
        labels = label_number(accuracy = 1)
        #minor_breaks = log10(1:10 %o% 10^(0:3))
      )+
      labs(
        title= "", 
        x="Diagnosis Year",
        y="Age-Standardised rate per 100,000"
      )+
      theme_minimal()+
      #comp_axis_theme+
      summary_axis_theme2+
      theme(
        panel.grid.major.y = element_line(color="#EBEBEB"),
        panel.grid.minor.y = element_blank() ) -> gp
    
    
  } else {
    
    model_ci_indx = rep(NA,Nm)
    
    for( i in 1:Nm ){
      model_ci_indx[i] = switch( model_fit[[i]]$fit,
                                 "jp_none" = 1,
                                 "jp_one"  = 2,
                                 "jp_two"  = 3)
    }
    
    ## Plot
    df_fit = x$df %>% left_join(dplyr::bind_rows(x$early$dat[[1]],
                                                 x$sub_early$dat[[1]],
                                                 x$late$dat[[1]]), 
                                by=c("diag_year", "group", "log_est")) %>%
      mutate(group = fct_relevel(group, "late", "sub_early", "early")) 
    
    df_fit %>%
      rename( est     = log_est) %>%
      mutate( est     = exp(est),
              jp_none = exp(jp_none),
              jp_one  = exp(jp_one),
              jp_two  = exp(jp_two)) %>% 
      mutate( fit = case_when( group == "early"     ~ !!sym(model_fit$early$fit), 
                               group == "sub_early" ~ !!sym(model_fit$sub_early$fit), 
                               group == "late"      ~ !!sym(model_fit$late$fit))) %>% 
      ggplot(., aes(x=diag_year, y=est, color=group))+
      geom_point()+ 
      geom_line(aes(y=fit))+ 
      {if(show_error_band) geom_ribbon(aes(ymin=asr_ci_lw, ymax=asr_ci_up), alpha=0.2, color=NA)}+
      scale_color_manual(values=c("early"     = "#FF0000", 
                                  "sub_early" = "#8B0000", 
                                  "late"      = "blue"),
                         labels=c("early"     = model_fit$early$lab, 
                                  "sub_early" = model_fit$sub_early$lab, 
                                  "late"      = model_fit$late$lab))+
      scale_y_log10(
        limits = range(y_axis_breaks), 
        breaks = y_axis_breaks,
        labels = label_number(accuracy = 1)
        #minor_breaks = log10(1:10 %o% 10^(0:3))
      )+
      labs(
        title= "", 
        x="Diagnosis Year",
        y="Age-Standardised rate per 100,000"
      )+
      theme_minimal()+
      #comp_axis_theme+
      summary_axis_theme2+
      theme(
        panel.grid.major.y = element_line(color="#EBEBEB"),
        panel.grid.minor.y = element_blank() ) -> gp
    
    
  }
  
  return(gp)
  
}


compute_aapc_linear <- function(fit){
  beta1    <- coef(fit)[2]
  se_beta1 <- summary(fit)$coefficients[2,"Std. Error"]
  aapc     <- (exp(beta1)-1) * 100
  se_appc  <- 100*exp(beta1) * se_beta1
  ci       <- confint(fit)[2,]
  aapc_l   <- (exp(ci[1])-1) * 100
  aapc_u   <- (exp(ci[2])-1) * 100
  return( c(aapc, se_appc, aapc_l, aapc_u))
}

compute_jps <- function(fit){
  years    = round(fit$psi[,2],1)
  se_years = round(fit$psi[,3],1)
  
  result = c(NA,NA)
  
  if( length(years)==1 ){
    result = c(years, se_years)
  } 
  
  if( length(years)==2 ){
    result[1] = paste0(years,    collapse=" ,")
    result[2] = paste0(se_years, collapse=" ,")
  }
  
  return(as.character(result))
  
}

## Modify R's segmented package aapc function 
compute_aapc_interval <- function(seg_obj, from_year, to_year, conf.level=0.95, compare=FALSE){
  
  compute_aapc <- function(slope){( exp(slope)-1) * 100 }
  
  blockdiag <- function(...) {
    args <- list(...)
    args <- lapply(args, function(x) {
      if( is.vector(x)) x <- matrix(x, nrow=1)
      if( is.null(dim(x))) stop("args must be matrix or convertible to matrix")
      x
    })
    
    nc <- sapply(args, ncol)
    cumnc <- cumsum(nc)
    NC <- sum(nc)
    rowfun <- function(m, zbefore, zafter) {
      cbind(matrix(0, ncol = zbefore, nrow = nrow(m)), 
            m, matrix(0, ncol = zafter, nrow = nrow(m)))
    }
    ret <- rowfun(args[[1]], 0, NC - ncol(args[[1]]))
    for (i in 2:length(args)) {
      ret <- rbind(ret, rowfun(args[[i]], cumnc[i - 1], 
                               NC - cumnc[i]))
    }
    ret
  }
  
  
  COV     <- vcov(seg_obj)
  estcoef <- coef(seg_obj)
  
  nomeZ       <- seg_obj$nameUV$Z
  term        <- nomeZ[1]
  nomi.psi    <- grep(paste("\\.", term, sep = ""), seg_obj$nameUV$V, value = TRUE)
  nomi.dslope <- grep(paste("\\.", term, sep = ""), seg_obj$nameUV$U, value = TRUE)
  
  null.left   <- TRUE
  if (term %in% names(estcoef)) {
    nomi.dslope <- c(term, nomi.dslope)
    null.left <- FALSE
  }
  
  
  ## Slopes and breakpionts 
  est.slopes  = slope(seg_obj)[[1]][,"Est."]
  break_years = seg_obj$psi[nomi.psi,2]
  all_years   = sort(unique(c(from_year, break_years, to_year)))
  seg_start   = head(all_years, -1)
  seg_end     = tail(all_years, -1)
  
  start_year  = min(seg_obj$rangeZ[,term])
  end_year    = max(seg_obj$rangeZ[,term])
  
  ## Identify overlapping segments
  in_window     <- (seg_start < to_year) & (seg_end > from_year)
  overlap_start <- pmax(seg_start[in_window], from_year)
  overlap_end   <- pmin(seg_end[in_window], to_year)
  
  weights       <- (overlap_end - overlap_start) / (to_year - from_year)
  slopes        <- est.slopes[in_window]
  psi           <- sum(weights * slopes)
  mu            <- compute_aapc(sum(weights * slopes))
  
  
  ## Delta methods for SE
  k           <- length(break_years)
  A           <- matrix(0, k + 1, k + 1)
  A[row(A) >= col(A)] <- 1
  B           <- diff(diag(k + 2), diff = 1)/(end_year - start_year)
  xsi         <- c(crossprod(weights, A[in_window,, drop=FALSE]), crossprod(slopes, B[in_window,, drop=FALSE]))
  
  cof         <- estcoef[nomi.dslope]
  v.DeltaPsi  <- COV[c(nomi.dslope, nomi.psi), c(nomi.dslope, nomi.psi)]
  
  rownames(v.DeltaPsi) <- colnames(v.DeltaPsi) <- c(names(cof), nomi.psi)
  
  v.delta <- v.DeltaPsi[names(cof), names(cof)]
  
  VC      <- v.DeltaPsi[nomi.psi, names(cof)]
  v.psi   <- as.matrix(COV[nomi.psi, nomi.psi])
  VV      <- blockdiag(v.delta, diag(1) * 0, v.psi, diag(1) * 0)
  
  id.cov1 <- 1:length(est.slopes)
  id.cov2 <- seq.int((length(est.slopes)+2), length.out=length(break_years))
  
  #VC <- cbind(0,VC)
  #VV <- rbind(0, cbind(0,VV))
  
  VV[id.cov2, id.cov1] <- VC
  VV[id.cov1, id.cov2] <- t(VC)
  
  se.mu  <- sqrt(drop(xsi %*% VV %*% xsi))
  se.psi <- sqrt(drop(xsi %*% VV %*% xsi))
  se.mu <- compute_aapc(se.mu)
  z     <- abs(qnorm((1 - conf.level)/2))
  r     <- c(Est = mu, St.Err = se.mu, mu + c(-z, z) * se.mu)
  cin   <- paste("CI", "(", conf.level * 100, "%", ")", c(".l", ".u"), sep = "")
  names(r) <- c("Est.", "St.Err", cin) 
  
  if(compare){
    r   <- c(Psi = psi, Psi.Err = se.psi, Est = mu, St.Err = se.mu, mu + c(-z, z) * se.mu)
    cin   <- paste("CI", "(", conf.level * 100, "%", ")", c(".l", ".u"), sep = "")
    names(r) <- c("Psi","Psi.Err","Est.", "St.Err", cin) 
  }
  
  return(r)
  
}



get_aapc <- function(fit, cancer="", country="", groups=c("early"="sg1", "sub_early"="sg1", "late"="sg1"),
                     age_bands=c("early"="20-39 years", "early"="40-49 years", "late"="50-85+ years"),
                     aapc_years=10){
  
  col_names = c("country", "cancer", "group", "age_band", "Joint-Point (Year)", "Joint-Point (Std. Err)", "from_year", "to_year",
                "AAPC (Est)", "AAPC (Std. Err)", "CI(95%).l", "CI(95%).u")
  
  df = setNames(data.frame(matrix(ncol=length(col_names), nrow=length(groups))), col_names)
  
  df$country  = country
  df$cancer   = cancer
  df$group    = names(groups)
  df$age_band = age_bands
  
  diag_years = unique(fit$df$diag_year)
  n_years    = length(diag_years)
  to_year    = diag_years[n_years]
  from_year  = to_year-aapc_years #10
  
  df$"from_year" = from_year
  df$"to_year"   = to_year
  
  for( i in 1:length(groups) ){
    
    model = groups[i]
    
    indx = which(names(fit)==names(groups)[i])
    
    jps <- switch(model[[1]],
                  "ll"  = {c("-","-")},
                  "sg1" = {compute_jps(fit[[indx]]$sg1[[1]])},
                  "sg2" = {compute_jps(fit[[indx]]$sg2[[1]])}
    )
    
    result <- switch(model[[1]],
                     "ll"  = {compute_aapc_linear(fit[[indx]]$ll[[1]])},
                     "sg1" = {compute_aapc_interval(fit[[indx]]$sg1[[1]], to_year=to_year, from_year=from_year)},
                     "sg2" = {compute_aapc_interval(fit[[indx]]$sg2[[1]], to_year=to_year, from_year=from_year)}
    )
    
    grp_indx = which(df$group==names(groups)[i])
    
    df$"Joint-Point (Year)"[grp_indx]     = jps[1]
    df$"Joint-Point (Std. Err)"[grp_indx] = jps[2]
    df$"AAPC (Est)"[grp_indx]             = result[1]
    df$"AAPC (Std. Err)"[grp_indx]        = result[2]
    df$"CI(95%).l"[grp_indx]              = result[3]
    df$"CI(95%).u"[grp_indx]              = result[4]
    
    
  }
  
  df$group = names(age_bands)
  
  return(df)
  
}



## Compute AAPC Comparison 
aapc_comparison <- function(x, y, from_year, to_year, conf.level=0.95){
  
  ## x = sg1[[1]]
  ## y = sg1[[2]]
  
  A = compute_aapc_interval(x, from_year=from_year, to_year=to_year, conf.level = conf.level, compare=TRUE)
  B = compute_aapc_interval(y, from_year=from_year, to_year=to_year, conf.level = conf.level, compare=TRUE)
  D = A[3] - B[3]
  
  p.val = 2*pnorm(q=D, lower.tail=FALSE)
  z     = abs(qnorm((1 - conf.level)/2))
  se.D  = sqrt(exp(2*A[1]) * A[4] + exp(2*B[1]) * B[4])
  
  
  r     <- c(D = D,  D + c(-z, z) * se.D, p.val)
  cin   <- paste("CI", "(", conf.level * 100, "%", ")", c(".l", ".u"), sep = "")
  names(r) <- c("D",cin, "p.val")
  
  return(r)
  
}


joint_point_model <- function(x, conf_level=0.95){
  dat  = list()
  ll   = list()
  sg1  = list()
  sg2  = list()
  sp   = list()
  
  ## x ==> diag_year
  ## y ==> log_est
  dat[[1]] = x
  
  ## Compute 0-2 joint points
  ll[[1]]  = lm( log_est ~ diag_year, data=dat[[1]])
  sg1[[1]] = segmented(ll[[1]], npsi = 1) 
  sg2[[1]] = segmented(ll[[1]], npsi = 2) 
  
  ## Spline fit
  sp[[1]]  = gam(log_est ~ s(diag_year, k=10, bs="cs"), data=dat[[1]], method="REML", select=TRUE)
  
  dat[[1]]$jp_none = predict(ll[[1]])
  dat[[1]]$jp_one  = predict(sg1[[1]])
  dat[[1]]$jp_two  = predict(sg2[[1]])
  dat[[1]]$spline  = predict(sp[[1]])

  
  ## Error Bands
  CI = setNames(vector(mode="list", length=3), c("none","one","two"))
  
  ## model goodness-of-fit estimations 
  BIC = setNames(vector(mode="list", length=3), c("none","one","two"))
  
  ci             = data.frame(diag_year=dat[[1]]$diag_year)
  ci$total_cases = as.numeric(dat[[1]]$total_cases)
  ci$total_pop   = as.numeric(dat[[1]]$total_pop)
  ci$group       = dat[[1]]$group
  
  ci_for_fit <- function(ci, pred, conf_level=0.95){
    
    z = qnorm(1-(1-0.95)/2)
    ci$est         = exp(pred$fit)
    ci$se_model    = exp(pred$se)
    ci$total_cases = as.numeric(dat[[1]]$total_cases)
    ci$se_poisson  = ci$est / sqrt(ci$total_cases)
    ci$se_comb     = sqrt(ci$se_model^2 + ci$se_poisson^2)
    # ci$total_pop   = as.numeric(dat[[1]]$total_pop)
    # ci$obs         = ci$total_cases / ci$total_pop * 1e5
    # ci$se_obs      = sqrt(ci$total_cases) / ci$total_pop * 1e5
    # ci$relative_se = ci$se_obs/ ci$obs
    # ci$se_pred     = ci$est * ci$relative_se
    # ci$scale       = ci$est / ci$obs
    # ci$ci_gamma_lw = if_else( ci$total_cases > 0, 
    #                           ci$scale * ci$obs * qgamma(0.025, shape=ci$total_cases, rate=1)/ci$total_cases,
    #                           0)
    # ci$ci_gamma_up = if_else( ci$total_cases > 0, 
    #                           ci$scale * ci$obs * qgamma(0.975, shape=ci$total_cases+1, rate=1)/ci$total_cases,
    #                           qgamma(0.975, shape=1, rate=1))
    ci$ci_norm_lw  = ci$est - z * ci$se_comb
    ci$ci_norm_up  = ci$est + z * ci$se_comb
    ci$ci_norm_lw  = if_else( ci$ci_norm_lw > 0, ci$ci_norm_lw, 0)
    
    return(ci)
  }
  
  
  partial_r2_two_segments <- function(fit, data){
    x <- data[["diag_year"]]
    y <- data[["log_est"]]
    
    tau <- sort(drop(fit$psi[,"Est."])) ## get joint-points
    
    k   <- length(tau)
    if( k == 0 ) return(numeric(0))
    
    ## tau extended
    ## tau_0 and tau_{k+1}
    tau_ext <- c(min(x)-1, tau, max(x))
    
    r2 <- as.numeric(k)
    
    for( i in seq_len(k) ){
      
      tau_i   <- tau_ext[i+1]
      left_b  <- tau_ext[i]
      right_b <- tau_ext[i+2]
      
      idx     <- (x > left_b) & (x <= right_b)
      subx    <- x[idx]
      suby    <- y[idx]
      
      z       <- pmax(subx - tau_i, 0)
      
      ## residualise on straight line (intercept + x) within this window
      ytil    <- residuals(lm(suby ~ subx))
      ztil    <- residuals(lm(   z ~ subx))
      
      if( sd(ztil) == 0 || sd(ytil) == 0 ){
        r2[i] <- NA
      } else {
        r2[i] <- cor(ytil, ztil)^2
      }
      
      return(r2)
    }
  }
  
 
  compute_bics <- function(x, data){
    ## bic  = -2*log(L)/n + p*log(n)/n
    ## bic3 = -2*log(L)/n + (3*k/n)*log(n)
    ## logL = log-likelihood of the fitted model
    ## p    = number of estimated parameters (including breakpoints)
    ## k    = number of breakpoints
    ## n    = number of observations
    ll0     <- logLik(x)
    k0      <- 0
    r2max   <- 0
    p0      <- attr(ll0, "df")
    n0      <- nobs(x)
    ll_term <- -2 * as.numeric(ll0) / n0
    logn    <- log(n0) / n0
    
    if( inherits(x, "segmented")){
      k0    <- nrow(x$psi)
      r2    <- suppressWarnings(partial_r2_two_segments(fit=x, data))
      r2max <- max(r2, na.rm=TRUE)
    }
    
    bic  <- ll_term +               p0 * logn
    bic3 <- ll_term +           3 * k0 * logn  
    wbic <- ll_term + (2 + r2max) * k0 * logn
    return(list(bic=bic, bic3=bic3, wbic=wbic))
  } 
  
  CI[[1]] = ci_for_fit(ci, pred=predict(ll[[1]], se.fit=TRUE),  conf_level = conf_level)
  CI[[2]] = ci_for_fit(ci, pred=predict(sg1[[1]], se.fit=TRUE), conf_level = conf_level)
  CI[[3]] = ci_for_fit(ci, pred=predict(sg2[[1]], se.fit=TRUE), conf_level = conf_level)
  
  BIC[[1]] = compute_bics(x=ll[[1]],  data=dat[[1]])
  BIC[[2]] = compute_bics(x=sg1[[1]], data=dat[[1]])
  BIC[[3]] = compute_bics(x=sg2[[1]], data=dat[[1]])
  
  pred           = predict(sg1[[1]], se.fit=TRUE)
  
  # ci             = data.frame(diag_year=dat[[1]]$diag_year)
  # ci$est         = exp(pred$fit)
  # ci$total_cases = as.numeric(dat[[1]]$total_cases)
  # ci$total_pop   = as.numeric(dat[[1]]$total_pop)
  # ci$se_obs      = sqrt(ci$total_cases) / ci$total_pop * 1e5
  # ci$ci_gamma_lw = if_else( ci$total_cases > 0, 
  #                           ci$est * qgamma(0.025, shape=ci$total_cases, rate=1)/ci$total_cases,
  #                           0)
  # ci$ci_gamma_up = if_else( ci$total_cases > 0, 
  #                           ci$est * qgamma(0.975, shape=ci$total_cases+1, rate=1)/ci$total_cases,
  #                           ci$est * qgamma(0.975, shape=1, rate=1))
  # ci$lower = ci$est - 1.96 * ci$se_obs #exp(pred$se)
  # ci$upper = ci$est + 1.96 * ci$se_obs #exp(pred$se)
  # ci$group = dat[[1]]$group
  
  return(list(dat=dat, ll=ll, sg1=sg1, sg2=sg2, sp=sp, pred=pred, ci=CI, bic=BIC))
  
}


##----------------------------------------------------------------------------
## Plot Functions
##----------------------------------------------------------------------------
base_incidence_plot <- function(x, est_var="asr", plot_type=1){
  
  y_axis_breaks <- c(1,2,3,4,5,7,10,20,30,50,70,100,200)

  est_name <- names(est_var)
  
  df <- x %>% 
        mutate( year = diag_year) %>%
        mutate( group = as.factor(group)) %>%
        mutate(y = as.numeric(!!sym(est_var))) %>%
        filter(y > 0 ) %>%
        dplyr::select(-!!sym(est_var))

  
 
  base_layers <- list(
    geom_point(size=2),
      geom_line(linewidth=1),
      scale_color_manual(values=c("early"="red", "late"="blue")),
      ## Year Bowel Cancer Screening was introduced in Scotland
      geom_vline(xintercept = 2007, color = "grey50", linetype = "dashed"),
      labs(
        title= sprintf("SMR06 (%s)", est_name),
        x="Diagnosis Year",
        y="Age-Standardised rate per 100,000"
      ),
    comp_axis_theme
  )
  
  if( plot_type == 1){
    gp <- ggplot(df, aes(x=year, y=y, color=group))+
          base_layers+
          scale_y_log10()
  } else {
  gp <- ggplot(df, aes(x=year, y=y, color=group))+
        base_layers +
    scale_y_log10(
      limits = range(y_axis_breaks),
      breaks = y_axis_breaks,
      labels = label_number(accuracy = 1)
    )+
    #annotation_logticks(sides="l")+
    theme(
      panel.grid.major.y = element_line(color="#EBEBEB"),
      panel.grid.minor.y = element_blank() )
  }

  return(gp)
  
}
  

base_incidence_plot2 <- function(x, est_var="asr", plot_type=1,
                                 y_axis_breaks = c(0.1,1,2,3,4,5,7,10,20,30,50,70,100,200),
                                 legend.labs = c("late"="late, > 50", "sub_early"="early, 40-49", "early"="early, 20-39")){
  
  est_name <- names(est_var)
  
  df <- x %>% 
    mutate( year  = diag_year) %>%
    mutate( group = as.factor(group)) %>%
    mutate(y = as.numeric(!!sym(est_var))) %>%
    filter(y > 0 ) %>%
    dplyr::select(-!!sym(est_var))
  
  
  base_layers <- list(
    geom_point(size=2),
    geom_line(linewidth=1),
    {if(is.na(legend.labs[which(names(legend.labs)=="sub_early")][[1]])){
      scale_color_manual(values=c("early"     = "#FF0000", 
                                  "late"      = "blue"),
                         labels=c("early"     = legend.labs[which(names(legend.labs)=="early")][[1]], 
                                  "late"      = legend.labs[which(names(legend.labs)=="late")][[1]]))
      } else {
      scale_color_manual(values=c("early"     = "#FF0000", 
                                  "sub_early" = "#8B0000", 
                                  "late"      = "blue"),
                         labels=c("early"     = legend.labs[which(names(legend.labs)=="early")][[1]], 
                                  "sub_early" = legend.labs[which(names(legend.labs)=="sub_early")][[1]], 
                                  "late"      = legend.labs[which(names(legend.labs)=="late")][[1]])) 
      }
    },
    labs(
      title= sprintf("SMR06 (%s)", est_name),
      x="Diagnosis Year",
      y="Age-Standardised rate per 100,000"
    ),
    #summary_axis_theme2
    comp_axis_theme
  )
  
  if( plot_type == 1){
    gp <- ggplot(df, aes(x=year, y=y, color=group))+
      base_layers+
      scale_y_log10()
  } else {
    gp <- ggplot(df, aes(x=year, y=y, color=group))+
      base_layers +
      scale_y_log10(
        limits = range(y_axis_breaks),
        breaks = y_axis_breaks,
        labels = label_number(accuracy = 1)
      )+
      theme(
        panel.grid.major.y = element_line(color="#EBEBEB"),
        panel.grid.minor.y = element_blank() )
  }
  
  return(gp)
  
}
