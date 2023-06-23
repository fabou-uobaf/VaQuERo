

##################
#generic helper functions

`%notin%` <- Negate(`%in%`)




##################
#statistics helper functions

# function to calculate shape2 pamerged_samplerameter of a betadistribution, given shape1 and expected value
betaParamFromMean <- function(mean, shape1){
    shape2 <- (shape1 - mean*shape1)/mean
    return(shape2)
}


##################
#time/date helper functions

leapYear <- function(y){
  if((y %% 4) == 0) {
    if((y %% 100) == 0) {
      if((y %% 400) == 0) {
        return(366)
      } else {
        return(365)
      }
    } else {
      return(366)
    }
  } else {
    return(365)
  }
}
decimalDate <- function(x, d){
  strftime(as.POSIXct(as.Date(x)+(d/96)), format="%Y-%m-%d %H:%M", tz="UCT") -> dx
  daysinyear <- leapYear(year(dx))
  year(dx) + yday(dx)/daysinyear + hour(dx)/(daysinyear*24) + minute(dx)/(daysinyear*24*60) -> ddx
  return(ddx)
}


##################
#optim optimisation

objfnct <- function(data, par) {
  varN <- length(par)
  rs <-rowSums(as.matrix(data[,2:(varN+1)]) * matrix(rep(par, each = dim(data)[1]), ncol = dim(data)[2]-2))
  re <- (rs - data$fit1)^2
  #re <- (rs - data$value.freq)^2
  #re <- (rs - (data$fit1 + data$value.freq)/2)^2
  sum(re)
}

starterV <- function(data){
  data <- as.data.frame(data)
  data[,2:(dim(data)[2])] -> data
  if(dim(data)[2]-1 > 1){
    data[which(rowSums(data[,1:(dim(data)[2]-1)]) == 1),] -> data
  }

  rep(NA, length(data)-1) -> startV

  for (i in seq_along(startV)){
    startV[i] <-  mean(data[which(data[,i] == 1),]$fit1)
  }
  (1-sum(startV, na.rm = TRUE))/sum(is.na(startV)) -> prop.rest
  startV[is.na(startV)] <- prop.rest
  return(startV)
}


##################
#NLoptr inequality constrainted optimisation

# Objective function
eval_f <- function(par, data){
        varN <- length(par)
        rs <-rowSums(as.matrix(data[,2:(varN+1)]) * matrix(rep(par, each = dim(data)[1]), ncol = dim(data)[2]-2))
        #re <- (data$fit1 + data$value.freq)/2
        #re <- data$value.freq
        re <- data$fit1
        ev <- (rs-re)^2
        return(sum(ev))
}


# Inequality constraints
eval_g_ineq <- function (x, data) {
        constr <- c(sum(x)-1)
        return (constr)
}

# Equality constraints
eval_g_eq <- function (x, data) {
        constr <- c(sum(x)-1)
        return (constr)
}

opts <- list(
            "algorithm" = "NLOPT_LN_COBYLA",
            "maxeval"   = 2000,
            "tol_constraints_ineq" = 0.0001
        )

##################
#Variant name conversion


## dealias variants
dealias <- function(x){
  base <- strsplit(x, split="\\.")[[1]][1]

  if(!any(names(aliases) == base)){
    y <- base
  } else if(nchar(aliases[names(aliases) == base]) == 0){
    y <- base
  } else if(grepl("^X", base)){
    y <- base
  } else {
    y <- aliases[names(aliases) == base]
  }
  dealiased <- gsub(base, y, x)
  return(dealiased)
}

## re-alias dealiased variant names
realias <- function(x){
  seq <- strsplit(x, split="\\.")[[1]]
  seq <- seq[-length(seq)]

  if(length(seq) > 0){
    for (i in rev(seq_along(seq))){
      trunc <- paste(seq[1:i], collapse = ".")
      if ( any(names(dealiases) == trunc) ){
        alias <- dealiases[names(dealiases) == trunc]
        realiased <- gsub(trunc, alias, x)
                  #print(paste("realiased:", x, " <=> ", realiased))
        break
      }
    }
  } else{
    realiased <- x
  }
  if(exists("realiased")){
    return(realiased)
  } else{
    return(x)
  }
}


# extract from ";"-separated list x all occurences of elements in y
extract_loi <- function(x, y ){
  llisty <- unlist(str_split(x, ";"))
  if(any(llisty %in% y)){
    treffer <- paste(llisty[llisty %in% y], collapse = ";", sep=";")
    restl   <- paste(llisty[llisty %notin% y], collapse = ";", sep=";")
    return(c(treffer))
  } else{
    return(c("NA"))
  }
}
extract_rest <- function(x, y ){
  llisty <- unlist(str_split(x, ";"))
  if(any(llisty %in% y)){
    treffer <- paste(llisty[llisty %in% y], collapse = ";", sep=";")
    restl   <- paste(llisty[llisty %notin% y], collapse = ";", sep=";")
    return(c(restl))
  } else{
    return(c(x))
  }
}


##################
#Plotting functions

getColor <- function(n, id, i){
  colorSet <- colorSets[id]
  cols <- brewer.pal(9, colorSet)
  col.palette <- colorRampPalette(c(cols[4], cols[5]), space = "Lab")
  cols <- col.palette(n)
  cols[i]
}

# functions for sankey plot
expandVar <- function(var, L){
  l = length(unlist(strsplit(var, split = "\\.")))
  while(l < L){
    paste0(var, ".0") -> var
    l = length(unlist(strsplit(var, split = "\\.")))
  }
  return(var)
}
getTree <- function(x){
  x$variant_dealias -> vars
  x$long_variant -> long_var
  list() -> ll
  for ( i in seq_along(vars)){
    ancestor <- c()
    if(grepl("^B", vars[i]) & "B" %notin% vars) {
      ancestor <- c("B")
    }
    for ( j in seq_along(vars)){
      if(grepl(vars[j], long_var[i], fixed = TRUE) & (vars[i] != vars[j])){
        ancestor <- c(ancestor, vars[j])
      }
    }
    ll[[long_var[i]]] <- ancestor
  }
  return(ll)
}
makeObs <- function(x, n, ll, p){
  x %>% mutate(freq = floor(0.5+p*freq)) -> x
  Llevels <- paste0("level", 0:n)
  as.data.frame(matrix(ncol = length(Llevels))) -> res
  colnames(res) <- Llevels

  for (i in 1:dim(x)[1]){
    y <- x[i,]
    var <- y$long_variant
    freq <- y$freq
    ancestors <- sort(unlist(ll[names(ll) == var]))
    Llevel = length(ancestors)+1
    var1 <- gsub("(\\.0)*$", "", var)
    var1 <- realias(var1)
    ancestors <- unlist(lapply(as.list(ancestors), realias))
    z <- c(as.character(ancestors), rep(var1, length(Llevels) - length(ancestors) ))
    z <- as.data.frame(t(z))[rep(1, freq),]
    colnames(z) <- Llevels
    rbind(res, z) -> res
  }
  return(res)
}
# function to generate TeX table from data.frame


##################
#TeX

makeTexTab <- function(filename, TAB, legendTxt){
    tabheaddef = paste0("\\begin{longtable}{", paste(rep(paste0("p{", 0.8/length(colnames(TAB)), "\\textwidth}"), length(colnames(TAB))), collapse = " | "), "}")
    tabhead    = paste(paste(colnames(TAB), collapse = " & "), "\\\\")
    legend     = paste0("\\caption{", legendTxt, "}\\\\")

    tabheaddef = gsub("%", "\\%", tabheaddef, fixed = TRUE)
    tabheaddef = gsub("_", "\\_", tabheaddef, fixed = TRUE)
    tabhead = gsub("%", "\\%", tabhead, fixed = TRUE)
    tabhead = gsub("_", "\\_", tabhead, fixed = TRUE)
    legend = gsub("%", "\\%", legend, fixed = TRUE)
    legend = gsub("_", "\\_", legend, fixed = TRUE)

    write("\\begin{footnotesize}", file = filename, append = TRUE)
    write(tabheaddef, file = filename, append = TRUE)
    write(legend, file = filename, append = TRUE)
    write("\\hline", file = filename, append = TRUE)
    write(tabhead, file = filename, append = TRUE)
    write("\\hline", file = filename, append = TRUE)
    for (i in 1:dim(TAB)[1] ){
      line = paste(paste(TAB[i,], collapse = " & "), "\\\\")
      line = gsub("%", "\\%", line, fixed = TRUE)
      line = gsub("_", "\\_", line, fixed = TRUE)
      write(line, file = filename, append = TRUE)
    }
    write("\\hline", file = filename, append = TRUE)

    write("\\end{longtable}", file = filename, append = TRUE)
    write("\\label{tab:synopsis}", file = filename, append = TRUE)
    write("\\end{footnotesize}", file = filename, append = TRUE)
}



##################
#Function to detect lineages in mutation data based on (mutation) counts

detect_lineages <- function(DT_, timepoint_, verbose_ = FALSE){
    ## define variants for which at least minMarker and at least minMarkerRatio are found
    ## define variants which could be detected by unique markers
    ## define variants with at least minUniqMarkerRatio and minUniqMarker of all uniq markers are found
    ## remove all variants not detected based on minMarker and all which could be, but were not detected, based on uniq markers
    ## repeat above step until no change (or 10 times); by removing variants, marker can become unique for other variants
    ## define variants for which at least minUniqMarker and minUniqMarkerRatio are found in the reduced marker definition

    detectedLineages <- c()

    DT_ %>% filter(sample_date_decimal %in% timepoint_) %>%
          mutate(Variants = strsplit(as.character(Variants), ";")) %>%
          unnest(Variants) %>% group_by(Variants, ID, sample_date_decimal) %>%
          left_join(y = moi_marker_count, by = "Variants", multiple = "all") %>%
          ungroup() %>%
          filter(value.freq>zeroo) %>%
          filter(value.depth > min.depth) %>%
          group_by(Variants, ID, sample_date_decimal) %>%
          mutate(N=n()) %>%
          mutate(r=N/n) %>%
          filter(r >= minMarkerRatio) %>%
          filter(N >= minMarker) %>%
          summarize(.groups = "drop_last") %>%
          pull(Variants) %>%
          unique() -> variants_passed_markerCount_filter

    if(TRUE){
        DT_ %>% filter(sample_date_decimal %in% timepoint_) %>%
            filter(!grepl(";",Variants)) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            left_join(y = moi_uniq_marker_count, by = "Variants", multiple = "all") %>%
            ungroup() %>% filter(value.freq>zeroo) %>%
            filter(value.depth > min.depth) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            mutate(N=n()) %>%
            mutate(r=N/n) %>%
            filter(r > minUniqMarkerRatio) %>%
            filter(N >= minUniqMarker) %>%
            summarize(NUCS = paste(NUC, collapse = "; "), .groups = "drop_last") %>%
            ungroup() %>% dplyr::select(Variants, NUCS) %>%
            pull(Variants) %>% unique() -> variants_passed_uniqMarkerCount_filter

        if(verbose_){
            print(paste("LOG: variants detected based on unique markers as following:"))
            print(variants_passed_uniqMarkerCount_filter)
        }

        DT_ %>% filter(sample_date_decimal %in% timepoint_) %>%
            filter(!grepl(";",Variants)) %>%
            filter(value.depth > min.depth) %>%
            ungroup() %>%
            pull(Variants) %>% unique()  -> could_be_detected_based_on_uniqMarkers
        could_be_but_not_have_been_detected_based_on_uniqMarkers <- could_be_detected_based_on_uniqMarkers[could_be_detected_based_on_uniqMarkers %notin% variants_passed_uniqMarkerCount_filter]
        DT_ %>%
            filter(sample_date_decimal %in% timepoint_) %>%
            mutate(Variants = strsplit(as.character(Variants), ";")) %>%
            unnest(Variants) %>%
            filter(Variants %in% variants_passed_markerCount_filter) %>%
            filter(Variants %notin% could_be_but_not_have_been_detected_based_on_uniqMarkers) %>%
            group_by(ID, sample_date, sample_date_decimal, value.freq, value.depth, LocationID, LocationName, NUC) %>%
            summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>%
            ungroup() %>%
            dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC")-> DT_reduced
        DT_reduced %>%
            filter(sample_date_decimal %in% timepoint_) %>%
            filter(!grepl(";",Variants)) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            left_join(y = moi_uniq_marker_count, by = "Variants", multiple = "all") %>%
            ungroup() %>%
            filter(value.freq>zeroo) %>%
            filter(value.depth > min.depth) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            mutate(N=n()) %>%
            mutate(r=N/n) %>%
            filter(r > minUniqMarkerRatio) %>%
            filter(N >= minUniqMarker) %>%
            summarize(.groups = "drop_last") %>%
            pull(Variants) %>% unique() -> variants_passed_uniqMarkerCount_filter

    }
    detectedLineages <- unique(c(detectedLineages, variants_passed_uniqMarkerCount_filter[variants_passed_uniqMarkerCount_filter %in% variants_passed_markerCount_filter]))

    # repeat above variant detection per unique markers
    # removing each time variants which could have beed detected based on uniq markers but did not
    # make it 20 times or until no more change
    c=1
    list() -> variants_passed_uniqMarkerCount_filter.iterative
    variants_passed_uniqMarkerCount_filter.iterative[[c]] <- variants_passed_uniqMarkerCount_filter
    list() -> could_be_detected_based_on_uniqMarkers.iterative
    could_be_detected_based_on_uniqMarkers.iterative[[c]] <- could_be_detected_based_on_uniqMarkers
    while(c < 20){
        c <- c + 1
        DT_reduced %>%
            filter(sample_date_decimal %in% timepoint_) %>%
            filter(!grepl(";",Variants)) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            left_join(y = moi_uniq_marker_count, by = "Variants", multiple = "all") %>%
            ungroup() %>%
            filter(value.freq>zeroo) %>%
            filter(value.depth > min.depth) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            mutate(N=n()) %>% mutate(r=N/n) %>%
            filter(r > minUniqMarkerRatio) %>%
            filter(N >= minUniqMarker) %>%
            summarize(NUCS = paste(NUC, collapse = "; "), .groups = "drop_last") %>%
            ungroup() %>% dplyr::select(Variants, NUCS) %>% distinct() -> variants_passed_uniqMarkerCount_filter

        if(verbose_){
            print(paste("LOG: [", c, "] variants detected based on unique markers after excluding detectable-but-not-detected lineages as following:"))
            print(variants_passed_uniqMarkerCount_filter %>% filter(Variants %notin% variants_passed_uniqMarkerCount_filter))
        }

        variants_passed_uniqMarkerCount_filter %>%
            pull(Variants) %>% unique() -> variants_passed_uniqMarkerCount_filter_list

        DT_reduced %>%
            filter(sample_date_decimal %in% timepoint_) %>%
            filter(!grepl(";",Variants)) %>%
            filter(value.depth > min.depth) %>%
            ungroup() %>%
            pull(Variants) %>%
            unique() -> could_be_detected_based_on_uniqMarkers
        could_be_but_not_have_been_detected_based_on_uniqMarkers <- could_be_detected_based_on_uniqMarkers[could_be_detected_based_on_uniqMarkers %notin% variants_passed_uniqMarkerCount_filter_list]
        DT_reduced %>%
            filter(sample_date_decimal %in% timepoint_) %>%
            mutate(Variants = strsplit(as.character(Variants), ";")) %>%
            unnest(Variants) %>%
            filter(Variants %in% variants_passed_markerCount_filter) %>%
            filter(Variants %notin% could_be_but_not_have_been_detected_based_on_uniqMarkers) %>%
            group_by(ID, sample_date, sample_date_decimal, value.freq, value.depth, LocationID, LocationName, NUC) %>%
            summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>%
            ungroup() %>%
            dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC")-> DT_reduced
        DT_reduced %>%
            filter(sample_date_decimal %in% timepoint_) %>%
            filter(!grepl(";",Variants)) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            left_join(y = moi_uniq_marker_count, by = "Variants", multiple = "all") %>%
            ungroup() %>%
            filter(value.freq>zeroo) %>%
            filter(value.depth > min.depth) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            mutate(N=n()) %>%
            mutate(r=N/n) %>%
            filter(r > minUniqMarkerRatio) %>%
            filter(N >= minUniqMarker) %>%
            summarize(.groups = "drop_last") %>%
            pull(Variants) %>% unique() -> variants_passed_uniqMarkerCount_filter.iterative[[c]]
        could_be_detected_based_on_uniqMarkers.iterative[[c]] <- could_be_detected_based_on_uniqMarkers
        if(all(could_be_detected_based_on_uniqMarkers.iterative[[c]] %in% could_be_detected_based_on_uniqMarkers.iterative[[c-1]])){
            c = 100
        }
    }
    detectedLineages <- unique(c(detectedLineages, unique(unlist(variants_passed_uniqMarkerCount_filter.iterative))[unique(unlist(variants_passed_uniqMarkerCount_filter.iterative)) %in% variants_passed_markerCount_filter]))


    # call variant detected if AF of enough marker shared with a already detected variant Z can not be
    # explained by Z using a beta distribution
    beta_shape1 <- 2.2
    pval_th <- 0.05

    c = 1
    overhang.added.lineages = c()
    while(c < 20){
        c <- c + 1
        DT_reduced %>%
            filter(!grepl(";", Variants)) %>%
            group_by(Variants) %>%
            summarize(mean_overhang = mean(value.freq), .groups = "drop") -> expected_value_uniqSupported_lineages
        DT_reduced %>%
            rowwise() %>%
            filter(grepl(";", Variants)) %>%
            mutate(overhang = extract_loi( x = Variants, y = expected_value_uniqSupported_lineages$Variants), rest = extract_rest(x = Variants, y = expected_value_uniqSupported_lineages$Variants)) %>%
            filter(!is.na(overhang)) %>%
            filter(!grepl(";", rest)) %>%
            dplyr::select(Variants, overhang, rest, value.freq, NUC) %>%
            left_join(y = expected_value_uniqSupported_lineages, by = c("overhang" = "Variants")) %>%
            mutate(pval = pbeta(value.freq, beta_shape1, betaParamFromMean(mean_overhang, beta_shape1), ncp = 0, lower.tail = FALSE, log.p = FALSE)) %>%
            filter(pval <= pval_th) -> iteratively.added
        if(dim(iteratively.added)[1] > 0){
            if(verbose_){
                print(paste("LOG: lineages [rest] detected due to mutations [NUC] observed in excess [value.freq] to uniquely detected lineages [overhang]. only overhang would explain mean_overhang of AF"))
                print(iteratively.added %>% dplyr::select(rest, overhang, NUC, value.freq, mean_overhang, pval))
            }
        }
        iteratively.added %>% pull(rest) %>% unique() -> iteratively.overhang.added.lineages
        DT_reduced %>%
            mutate(Variants = strsplit(as.character(Variants), ";")) %>%
            unnest(Variants) %>%
            filter(Variants %in% variants_passed_markerCount_filter) %>%
            filter(Variants %notin% iteratively.overhang.added.lineages) %>%
            group_by(ID, sample_date, sample_date_decimal, value.freq, value.depth, LocationID, LocationName, NUC) %>%
            summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>%
            ungroup() %>%
            dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC") -> DT_reduced
        overhang.added.lineages <- c(overhang.added.lineages, iteratively.overhang.added.lineages)
        if(length(iteratively.overhang.added.lineages) == 0){
            c = 100
        }
    }
    detectedLineages <- unique(c(detectedLineages, overhang.added.lineages))


    ## accept all variants_passed_markerCount_filter if ancestors to accepted by variants_passed_uniqMarkerCount_filter
    detectedLineages_dealiased <- unlist(lapply(as.list(detectedLineages), dealias))
    variants_passed_markerCount_filter_dealiased     <- unlist(lapply(as.list(variants_passed_markerCount_filter), dealias))
    ancestors_to_add <- c()

    if(length(unique(variants_passed_markerCount_filter_dealiased)) > 1){
      for (ii in 1:length(variants_passed_markerCount_filter_dealiased)){
        if(variants_passed_markerCount_filter_dealiased[ii] %in% detectedLineages_dealiased){
          next
        } else if (any(grepl(variants_passed_markerCount_filter_dealiased[ii], detectedLineages_dealiased, fixed = TRUE))){
          ancestors_to_add <- unique(c(ancestors_to_add, realias(variants_passed_markerCount_filter_dealiased[ii])))
        }
      }
    }
    detectedLineages <- unique(c(detectedLineages, ancestors_to_add))

    return(detectedLineages)
}

##############################################
## function to detect lineage in mutation data by phylogenetic tree guided stepwise regression


getPhyloGuideTree <- function(x){
  x$variant         -> variant
  x$variant_dealias -> variant_dealias

  res.dt <- data.table(basis = character(), add1 = character(), status = character())

  for ( i in seq_along(variant_dealias)){
    descendant <- variant_dealias[i]
    potential_ancestors <- c()
    for ( j in seq_along(variant_dealias)){
      potential_ancestor <- variant_dealias[j]
      if ( potential_ancestor == descendant){
        next
      }
      if(grepl(paste0("^", potential_ancestor), descendant, fixed = FALSE, perl = TRUE)){
        potential_ancestors <- c(potential_ancestors, potential_ancestor)
      }
    }
    if(is.null(potential_ancestors)){
      res.dt <- rbind(res.dt, data.table(basis = "1", add1 = descendant, status = "to_be_considered"))
    } else{
      closest_ancestor <- data.table(potential_ancestors = potential_ancestors, descendant = descendant) %>%
          rowwise() %>% mutate(descendant_level = length(unlist(str_split(pattern="\\.", descendant)))) %>%
          mutate(ancestors_level = length(unlist(str_split(pattern="\\.", potential_ancestors)))) %>%
          mutate(diff_level = ancestors_level-descendant_level) %>%
          group_by(1) %>% filter(diff_level == max(diff_level))
      res.dt <- rbind(res.dt, data.table(basis = closest_ancestor$potential_ancestors, add1 = descendant, status = "to_be_considered"))
    }
  }

  return(res.dt)
}


drop1All.wrapper <- function(dt, all_vars, scope_vars){
    formula_full <- as.formula(paste0("value.freq ~ ", paste(all_vars, collapse = " + ")))
    mod0 <- gamlss(formula = formula_full, data = dt, family = "NO", trace = FALSE)
    droper <- drop1All(object = mod0, scope = scope_vars, print = TRUE, parallel = "no")
    return(droper)
}

treeguided_stepwise_regression_dataprep <- function(dt, dm, tree, time, verbose = 0, cores=1, q=.1, f=.02){

    ## remove data point not in current time point
    dt %>% filter(sample_date_decimal == time) -> dt

    ### transform to avoid 0|1
    dt %>% mutate(value.freq = asin(sqrt(value.freq))) -> dt

    # remove "OTHERS"
    dt %>% mutate(Variants = gsub("other;*", "", Variants)) -> dt

    # join desing matrix to data matrix
    left_join(x = dt, y = dm, by = c("NUC"), multiple = "all") -> dt

    ## remove mutations which are not marker of any of the tree$add1
    ## keep mutations which are marker of all of the tree$add1
    if( sum(colnames(dt) %in% tree$add1) > 1){
      dt[( rowSums(as.data.frame(dt)[,colnames(dt) %in% tree$add1 ]) > 0 & rowSums(as.data.frame(dt)[,colnames(dt) %in% tree$add1 ]) <= length(tree$add1)),] -> dt
    } else{
      dt[as.data.frame(dt)[,colnames(dt) %in% tree$add1 ] > 0,] -> dt
    }

    ## select only essential columns
    dt[,colnames(dt) %in% c(tree$add1, "value.freq")] -> dt

    return(dt)
}


treeguided_stepwise_regression_backward <- function(dt, dm, tree, time, verbose = 0, cores=1, q=.1, f=.02){
    cur_included <- c()
    cidx <- 1

    while( cidx < 100 & any(grepl("consider", tree$status)) ){

        tree %>% filter(status == "to_be_considered") %>% pull(basis) %>% unique() -> to_be_considered_prenodes

        tree <-  tree %>% rowwise() %>% mutate(status = ifelse( (status == "to_be_considered" & (add1 %notin% to_be_considered_prenodes | basis %in% cur_included)), "in_consideration", status))

        ## define regression formula
        scope_vars <- tree %>% filter(status == "in_consideration") %>% pull(add1) %>% unique()
        all_vars  <- tree %>% filter(status != "excluded") %>% pull(add1) %>% unique() ##!!

        dropMod <- drop1All.wrapper(dt = dt, all_vars = all_vars, scope_vars = scope_vars)

        mod0_aic  <- dropMod$AIC[1]
        dropMod_aic  <- dropMod$AIC[2:length(dropMod$AIC)]
        dropMod_pval <- p.adjust(dropMod$`Pr(Chi)`[2:length(dropMod$`Pr(Chi)`)], method="fdr")
        droppedVar <- data.table(var = scope_vars, aic = dropMod_aic, pval = dropMod_pval) %>% filter(aic < mod0_aic) %>% filter(pval > q)
        fixedVar  <- data.table(var = scope_vars, aic = dropMod_aic, pval = dropMod_pval) %>% filter(aic > mod0_aic) %>% filter(pval < q) %>% ungroup() %>% filter(aic == max(aic))

        tree <- tree %>% rowwise() %>% mutate(status = ifelse(add1 %in% droppedVar$var, "excluded", status))
        tree <- tree %>% rowwise() %>% mutate(status = ifelse(add1 %in% fixedVar$var, "included", status))

        if(verbose){
            left_join(x = tree, y = data.table(var = scope_vars, aic = dropMod_aic, pval = dropMod_pval), by = c("add1" = "var")) %>% mutate(pval_crit = ifelse(pval<q, "***", "")) %>% mutate(aic_crit = ifelse(aic > mod0_aic, "***", "")) -> log_overview
            print(paste("LOG [", cidx, "iteration ] | included: ", paste(fixedVar$var, collapse = ",")))
            print(paste("LOG [", cidx, "iteration ] | excluded: ", paste(droppedVar$var, collapse = ",")))
            print(paste("LOG [", cidx, "iteration ] | mod0_aic: ", mod0_aic))
            print(log_overview, n = 100)
        }

        tree <- tree %>% rowwise() %>% filter(status != "excluded")
        table(tree$status)

        if(length(fixedVar$var) == 0 & length(droppedVar$var) == 0){
          break
        }
        cidx <- cidx+1
    }

    cur_included <- tree %>% filter(status == "included") %>% pull(add1) %>% unique()

    formula_full <- as.formula(paste( "value.freq", "~", paste(cur_included, collapse=" + ")))
    for (method in c("NO")){
        mod1 <- tryCatch(gamlss(formula_full, data = dt, family = method, trace = FALSE),error=function(e) e, warning=function(w) w)
        if(!any(grepl("Error|Warning", mod1))){
            if(verbose){
               print(paste("LOG: gamlss formula_full:", method, "successful"))
            }
            break
        }
    }
    if(any(grepl("Error|Warning", mod1))){
      print(paste("WARNING: no successful regression with any of the tried methods: [", cidx, "iteration ]", methods))
      return(cur_included)
    }
    predicted_freqs_fullModel <- sin(coef(mod1)[names(coef(mod1)) != "(Intercept)"] + coef(mod1)["(Intercept)"])^2
    cur_included <- names(predicted_freqs_fullModel[predicted_freqs_fullModel > f])

    return(cur_included)

}

stepwise_regression <- function(data, design_matrix, voi, direction = "both", time, cores=1, verbose=0){

    methods <- c("NO")

    if(any("Variants" == colnames(data))){
        ## remove data point not in current time point
        data %>% filter(sample_date_decimal == time) -> data

        ### transform to avoid 1
        data %>% mutate(value.freq = asin(sqrt(value.freq))) -> data

        # remove "OTHERS"
        data %>% mutate(Variants = gsub("other;*", "", Variants)) -> data

        # join desing matrix to data matrix
        left_join(x = data, y = design_matrix, by = c("NUC"), multiple = "all") -> data

        ## remove mutations which are not marker of any of the voi
        ## keep mutations which are marker of all of the voi
        if( sum(colnames(data) %in% voi) > 1){
          data[( rowSums(as.data.frame(data)[,colnames(data) %in% voi ]) > 0 & rowSums(as.data.frame(data)[,colnames(data) %in% voi ]) <= length(voi)),] -> data
        } else{
          data[as.data.frame(data)[,colnames(data) %in% voi ] > 0,] -> data
        }

        ## select only essential columns
        data[,colnames(data) %in% c(voi, "value.freq")] -> data
    }

    ## define basic formula
    if(direction == "forward"){
      formula_basis <- as.formula(paste("value.freq", "~", "1"))
    } else{
      formula_basis <- as.formula(paste("value.freq", "~", paste(unique(colnames(data)[colnames(data) != "value.freq"]), collapse = " + ")))
    }

    ## run basic model
    for (method in methods){
        mod0 <- tryCatch(gamlss(formula = formula_basis, data = data, family = method, trace = FALSE),error=function(e) e, warning=function(w) w)
        if(!any(grepl("Error|Warning", mod0))){
            if(verbose){
               print(paste("LOG: gamlss formula_basis: <", method, "> successful"))
            }
            break
        }
    }

    ## run stepwise regression
    scope_lower <- as.formula(paste("~", 1))
    scope_upper <- as.formula(paste("~", paste(unique(colnames(data)[colnames(data) != "value.freq"]), collapse = " + ")))
    mod1.log <- capture.output(tryCatch(mod1 <- stepGAICAll.B(mod0, scope=list(lower=scope_lower, upper=scope_upper), what="mu", direction = direction, trace = FALSE),error=function(e) e, warning=function(w) w))

    if(exists("mod1.log") & any(grepl("Model with term", mod1.log))){
      if(any(grepl(paste(voi, collapse = "|"), mod1.log))){
          ignored_terms <- c()
          for(voi_ in voi){
             if(any(grepl(paste("", voi_, ""), mod1.log, fixed = TRUE))){
                ignored_terms <- c(ignored_terms, voi_)
             }
          }
          voi <- voi[voi %notin% ignored_terms]
          data[,colnames(data) %in% c(voi, "value.freq")] -> data
          design_matrix[,colnames(data) %in% c(voi, "NUC")] -> design_matrix

          print(paste("WARNING: gamlss stepGAICAll.B:", method, "failed; retry without following term(s): ", paste(collapse=", ", ignored_terms)))

          included <- stepwise_regression(data = data, design_matrix = design_matrix, voi = voi, direction = direction, time = time, cores=cores, verbose=verbose)
          return(included)
      }
    }


    if(exists("mod1.log") & any(grepl("Error|Warning", mod1.log))){
      if(verbose){
        print(paste("LOG: gamlss stepGAICAll.B:", method, "failed; return NA"))
      }
      return(NA)
    }

    if(!exists("mod1")){
      print(paste("WARNING: gamlss stepGAICAll.B:", method, "failed; return NA"))
      return(NA)
    }


    ## deduce coefficients and select coefficients > 0 z-score
    deduded_freqs <- coef(mod1)[names(coef(mod1)) != "(Intercept)"] + coef(mod1)["(Intercept)"]
    deduded_freqs_zscore <- (deduded_freqs - mean(deduded_freqs))/sd(deduded_freqs)
    included_freqs <- deduded_freqs[deduded_freqs_zscore > .0]
    included <- names(included_freqs)

    return(included)

}
