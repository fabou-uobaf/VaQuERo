`%notin%` <- Negate(`%in%`)
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

## dealias variants
dealias <- function(x){
  ll <- strsplit(x, split="\\.")
  if(length(ll) >= 1 & length(ll[[1]]) >= 1){
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
  } else{
    return(NA)
  }
}
objfnct <- function(data, par) {
  varN <- length(par)
  rs <-rowSums(as.matrix(data[,2:(varN+1)]) * matrix(rep(par, each = dim(data)[1]), ncol = dim(data)[2]-2))
  re <- (rs - data$fit1)^2
  #re <- (rs - data$value.freq)^2
  re <- unique(re)
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
getColor <- function(n, id, i){
  colorSet <- colorSets[id]
  cols <- brewer.pal(9, colorSet)
  col.palette <- colorRampPalette(c(cols[4], cols[5]), space = "Lab")
  cols <- col.palette(n)
  cols[i]
}

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

# function to calculate shape2 pamerged_samplerameter of a betadistribution, given shape1 and expected value
betaParamFromMean <- function(mean, shape1){
    shape2 <- (shape1 - mean*shape1)/mean
    return(shape2)
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


# function to label a variant if it is already represented by an ancestor in a list of variants
collapse2mother <- function(x, xset){
    xxset <- unlist(strsplit(xset, split=";"))
    if(length(xxset) == 1){
      return(TRUE)
    }
    xxset <- xxset[xxset != x]
    pattern <- paste(xxset, collapse = "|")
    !grepl(pattern, x)
}


detect_lineages <- function(DT_, timepoint_){
    ## define variants for which at least minMarker and at least minMarkerRatio are found
    ## define variants which could be detected by unique markers
    ## define variants with at least minUniqMarkerRatio and minUniqMarker of all uniq markers are found
    ## remove all variants not detected based on minMarker and all which could be, but were not detected, based on uniq markers
    ## repeat above step until no change (or 10 times); by removing variants, marker can become unique for other variants
    ## define variants for which at least minUniqMarker and minUniqMarkerRatio are found in the reduced marker definition
    ## define variants whos AF can not be explained by above detected markers assuming beta distribution (unless they could have been detected based on uniq markers)

    detectedLineages <- c()
    all_could_be_but_not_have_been_detected_based_on_uniqMarkers <- c()

    # detect variants with enough marker mutations
    # this is a necessary but not sufficient condition
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

    if(length(variants_passed_markerCount_filter) <= 0){
      return(detectedLineages)
    }

    # detect least common ancestor of variant combination mutations
    DT_ %>% filter(sample_date_decimal %in% timepoint_) %>%
          group_by(Variants, ID, sample_date_decimal) %>%
          left_join(y = moi_markerCombination_count, by = "Variants", multiple = "all") %>%
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
          unique() -> variants_passed_markerCombinationCount_filter
    mothers_variants_passed_markerCombinationCount_filter <- c()
    if(length(variants_passed_markerCombinationCount_filter) > 0){
        data.table(Variants = variants_passed_markerCombinationCount_filter, line = 1:length(variants_passed_markerCombinationCount_filter)) %>%
              mutate(Variants = strsplit(as.character(Variants), ";")) %>%
              unnest(Variants) %>%
              rowwise() %>%
              mutate(dealias = dealias(Variants)) %>%
              group_by(line) %>%
              mutate(all_variants = paste(sep=";", dealias, collapse = ";")) %>%
              rowwise() %>%
              mutate(keep = collapse2mother(dealias, all_variants)) %>%
              filter(keep == TRUE) %>% pull(Variants) %>%
              unique() -> mothers_variants_passed_markerCombinationCount_filter
        mothers_variants_passed_markerCombinationCount_filter <- get_LCA(mothers_variants_passed_markerCombinationCount_filter)
        mothers_variants_passed_markerCombinationCount_filter <- mothers_variants_passed_markerCombinationCount_filter[mothers_variants_passed_markerCombinationCount_filter %in% moi_marker_count$Variants]
    }
    print(paste("LOG: variants detected based on lca marker combinations as following:"))
    print(mothers_variants_passed_markerCombinationCount_filter)


    # detect variants with enough uniq marker mutations
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

    print(paste("LOG: variants detected based on unique markers as following:"))
    print(variants_passed_uniqMarkerCount_filter)

    # detect variants with enough uniq marker mutations
    # after variants which could have been detect by uniq markers, but have not, were ignored
    #DT_ %>% filter(sample_date_decimal %in% timepoint_) %>%
    #      filter(!grepl(";",Variants)) %>%
    #      filter(value.depth > min.depth) %>%
    #      ungroup() %>%
    #      pull(Variants) %>% unique()  -> could_be_detected_based_on_uniqMarkers
    moi %>% filter(!grepl(";", Variants)) %>% pull(Variants) %>% unique() -> could_be_detected_based_on_uniqMarkers
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
          filter(N >= minUniqMarker) %>%
          summarize(.groups = "drop_last") %>%
          pull(Variants) %>% unique() -> variants_passed_uniqMarkerCount_iterative_filter
    all_could_be_but_not_have_been_detected_based_on_uniqMarkers <- unique(c(all_could_be_but_not_have_been_detected_based_on_uniqMarkers, could_be_but_not_have_been_detected_based_on_uniqMarkers))

    print(paste("LOG: variants detected based on iterative unique markers as following (1):"))
    print(variants_passed_uniqMarkerCount_iterative_filter)

    # merge above detected lineages
    detectedLineages <- unique(c(detectedLineages, variants_passed_uniqMarkerCount_filter[variants_passed_uniqMarkerCount_filter %in% variants_passed_markerCount_filter], variants_passed_uniqMarkerCount_iterative_filter[variants_passed_uniqMarkerCount_iterative_filter %in% variants_passed_markerCount_filter]))

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
            filter(N >= minUniqMarker) %>%
            summarize(NUCS = paste(NUC, collapse = "; "), .groups = "drop_last") %>%
            ungroup() %>% dplyr::select(Variants, NUCS) %>% distinct() -> variants_passed_uniqMarkerCount_filter

        print(paste("LOG: [", c, "] variants detected based on unique markers after excluding detectable-but-not-detected lineages as following:"))
        print(variants_passed_uniqMarkerCount_filter %>% filter(Variants %notin% variants_passed_uniqMarkerCount_filter))

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
            filter(N >= minUniqMarker) %>%
            summarize(.groups = "drop_last") %>%
            pull(Variants) %>% unique() -> variants_passed_uniqMarkerCount_filter.iterative[[c]]
        could_be_detected_based_on_uniqMarkers.iterative[[c]] <- could_be_detected_based_on_uniqMarkers
        all_could_be_but_not_have_been_detected_based_on_uniqMarkers <- unique(c(all_could_be_but_not_have_been_detected_based_on_uniqMarkers, could_be_but_not_have_been_detected_based_on_uniqMarkers))
        if(all(could_be_detected_based_on_uniqMarkers.iterative[[c]] %in% could_be_detected_based_on_uniqMarkers.iterative[[c-1]])){
            break
        }
    }

    print(paste0("LOG: variants detected based on iterative unique markers as following (2-",c, "):"))
    print(variants_passed_uniqMarkerCount_iterative_filter)

    detectedLineages <- unique(c(detectedLineages, mothers_variants_passed_markerCombinationCount_filter, unique(unlist(variants_passed_uniqMarkerCount_filter.iterative))[unique(unlist(variants_passed_uniqMarkerCount_filter.iterative)) %in% variants_passed_markerCount_filter]))


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
                print(paste("LOG: lineages [rest] detected due to mutations [NUC] observed in excess [value.freq] to uniquely detected lineages [overhang]. only overhang would explain mean_overhang of AF"))
                print(iteratively.added %>% dplyr::select(rest, overhang, NUC, value.freq, mean_overhang, pval))
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
            break
        }
    }
    # remove overhang.added.lineages if they could have been detected by uniq markers
    overhang.added.lineages <- overhang.added.lineages[overhang.added.lineages %notin% all_could_be_but_not_have_been_detected_based_on_uniqMarkers]

    print(paste0("LOG: variants detected based on iterative overhang mutations as following (1-",c, "):"))
    print(overhang.added.lineages)


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

get_pairwise_LCA <- function(x,y){
   X <- unlist(str_split(pattern = "\\.", dealias(x)))
   Y <- unlist(str_split(pattern = "\\.", dealias(y)))
   Z <- c()
   for (i in 1:min(length(X), length(Y))){
     if(X[i] == Y[i]){
       Z <- append(Z, X[i])
     } else{
       break
     }
   }
   z <- realias(paste(Z, collapse = "."))
   return(z)
 }



get_LCA <- function(z){
    lcas <- c()
    for (i in seq_along(z)){
      for (j in seq_along(z)){
        lca <- get_pairwise_LCA(z[i], z[j])
        if(nchar(lca) > 0 & lca %notin% z){
          lcas <- unique(c(lcas, lca))
        }
      }
    }
    return(lcas)
}
