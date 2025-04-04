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
  if(is.factor(x)){
    x <- as.character(x)
  }
  ll <- strsplit(x, split="\\.")
  if(length(ll) >= 1){
    if(length(ll[[1]]) >= 1){
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
  } else{
    return(NA)
  }
}
# objectiv function for optimization
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
                  #writeLines(paste("realiased:", x, " <=> ", realiased))
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
    ancestor <- c("0")
    if(grepl("^B\\.", long_var[i]) & "A" %notin% vars) {
      ancestor <- c(ancestor, "B")
    }
    if(grepl("^A\\.", long_var[i]) & "A" %notin% vars) {
      ancestor <- c(ancestor, "A")
    }
    for ( j in seq_along(vars)){
      if(grepl(paste0(vars[j],"."), long_var[i], fixed = TRUE) & (vars[i] != vars[j])){
        ancestor <- c(ancestor, vars[j])
      }
    }
    ll[[long_var[i]]] <- unique(ancestor)
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
    tabheaddef = paste0("\\begin{longtable}{", paste(rep(paste0("p{", 0.90/length(colnames(TAB)), "\\textwidth}"), length(colnames(TAB))), collapse = " | "), "}")
    tabhead    = paste(paste(colnames(TAB), collapse = " & "), "\\\\")
    legend     = paste0("\\caption{", legendTxt, "}")

    tabheaddef = gsub("%", "\\%", tabheaddef, fixed = TRUE)
    tabheaddef = gsub("_", "\\_", tabheaddef, fixed = TRUE)
    tabhead = gsub("%", "\\%", tabhead, fixed = TRUE)
    tabhead = gsub("_", "\\_", tabhead, fixed = TRUE)
    legend = gsub("%", "\\%", legend, fixed = TRUE)
    legend = gsub("_", "\\_", legend, fixed = TRUE)

    write("\\begin{footnotesize}", file = filename, append = FALSE)
    write(tabheaddef, file = filename, append = TRUE)
    write("\\captionsetup{font=normalsize}", file = filename, append = TRUE)
    write("\\hline", file = filename, append = TRUE)
    write(tabhead, file = filename, append = TRUE)
    write("\\hline", file = filename, append = TRUE)
    for (i in 1:dim(TAB)[1] ){
      line = paste(do.call(paste, c(TAB[i,], collapse = " & ", sep = " & ")), "\\\\")
      line = gsub("%", "\\%", line, fixed = TRUE)
      line = gsub("_", "\\_", line, fixed = TRUE)
      write(line, file = filename, append = TRUE)
    }
    write("\\hline", file = filename, append = TRUE)
    write(legend, file = filename, append = TRUE)
    write("\\end{longtable}", file = filename, append = TRUE)
    write("\\label{tab:synopsis}", file = filename, append = TRUE)
    write("\\end{footnotesize}", file = filename, append = TRUE)
}

# function to calculate shape2 parameter of a betadistribution, given shape1 and expected value
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


detect_lineages <- function(DT_, timepoint_, minInfoRatio_){
    ## define variants for which at least minInfoRatio_ is found cummulatively over all detected marker
    ## define variants which could be detected by unique markers
    ## define variants with at least minUniqMarkerRatio and minUniqMarker of all uniq markers are found
    ## remove all variants not detected based on minMarker and all which could be, but were not detected, based on uniq markers
    ## repeat above step until no change (or 10 times); by removing variants, marker can become unique for other variants
    ## define variants for which at least minUniqMarker and minUniqMarkerRatio are found in the reduced marker definition
    ## define variants whos AF can not be explained by above detected markers assuming beta distribution (unless they could have been detected based on uniq markers)

    mutlog <- data.table(indicatingMutation = character(), ID = character(), LocationID = character(), sample_date = character(), sample_date_decimal = character(), Variants = character(), mutation = character(), meth = character(), iter = numeric())


    id_inuse <- DT_ %>% filter(sample_date_decimal %in% timepoint_) %>% pull(ID) %>% unique() %>% paste(collapse = ",")
    LocationID_inuse <- DT_ %>% filter(sample_date_decimal %in% timepoint_) %>% pull(LocationID) %>% unique() %>% paste(collapse = ",")
    sample_date_inuse <- DT_ %>% filter(sample_date_decimal %in% timepoint_) %>% pull(sample_date) %>% unique() %>% paste(collapse = ",") %>% as.character()
    sample_date_decimal_inuse <- DT_ %>% filter(sample_date_decimal %in% timepoint_) %>% mutate(sample_date_decimal = formatC(sample_date_decimal, digits = 7)) %>% pull(sample_date_decimal)  %>% unique() %>% paste(collapse = ",") %>% as.character()

    if(grepl(",", sample_date_decimal_inuse)){sample_date_decimal_inuse <- "multi"}
    if(grepl(",", sample_date_inuse)){sample_date_inuse <- "multi"}

    detectedLineages <- c()
    all_could_be_but_not_have_been_detected_based_on_uniqMarkers <- c()

    # collapse timepoints if several are specified
    DT_ <- DT_ %>% filter(sample_date_decimal %in% timepoint_) %>%
        group_by(Variants, LocationID, LocationName, NUC) %>%
        summarize(ID = DT_$ID[1], sample_date = DT_$sample_date[1], sample_date_decimal = timepoint_[1], value.freq = mean(value.freq, na.rm=TRUE), value.depth = mean(value.depth, na.rm = TRUE), .groups = "drop_last")


    # detect variants with enough information content given the detected marker mutations and minInfoRatio_
    # this is a necessary but not sufficient condition
    variants_markerInfo <- DT_ %>% filter(sample_date_decimal %in% timepoint_) %>%
          mutate(Variants = strsplit(as.character(Variants), ";")) %>%
          unnest(Variants) %>% group_by(Variants, ID, sample_date_decimal) %>%
          left_join(y = moi_marker_information, by = "NUC", multiple = "all") %>%
          group_by(Variants) %>%
          summarize(I = sum(i), .groups = "drop") %>%
          arrange(desc(I)) %>%
          left_join(y = moi_variant_information, by = "Variants", multiple = "all") %>%
          mutate(R = I/i)
    variants_passed_markerInfo_filter <- variants_markerInfo %>% filter(R >= minInfoRatio_) %>%
          pull(Variants) %>%
          unique()

    print(knitr::kable(variants_markerInfo %>% arrange(desc(R)) %>% mutate(Filter = ifelse(R > minInfoRatio_, "ACCEPTED due to Information", "! declined due to Information")) %>% filter(R > .75*minInfoRatio_), format = "rst", digits = 6, row.names = FALSE))


    if(length(variants_passed_markerInfo_filter) <= 0){
      return(detectedLineages)
    }

    # detect variants with enough uniq marker mutations
    variants_passed_uniqMarkerCount_filter.data <- DT_ %>% filter(sample_date_decimal %in% timepoint_) %>%
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
          ungroup() %>% dplyr::select(Variants, NUCS)
    variants_passed_uniqMarkerCount_filter <- variants_passed_uniqMarkerCount_filter.data %>% pull(Variants) %>% unique()

    #writeLines(paste("LOG: variants detected based on unique markers as following:"))
    #print(variants_passed_uniqMarkerCount_filter)


    if( length(variants_passed_uniqMarkerCount_filter) > 0) {
      mutlog <- rbind(mutlog, data.table(indicatingMutation = "indicatingMutation", ID = id_inuse, LocationID = LocationID_inuse, sample_date = sample_date_inuse, sample_date_decimal = sample_date_decimal_inuse, Variants = variants_passed_uniqMarkerCount_filter.data$Variants, mutation = variants_passed_uniqMarkerCount_filter.data$NUCS, meth = "UNIQ", iter = 0))
    }


    # detect variants with enough uniq marker mutations
    # after variants which could have been detect by uniq markers, but have not, were ignored
    could_be_detected_based_on_uniqMarkers <- moi %>% filter(!grepl(";", Variants)) %>% pull(Variants) %>% unique()
    could_be_but_not_have_been_detected_based_on_uniqMarkers <- could_be_detected_based_on_uniqMarkers[could_be_detected_based_on_uniqMarkers %notin% variants_passed_uniqMarkerCount_filter]
    DT_reduced <- DT_ %>%
          filter(sample_date_decimal %in% timepoint_) %>%
          mutate(Variants = strsplit(as.character(Variants), ";")) %>%
          unnest(Variants) %>%
          filter(Variants %in% variants_passed_markerInfo_filter) %>%
          filter(Variants %notin% could_be_but_not_have_been_detected_based_on_uniqMarkers) %>%
          group_by(ID, sample_date, sample_date_decimal, value.freq, value.depth, LocationID, LocationName, NUC) %>%
          summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>%
          ungroup() %>%
          dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC")
    variants_passed_uniqMarkerCount_iterative_filter.data <- DT_reduced %>%
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
          filter(N >= minUniqMarker)
    variants_passed_uniqMarkerCount_iterative_filter <- variants_passed_uniqMarkerCount_iterative_filter.data %>%
          summarize(.groups = "drop_last") %>%
          pull(Variants) %>% unique()
    all_could_be_but_not_have_been_detected_based_on_uniqMarkers <- unique(c(all_could_be_but_not_have_been_detected_based_on_uniqMarkers, could_be_but_not_have_been_detected_based_on_uniqMarkers))


    if( length(variants_passed_uniqMarkerCount_iterative_filter) > 0) {
      mutlog <- rbind(mutlog, data.table(indicatingMutation = "indicatingMutation", ID = id_inuse, LocationID = LocationID_inuse, sample_date = sample_date_inuse, sample_date_decimal = sample_date_decimal_inuse, Variants = variants_passed_uniqMarkerCount_iterative_filter.data$Variants, mutation = variants_passed_uniqMarkerCount_iterative_filter.data$NUC, meth = "iterUNIQ", iter = 1))
    }

    # merge above detected lineages
    detectedLineages <- unique(c(detectedLineages, variants_passed_uniqMarkerCount_filter[variants_passed_uniqMarkerCount_filter %in% variants_passed_markerInfo_filter], variants_passed_uniqMarkerCount_iterative_filter[variants_passed_uniqMarkerCount_iterative_filter %in% variants_passed_markerInfo_filter]))

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
        variants_passed_uniqMarkerCount_filter <- DT_reduced %>%
            filter(sample_date_decimal %in% timepoint_) %>%
            filter(!grepl(";",Variants)) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            left_join(y = moi_uniq_marker_count, by = "Variants", multiple = "all") %>%
            ungroup() %>%
            filter(value.freq>zeroo) %>%
            filter(value.depth > min.depth) %>%
            group_by(Variants, ID, sample_date_decimal) %>%
            mutate(N=n()) %>% mutate(r=N/n) %>%
            filter(N >= minUniqMarker)  %>%
            summarize(NUCS = paste(NUC, collapse = "; "), .groups = "drop_last") %>%
            ungroup() %>% dplyr::select(Variants, NUCS) %>% distinct()

        #writeLines(paste("LOG: [", c, "] variants detected based on unique markers after excluding detectable-but-not-detected lineages as following:"))
        #print(variants_passed_uniqMarkerCount_filter %>% filter(Variants %notin% variants_passed_uniqMarkerCount_filter))


        variants_passed_uniqMarkerCount_filter_list <- variants_passed_uniqMarkerCount_filter %>%
            pull(Variants) %>% unique()

        could_be_detected_based_on_uniqMarkers <- DT_reduced %>%
            filter(sample_date_decimal %in% timepoint_) %>%
            filter(!grepl(";",Variants)) %>%
            filter(value.depth > min.depth) %>%
            ungroup() %>%
            pull(Variants) %>%
            unique()
        could_be_but_not_have_been_detected_based_on_uniqMarkers <- could_be_detected_based_on_uniqMarkers[could_be_detected_based_on_uniqMarkers %notin% variants_passed_uniqMarkerCount_filter_list]
        DT_reduced <- DT_reduced %>%
            filter(sample_date_decimal %in% timepoint_) %>%
            mutate(Variants = strsplit(as.character(Variants), ";")) %>%
            unnest(Variants) %>%
            filter(Variants %in% variants_passed_markerInfo_filter) %>%
            filter(Variants %notin% could_be_but_not_have_been_detected_based_on_uniqMarkers) %>%
            group_by(ID, sample_date, sample_date_decimal, value.freq, value.depth, LocationID, LocationName, NUC) %>%
            summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>%
            ungroup() %>%
            dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC")
        variants_passed_uniqMarkerCount_filter.iterative.data <- DT_reduced %>%
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
            filter(N >= minUniqMarker)
        variants_passed_uniqMarkerCount_filter.iterative[[c]] <- variants_passed_uniqMarkerCount_filter.iterative.data %>%
            summarize(.groups = "drop_last") %>%
            pull(Variants) %>% unique()
        could_be_detected_based_on_uniqMarkers.iterative[[c]] <- could_be_detected_based_on_uniqMarkers
        all_could_be_but_not_have_been_detected_based_on_uniqMarkers <- unique(c(all_could_be_but_not_have_been_detected_based_on_uniqMarkers, could_be_but_not_have_been_detected_based_on_uniqMarkers))

        if( length(variants_passed_uniqMarkerCount_filter.iterative[[c]]) > 0) {
          mutlog <- rbind(mutlog, data.table(indicatingMutation = "indicatingMutation", ID = id_inuse, LocationID = LocationID_inuse, sample_date = sample_date_inuse, sample_date_decimal = sample_date_decimal_inuse, Variants = variants_passed_uniqMarkerCount_filter.iterative.data$Variants, mutation = variants_passed_uniqMarkerCount_filter.iterative.data$NUC, meth = "iterativeUNIQ", iter = c))
        }


        if(all(could_be_detected_based_on_uniqMarkers.iterative[[c]] %in% could_be_detected_based_on_uniqMarkers.iterative[[c-1]])){
            break
        }
    }

    #writeLines(paste0("LOG: variants detected based on iterative unique markers as following (2-",c, "):"))
    #print(variants_passed_uniqMarkerCount_iterative_filter)

    detectedLineages <- unique(c(detectedLineages, unique(unlist(variants_passed_uniqMarkerCount_filter.iterative))[unique(unlist(variants_passed_uniqMarkerCount_filter.iterative)) %in% variants_passed_markerInfo_filter]))


    # call variant detected if AF of enough marker shared with a already detected variant Z can not be
    # explained by Z using a beta distribution
    pval_th <- 0.05

    c = 20
    overhang.added.lineages = c()
    while(c < 40){
        c <- c + 1
        expected_value_uniqSupported_lineages <- DT_reduced %>%
            filter(!grepl(";", Variants)) %>%
            group_by(Variants) %>%
            summarize(mean_overhang = mean(value.freq), .groups = "drop")
        iteratively.added <- DT_reduced %>%
            rowwise() %>%
            filter(grepl(";", Variants)) %>%
            mutate(overhang = extract_loi( x = Variants, y = expected_value_uniqSupported_lineages$Variants), rest = extract_rest(x = Variants, y = expected_value_uniqSupported_lineages$Variants)) %>%
            filter(!is.na(overhang)) %>%
            filter(!grepl(";", rest)) %>%
            dplyr::select(Variants, overhang, rest, value.freq, NUC) %>%
            left_join(y = expected_value_uniqSupported_lineages, by = c("overhang" = "Variants")) %>%
            mutate(pval = pbeta(value.freq, alphaprime, betaParamFromMean(mean_overhang, alphaprime), ncp = 0, lower.tail = FALSE, log.p = FALSE)) %>%
            filter(pval <= pval_th)
        if(dim(iteratively.added)[1] > 0){
                writeLines(paste("LOG: lineages [rest] detected due to mutations [NUC] observed in excess [value.freq] to uniquely detected lineages [overhang]. only overhang would explain mean_overhang of AF"))
                print(iteratively.added %>% dplyr::select(rest, overhang, NUC, value.freq, mean_overhang, pval))
        }
        if( dim(iteratively.added)[1] > 0 ) {
          mutlog <- rbind(mutlog, data.table(indicatingMutation = "indicatingMutation", ID = id_inuse, LocationID = LocationID_inuse, sample_date = sample_date_inuse, sample_date_decimal = sample_date_decimal_inuse, Variants = iteratively.added$rest, mutation = iteratively.added$NUC, meth = "iterOVERHANG", iter = c))
        }

        iteratively.overhang.added.lineages <- iteratively.added %>% pull(rest) %>% unique()
        DT_reduced <- DT_reduced %>%
            mutate(Variants = strsplit(as.character(Variants), ";")) %>%
            unnest(Variants) %>%
            filter(Variants %in% variants_passed_markerInfo_filter) %>%
            filter(Variants %notin% iteratively.overhang.added.lineages) %>%
            group_by(ID, sample_date, sample_date_decimal, value.freq, value.depth, LocationID, LocationName, NUC) %>%
            summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>%
            ungroup() %>%
            dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC")
        overhang.added.lineages <- c(overhang.added.lineages, iteratively.overhang.added.lineages)
        if(length(iteratively.overhang.added.lineages) == 0){
            break
        }
    }
    # remove overhang.added.lineages if they could have been detected by uniq markers
    overhang.added.lineages <- overhang.added.lineages[overhang.added.lineages %notin% all_could_be_but_not_have_been_detected_based_on_uniqMarkers]

    #writeLines(paste0("LOG: variants detected based on iterative overhang mutations as following (1-",c-20, "):"))
    #print(overhang.added.lineages)


    detectedLineages <- unique(c(detectedLineages, overhang.added.lineages))

    ## collaps in the DT_reduced all variants of one mutation to its least common ancestor
    ## accept this one if it is accepted by variants_passed_uniqMarkerCount_filter
    DT_simplified <- DT_reduced %>%
      filter(sample_date_decimal %in% timepoint_) %>%
      rowwise() %>%
      mutate(Variant_vector = strsplit(as.character(Variants), ";")) %>%
      mutate(Variants_LCA = get_LCA(unlist(Variant_vector))) %>%
      filter(Variants_LCA %in% variants_passed_markerInfo_filter) %>%
      dplyr::select("Variants", "Variants_LCA", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC")

    variants_passed_uniqLCAMarkerCount_filter.data <- DT_simplified %>%
      filter(sample_date_decimal %in% timepoint_) %>%
      filter(!grepl(";",Variants_LCA)) %>%
      filter(value.freq>zeroo) %>%
      filter(value.depth > min.depth) %>%
      group_by(Variants_LCA, ID, sample_date_decimal) %>%
      summarize(NUCS = paste(NUC, collapse = "; "), .groups = "drop_last") %>%
      ungroup() %>% dplyr::select(Variants_LCA, NUCS)
    variants_passed_uniqLCAMarkerCount_filter <- variants_passed_uniqLCAMarkerCount_filter.data %>% pull(Variants_LCA) %>% unique()
    variants_passed_uniqLCAMarkerCount_filter <- variants_passed_uniqLCAMarkerCount_filter[variants_passed_uniqLCAMarkerCount_filter %notin% detectedLineages]
    variants_passed_uniqLCAMarkerCount_filter.data <- variants_passed_uniqLCAMarkerCount_filter.data %>% filter(Variants_LCA %notin% detectedLineages)


    #writeLines(paste("LOG: variants detected based on unique markers as following:"))
    #print(variants_passed_uniqMarkerCount_filter)


    if( length(variants_passed_uniqLCAMarkerCount_filter) > 0) {
      print(knitr::kable(variants_passed_uniqLCAMarkerCount_filter.data %>% rowwise() %>% mutate(NUCS = ifelse(nchar(NUCS) > 25, paste0(substr(NUCS, 1, 23), ".."), NUCS)), format = "rst", digits = 6, row.names = FALSE))
      mutlog <- rbind(mutlog, data.table(indicatingMutation = "indicatingMutation", ID = id_inuse, LocationID = LocationID_inuse, sample_date = sample_date_inuse, sample_date_decimal = sample_date_decimal_inuse, Variants = variants_passed_uniqLCAMarkerCount_filter.data$Variants_LCA, mutation = ifelse(nchar(variants_passed_uniqLCAMarkerCount_filter.data$NUCS) > 25, paste0(substr(variants_passed_uniqLCAMarkerCount_filter.data$NUCS, 1, 23), ".."), variants_passed_uniqLCAMarkerCount_filter.data$NUCS), meth = "UNIQ_LCA", iter = 0))
      detectedLineages <- unique(c(detectedLineages, variants_passed_uniqLCAMarkerCount_filter))
    }




    ## accept all variants_passed_markerInfo_filter if ancestors to accepted by variants_passed_uniqMarkerCount_filter
    detectedLineages_dealiased <- unlist(lapply(as.list(detectedLineages), dealias))
    variants_passed_markerInfo_filter_dealiased     <- unlist(lapply(as.list(variants_passed_markerInfo_filter), dealias))
    ancestors_to_add <- c()

    if(length(unique(variants_passed_markerInfo_filter_dealiased)) > 1){
      for (ii in 1:length(variants_passed_markerInfo_filter_dealiased)){
        if(variants_passed_markerInfo_filter_dealiased[ii] %in% detectedLineages_dealiased){
          next
        } else if (any(grepl(variants_passed_markerInfo_filter_dealiased[ii], detectedLineages_dealiased, fixed = TRUE))){
          ancestors_to_add <- unique(c(ancestors_to_add, realias(variants_passed_markerInfo_filter_dealiased[ii])))
        }
      }
    }
    detectedLineages <- unique(c(detectedLineages, ancestors_to_add))

    if( length(ancestors_to_add) > 0 ) {
      mutlog <- rbind(mutlog, data.table(indicatingMutation = "indicatingMutation", ID = id_inuse, LocationID = LocationID_inuse, sample_date = sample_date_inuse, sample_date_decimal = sample_date_decimal_inuse, Variants = ancestors_to_add, mutation = NA, meth = "ANCESTOR", iter = c))
    }

    if(dim(mutlog)[1] > 0){
      mutlog <- mutlog %>% group_by(Variants, mutation) %>% mutate(minIter = min(iter)) %>% filter(iter == minIter) %>% dplyr::select(-"minIter")
      print(knitr::kable(mutlog, format = "rst", digits = 6, row.names = FALSE))
    }

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



get_LCA.bak <- function(z){
    lcas <- c()

    if(length(z) == 1){
      lcas = z;
    } else{
      for (i in seq_along(z)){
        for (j in seq_along(z)){
          lca <- get_pairwise_LCA(z[i], z[j])
          if(nchar(lca) > 0 & lca %notin% z){
            lcas <- unique(c(lcas, lca))
          }
        }
      }
    }
    return(lcas)
}


get_LCA <- function(z){
    lcas <- z

    if(length(z) <= 1){
      lcas = z;
    } else{
      while(length(lcas) > 1){
        for (i in seq_along(lcas)){
          for (j in seq_along(lcas)){
            if(j <= i){
              next;
            }
            if(lcas[i] != "" & lcas[j] != ""){
              lca <- get_pairwise_LCA(lcas[i], lcas[j])
              lcas[i] <- lca
              lcas[j] <- lca
            } else{
              lcas[i] <- ""
              lcas[j] <- ""
            }
          }
        }
        lcas <- unique(lcas)
        if(any(unique(lcas) == "")){
          lcas <- c("")
          last;
        }
      }
    }

    return(lcas)
}



# generate cov-spectrum query link for list of mutations
covspectrumLinkAdvanced <- function(x){
  x <- unlist(lapply(unique(x), dephaseNuc))
  s <- paste(unique(x), collapse = ", ")
  l <- length(x)
  s <- paste0("[", ceiling(l/2), "-of: ", s, "]")
  return(s)
}
covspectrumLinkSimple <- function(x){
  x <- unlist(lapply(unique(x), dephaseNuc))
  s <- paste(unique(x), collapse = "%2C")
  l <- length(x)
  s <- paste0("\\href{https://cov-spectrum.org/explore/Europe/AllSamples/Past6M/variants?nucMutations=", s, "&}{", paste(x, collapse = " ", sep = " "),"}")
  return(s)
}
dephaseNuc <- function(x){
  y <- strsplit(x, split = "")[[1]]
  p <- as.numeric(paste0(y[grep("\\d", y)], collapse = ""))
  ref <- grep("\\D", y[1:(ceiling(length(y)/2))], value = TRUE)
  alt <- grep("\\D", y[(ceiling(length(y)/2)):length(y)], value = TRUE)
  r <- c()
  for (c in seq_along(ref)){
    paste(ref[c], p+c-1, ifelse(is.na(alt[c]), "-", alt[c]), sep = "")
    r <- append(r, paste(ref[c], p+c-1, ifelse(is.na(alt[c]), "-", alt[c]), sep = ""))
  }
  return(r)
}

# z-score normalisation function
zscore_norm <- function(x) {
    (x-mean(x))/sd(x)
}

n_topest <- function(x, n, descending = TRUE) {
  #Sort the wages in descending (default) order
  x1 <- sort(x, decreasing = descending)
  return(x1[n])
}

resetfunction <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

# function which takes a date and returns the date of wednesday of that week
date2weekwednesdaydate <- function(d){
  wwd = as.Date(d, tryFormats = c("%Y-%m-%d")) + 3 - as.numeric(strftime(as.Date(d, tryFormats = c("%Y-%m-%d")), format = "%u"))
  return(wwd)
}

# function to get median-min of selected wwtp
fun_median_min <- function(M){
    dist(M %>% dplyr::select(dcpLatitude, dcpLongitude)) -> distm
    # get for each wwtp the smallest non-zero value
    # get the mean across all wwtp of above value
    vdist <- median(apply(as.matrix(distm), 1,  function(x) n_topest(x[x>0], 1, FALSE)))
    return(vdist)
}

# function to get observed median-min distance for set of wwtp
fun_observed_distance <- function(l, XXX){
    locations <- unlist(str_split(l, pattern=";"))
    XXX %>% filter(LocationID %in% locations) -> YYY
    vdist <- fun_median_min(YYY)
    return(vdist)
}

# function to get expected median-min distance distriubtion between NNN random wwtp
fun_expected_distance_distro <- function(NNN, XXX, KKK){
    expected_distance_distro <- rep(NA, KKK)
    for (kkk in seq_along(expected_distance_distro)){
       YYY <- XXX
       random_assignment <- sample(c(rep(1, NNN), rep(0, dim(YYY)[1]-NNN)))
       YYY$selection <- random_assignment
       YYY %>% filter(selection == 1) -> YYY
       vdist <- fun_median_min(YYY)
       expected_distance_distro[kkk] <- vdist
    }
    return(expected_distance_distro)
}
