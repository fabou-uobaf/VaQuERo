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
extract_loi <- function(x, y = expected_value_uniqSupported_lineages$Variants){
  llisty <- unlist(str_split(x, ";"))
  if(any(llisty %in% y)){
    treffer <- paste(llisty[llisty %in% y], collapse = ";", sep=";")
    restl   <- paste(llisty[llisty %notin% y], collapse = ";", sep=";")
    return(c(treffer))
  } else{
    return(c("NA"))
  }
}
extract_rest <- function(x, y = expected_value_uniqSupported_lineages$Variants){
  llisty <- unlist(str_split(x, ";"))
  if(any(llisty %in% y)){
    treffer <- paste(llisty[llisty %in% y], collapse = ";", sep=";")
    restl   <- paste(llisty[llisty %notin% y], collapse = ";", sep=";")
    return(c(restl))
  } else{
    return(c(x))
  }
}
