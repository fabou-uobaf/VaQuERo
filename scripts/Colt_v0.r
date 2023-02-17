## load libraries
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("gamlss"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("betareg"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("rjson"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggsankey"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("patchwork"))

timestamp()

# get Options
option_list = list(
  make_option(c("--dir"), type="character", default="ExampleOutput",
              help="Directory to write results [default %default]", metavar="character"),
  make_option(c("--metadata"), type="character", default="data/metaDataSub.tsv",
              help="Path to meta data input file [default %default]", metavar="character"),
  make_option(c("--marker"), type="character", default="resources/mutations_list.csv",
              help="Path to marker mutation input file [default %default]", metavar="character"),
  make_option(c("--mutstats"), type="character", default="resources/mutations_stats_pango.csv",
              help="Path to mutation statisitcs input file [default %default]", metavar="character"),
  make_option(c("--group2var"), type="character", default="resources/groupMembers.csv",
              help="Path to file linking group names (several variants) to exact pango variant IDs [default %default]", metavar="character"),
  make_option(c("--pmarker"), type="character", default="resources/mutations_problematic_all.csv",
              help="Path to problematic mutation input file, which will be ignored throughout the analysis [default %default]", metavar="character"),
  make_option(c("--data2"), type="character", default="data/mutationDataSub_sparse.tsv",
              help="Path to data input file in deprecated sparse matrix format [default %default]", metavar="character"),
  make_option(c("--data"), type="character", default="data/mutationDataSub.tsv.g",
              help="Path to data input file in tidy table format [default %default]", metavar="character"),
  make_option(c("--inputformat"), type="character", default="tidy",
              help="Input data format. Should be 'tidy' or 'sparse' [default %default]", metavar="character"),
  make_option(c("--detectmode"), type="character", default="umm",
              help="Mode to detect variants. Select from umm (unique marker mutation), cons (one consenus sequence) [default %default]", metavar="character"),
  make_option(c("--plotwidth"), type="double", default=8,
              help="Base size of plot width [default %default]", metavar="character"),
  make_option(c("--plotheight"), type="double", default=4.5,
              help="Base size of plot height [default %default]", metavar="character"),
  make_option(c("--ninconsens"), type="double", default="0.4",
              help="Minimal fraction of genome covered by reads to be considered (0-1) [default %default]", metavar="character"),
  make_option(c("--zero"), type="double", default=0.02,
              help="Minimal allele frequency to be considered [default %default]", metavar="double"),
  make_option(c("--depth"), type="integer", default=75,
              help="Minimal depth at mutation locus to be considered [default %default]", metavar="character"),
  make_option(c("--removeLongIndels"), type="logical", default=TRUE,
              help="Should indels with a length longer than 2 codons be ignored [default %default]", metavar="character"),
  make_option(c("--minuniqmark"), type="integer", default=1,
              help="Minimal absolute number of uniq markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--minuniqmarkfrac"), type="double", default=0.4,
              help="Minimal fraction of uniq markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--minqmark"), type="integer", default=3,
              help="Minimal absolute number of markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--minmarkfrac"), type="double", default=0.4,
              help="Minimal fraction of markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--colorBase"), type="character", default="B.1.617.2,BA.1,BA.2,BA.4,BA.5",
              help="List of variants whos decendences should be grouped by color in plots. Maximal 5 variants. List separated by semicolon [default %default]", metavar="character"),
  make_option(c("--debug"), type="logical", default="FALSE",
              help="Toggle to use run with provided example data [default %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#####################################################
####  parameter setting for interactive debugging
if(opt$debug){
  opt$dir = "sandbox1"
  opt$metadata = "data_wwtp/metaData_targeted.csv"
  opt$data="data_wwtp/mutationData_DB_TargetedMonitoringSites_phased.tsv.gz"
  opt$inputformat = "tidy"
  opt$marker="VaQuERo/resources/mutations_list_grouped_pango_codonPhased_2023-01-02_Europe.csv"
  opt$mutstats  = "VaQuERo/resources/mutations_stats_pango_codonPhased_2023-01-02.csv.gz"
  opt$group2var = "VaQuERo/resources/groupMembers_pango_codonPhased_2023-01-02_Europe.csv"
  opt$pmarker="VaQuERo/resources/mutations_problematic_vss1_v3.csv"
  opt$detectmode = "umm"
  opt$ninconsens = 0.1
  opt$zero=0.02
  opt$depth=5
  opt$minuniqmark=1
  opt$minuniqmarkfrac=0.4
  opt$minqmark=4
  opt$minmarkfrac=0.4
  opt$colorBase="B.1.617.2,BA.1,BA.2,BA.4,BA.5"
  print("Warning: command line option overwritten")
}
#####################################################

## define functions
options(warn=-1)
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

    write("\\label{tab:synopsis}", file = filename, append = TRUE)
    write("\\end{longtable}", file = filename, append = TRUE)
    write("\\end{footnotesize}", file = filename, append = TRUE)
}

# function to calculate shape2 pamerged_samplerameter of a betadistribution, given shape1 and expected value
betaParamFromMean <- function(mean, shape1){
    shape2 <- (shape1 - mean*shape1)/mean
    return(shape2)
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

# generate cov-spectrum query link for list of mutations
covspectrumLink <- function(x){
  x <- unlist(lapply(unique(x), dephaseNuc))
  s <- paste(unique(x), collapse = ", ")
  l <- length(x)
  s <- paste0("[", ceiling(l/2), "-of: ", s, "]")
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



## print parameter to Log

print("##~LOG~PARAMETERS~####################")
print(opt)
print("##~LOG~PARAMETERS~####################")
writeLines("\n\n\n")

## read in config file to overwrite all para below

## set variable parameters
summaryDataFile <- paste0(opt$dir, "/summary.csv")
markermutationFile  <- opt$marker
mutationstatsFile  <- opt$mutstats
group2variantFile  <- opt$group2var
problematicmutationFile <- opt$pmarker
sankeyPrecision <- 1000
plotWidth  <- opt$plotwidth
plotHeight <- opt$plotheight



## create directory to write plots
print(paste0("PROGRESS: create directory "))
timestamp()
outdir = opt$dir
if( ! dir.exists(outdir)){
  dir.create(outdir, showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs"))){
  dir.create(paste0(outdir, "/figs"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/sankey"))){
  dir.create(paste0(outdir, "/figs/sankey"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/excessmutation"))){
  dir.create(paste0(outdir, "/figs/excessmutation"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/scatter"))){
  dir.create(paste0(outdir, "/figs/scatter"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/overview"))){
  dir.create(paste0(outdir, "/figs/overview"), showWarnings = FALSE)
}




## read alias file
aliases <- fromJSON(file ="https://raw.githubusercontent.com/cov-lineages/pangolin-data/main/pangolin_data/data/alias_key.json")
as.list(names(aliases)) -> dealiases
names(dealiases) = as.character(aliases)

## set definition of failed/detected/passed
### how many none-N in consensus may be seen to be called detected
### genome_size - amplicon_covered_genome_size * 5%
N_in_Consensus_detection_filter <- 29903 - 29781 * 0.05
### how many non-N in consensus may be seen to be called "pass"
N_in_Consensus_filter <- 29903 - 29781 * opt$ninconsens

zeroo <- opt$zero      # marker below this value are considered to be zero
min.depth = opt$depth     # mutations with less read support are ignored
minUniqMarker <- opt$minuniqmark # minmal absolute number of uniq markers that variant is considered detected
minUniqMarkerRatio <- opt$minuniqmarkfrac # minmal fraction of uniq markers that variant is considered detected (0.25)
minMarker <- opt$minqmark # minmal absolute number of sensitive markers that variant is considered detected
minMarkerRatio <- opt$minmarkfrac # minmal fraction of sensitive markers that variant is considered detected (0.25)
baseForColorsVariants <- unlist(str_split(opt$colorBase, pattern=c(";",","))[[2]])
baseForColorsVariants <- unlist(lapply(as.list(baseForColorsVariants), dealias))
colorSets <- c("Blues", "Greens", "Purples", "Oranges", "Reds", "Greys")
removeLongIndels <- opt$removeLongIndels


## define global variables to fill
globalFittedData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_id = character(), sample_date = character(), value = numeric(), norm.value = numeric() )
globalFullData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_date = character(), value = numeric(), norm.value = numeric(), marker = character(), singlevalue = numeric() )
globalAFdata   <- data.table(sampleID = character(), LocationID = character(), LocationName  = character(), sample_date = character(), nuc_mutation = character(), aa_mutation = character(), observed = numeric(), expected = numeric(), excess = character(), pvalue = numeric() )
globalVarData   <- data.table(sampleID = character(), LocationID = character(), LocationName  = character(), sample_date = character(), variant = character(), deduced.freq = numeric() )

## read mutation stats
print(paste0("PROGRESS: read mutation per variant stats "))
mstat <- fread(file = mutationstatsFile)
group2var <- fread(file = group2variantFile)
left_join(x = mstat, y = group2var, by = c("ID" = "variants_alias")) %>% group_by(pos, ref, alt, AA_change, nucc, groupName) %>% summarize(value = sum(value, na.rm = TRUE), count = sum(count, na.rm = TRUE), .groups = "keep") %>% ungroup() %>% mutate(sensitivity = value/count) -> mstatgroup
if (removeLongIndels == TRUE){
  mstatgroup %>% filter(nchar(alt) <= 9 & nchar(ref) <= 9) -> mstatgroup
}

## read mutations of interest from file
print(paste0("PROGRESS: read mutation files "))
moi <- fread(file = markermutationFile)
unite(moi, NUC, c(4,3,5), ALT, sep = "", remove = FALSE) -> moi
moi %>% mutate(Variants = gsub("other;?", "", Variants)) %>% mutate(Variants = gsub(";;", ";", Variants)) -> moi   ### just there to fix GISAID issue with BQ.1.1 vs BQ.1
moi %>% dplyr::select("Variants") %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest() %>% group_by(Variants) %>% summarize(n = n()) -> moi_marker_count
moi %>% filter(!grepl(";", Variants)) %>% dplyr::select("Variants") %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest() %>% group_by(Variants) %>% summarize(n = n()) -> moi_uniq_marker_count

## read problematic mutations from file
poi <- fread(file = problematicmutationFile)
unite(poi, NUC, c(4,3,5), ALT, sep = "", remove = FALSE) -> poi
poi %>% dplyr::select("PrimerSet","NUC", "AA") -> problematic

## remove all problematic sites from soi and moi
moi %>% filter(NUC %notin% poi$NUC) -> moi

## define subset of marker and special mutations
moi %>% dplyr::select("Variants") %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest() %>% distinct() %>% filter(Variants != "other")-> variants_total
allmut   <- unique(c(moi$NUC))

## read in meta data
print(paste0("PROGRESS: read and process meta data "))
metaDT       <- fread(file = opt$metadata)
unique(metaDT$BSF_sample_name) -> sampleoi

## check if all voi can be detected in principle given used parameters
if (opt$detectmode == "umm"){
    moi %>% dplyr::select("Variants") %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest() %>% distinct() %>% filter(Variants != "other")-> variants_total
    if(TRUE){
        moi %>% filter(!is.na(Variants)) %>% dplyr::select(Variants, NUC) %>% distinct() %>% filter(!grepl(";",Variants))  %>% group_by(Variants) %>% summarise(nn = n()) %>% filter(nn >= minUniqMarker) %>% ungroup() %>% dplyr::select(Variants) %>% distinct() -> uniqMarkerPerVariant_could_be_detected
        moi %>% filter(!is.na(Variants)) %>% dplyr::select(Variants, NUC) %>% distinct() %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) %>% filter(Variants %notin% uniqMarkerPerVariant_could_be_detected$Variants) %>% group_by(NUC) %>% summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>% ungroup() -> moi_reduced
    }

    # repeat above bloc until no more change
    c=1
    list() -> uniqMarkerPerVariant_could_be_detected.iterative
    uniqMarkerPerVariant_could_be_detected.iterative[[c]] <- uniqMarkerPerVariant_could_be_detected$Variants
    while(c < 20){
        c <- c + 1
        moi_reduced %>% filter(!is.na(Variants)) %>% dplyr::select(Variants, NUC) %>% distinct() %>% filter(!grepl(";",Variants))  %>% group_by(Variants) %>% summarise(nn = n()) %>% filter(nn >= minUniqMarker) %>% ungroup() %>% dplyr::select(Variants) %>% distinct() -> uniqMarkerPerVariant_could_be_detected
        moi_reduced %>% filter(!is.na(Variants)) %>% dplyr::select(Variants, NUC) %>% distinct() %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) %>% filter(Variants %notin% uniqMarkerPerVariant_could_be_detected$Variants) %>% group_by(NUC) %>% summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>% ungroup() -> moi_reduced
        uniqMarkerPerVariant_could_be_detected.iterative[[c]] <- uniqMarkerPerVariant_could_be_detected$Variants
        if(all(uniqMarkerPerVariant_could_be_detected.iterative[[c]] %in% uniqMarkerPerVariant_could_be_detected.iterative[[c-1]])){
            c = 100
        }
    }
    variants_total$Variants[variants_total$Variants %in% sort(unique(unlist(uniqMarkerPerVariant_could_be_detected.iterative)))] -> variants_w_enough_uniq
    variants_total$Variants[variants_total$Variants %notin% variants_w_enough_uniq] -> variants_futile
    if(length(variants_futile) > 0){
      for (c in seq_along(variants_futile)){
        print(paste("Warning(s):", variants_futile[c], "has less than", minUniqMarker, "unique markers defined also considering descendant."))
      }
    }
    variants_total$Variants <- variants_w_enough_uniq
}

# check for same date, same location issue
# introduce artifical decimal date
metaDT %>% group_by(LocationID, sample_date) %>% mutate(n = rank(BSF_sample_name)-1)  %>% ungroup() %>% rowwise() %>% mutate(sample_date_decimal = decimalDate(sample_date, n)) %>% dplyr::select(-n) -> metaDT

# define which locations are used for the report
metaDT %>% dplyr::select("LocationID") %>% distinct() -> locationsReportedOn


## determine last sequencing batch ID
as.character(metaDT %>% filter(!is.na(BSF_run)) %>% filter(!is.na(BSF_start_date)) %>% group_by(BSF_run, BSF_start_date) %>% summarize(.groups = "drop") %>% ungroup() %>% filter(as.Date(BSF_start_date) == max(as.Date(BSF_start_date))) %>% dplyr::select("BSF_run")) -> last_BSF_run_id
print(paste("LOG: last_BSF_run_id", last_BSF_run_id))


## generate general model matrix
moi %>% dplyr::select("Variants","NUC") -> marker
marker %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) -> marker

matrix(rep(0,(length(unique(marker$Variants))*length(unique(marker$NUC)))), nrow = length(unique(marker$NUC))) -> mmat
as.data.table(mmat) -> mmat
colnames(mmat) <- unique(marker$Variants)
rownames(mmat) <- unique(marker$NUC)

for (i in 1:dim(marker)[1]){
  cc <- as.character(marker[i,1])
  rr <- as.character(marker[i,2])
  mmat[which(rownames(mmat) == rr), which(colnames(mmat) == cc)] <- 1
}
mmat$NUC <- rownames(mmat)


## read in mutations data
if(opt$inputformat == "sparse"){
  print(paste0("PROGRESS: read AF data (deprecated file format!) "))
  timestamp()
  sewage_samps <- fread(opt$data2 , header=TRUE, sep="\t" ,na.strings = ".", check.names=TRUE)

  ## remove mutation positions from input data which is not used further
  sewage_samps %>% filter(POS %in% unique(sort(c(moi$Postion)))) -> sewage_samps

  ## get sample names
  sample_names = grep("\\.AF$", names(sewage_samps), value = TRUE)
  sample_names = gsub("\\.AF$","",sample_names)

  ## generate nucc string
  unite(sewage_samps, NUC, c(3,2,4), ALT, sep = "", remove = FALSE) -> sewage_samps

  ## remove all positions which are problematic
  sewage_samps[sewage_samps$NUC %notin% poi$NUC,] -> sewage_samps

  ## set missing mutation frequencies to 0
  sewage_samps[is.na(sewage_samps)] <- 0

  ## tidy up sewage data table
  sewage_samps[,-c(10:18)] -> sewage_samps
  sewage_samps.dt <- reshape2::melt(sewage_samps, id.vars=1:10)
  rm(sewage_samps)
  sewage_samps.dt %>% separate( col = variable, into = c("ID", "VarType"), sep = "\\.") -> sewage_samps.dt
  sewage_samps.dt %>% filter(VarType == "AF") -> dt_AF
  sewage_samps.dt %>% filter(VarType == "DP") -> dt_DP
  inner_join(dt_AF, dt_DP, by =  c("CHROM", "NUC", "POS", "REF", "ALT", "ANN.GENE", "ANN.FEATUREID", "ANN.EFFECT", "ANN.AA", "EFF", "ID"), suffix = c(".freq", ".depth")) -> sewage_samps.dt

  ## remove all samples which are not specified in metadata
  sewage_samps.dt %>% filter(ID %in% sampleoi) -> sewage_samps.dt

} else{
  print(paste0("PROGRESS: read AF data "))
  timestamp()
  sewage_samps <- fread(opt$data , header=TRUE, sep="\t" ,na.strings = "NA", check.names=TRUE)
  colnames(sewage_samps)[colnames(sewage_samps) == "SAMPLEID"] <- "ID"
  colnames(sewage_samps)[colnames(sewage_samps) == "GENE"] <- "ANN.GENE"
  colnames(sewage_samps)[colnames(sewage_samps) == "AA"] <- "ANN.AA"
  colnames(sewage_samps)[colnames(sewage_samps) == "AF"] <- "value.freq"
  colnames(sewage_samps)[colnames(sewage_samps) == "DP"] <- "value.depth"
  colnames(sewage_samps)[colnames(sewage_samps) == "PQ"] <- "value.qual"

  sewage_samps %>% mutate(NUC = paste(REF, POS, ALT, sep = "")) -> sewage_samps.dt
  sewage_samps.dt$value.freq[is.na(sewage_samps.dt$value.freq)] <- 0
  sewage_samps.dt$value.depth[is.na(sewage_samps.dt$value.depth)] <- 0
  sewage_samps.dt$value.qual[is.na(sewage_samps.dt$value.qual)] <- 0

  ## remove all samples which are not specified in metadata
  sewage_samps.dt %>% filter(ID %in% sampleoi) -> sewage_samps.dt

  rm(sewage_samps)
}

#removeLongIndels
if (removeLongIndels == TRUE){
  sewage_samps.dt %>% filter(nchar(ALT) <= 9 & nchar(REF) <= 9) -> sewage_samps.dt
}

## add location to sewage_samps.dt
sewage_samps.dt %>% mutate(RNA_ID_int = gsub("_S\\d+$","", ID)) -> sewage_samps.dt
metaDT %>% filter(! is.na(LocationID) ) %>% dplyr::select("RNA_ID_int", "LocationID", "LocationName") -> sample_location
left_join(x = sewage_samps.dt, y = sample_location, by = "RNA_ID_int") -> sewage_samps.dt

## remove samples which are not include_in_report == TRUE
metaDT %>% filter(is.na(include_in_report) | include_in_report == TRUE ) %>% dplyr::select("RNA_ID_int", "BSF_sample_name") -> sample_includeInReport
sewage_samps.dt %>% filter(ID %in% sample_includeInReport$BSF_sample_name) -> sewage_samps.dt

## remove samples with >N_in_Consensus_filter N
metaDT  %>% filter(N_in_Consensus < N_in_Consensus_filter) %>% dplyr::select("RNA_ID_int", "BSF_sample_name") -> passed_samples
sewage_samps.dt <- sewage_samps.dt[sewage_samps.dt$ID %in% passed_samples$BSF_sample_name,]

## add variant of interests to table
left_join(x=sewage_samps.dt, y=moi, by = c("POS"="Postion", "REF", "ALT", "NUC")) -> sewage_samps.dt

## add sampling date
metaDT %>% dplyr::select("RNA_ID_int", "sample_date", "sample_date_decimal") -> sample_dates
left_join(x=sewage_samps.dt, y=sample_dates, by = "RNA_ID_int") -> sewage_samps.dt

## get most recent sampling date from last Run
metaDT %>% filter(BSF_start_date == sort(metaDT$BSF_start_date)[length(sort(metaDT$BSF_start_date))]) %>% dplyr::select("BSF_run", "RNA_ID_int", "BSF_sample_name", "sample_date") -> RNA_ID_int_currentRun

RNA_ID_int_currentRun %>% ungroup() %>% summarize(latest = max(as.Date(sample_date))) -> latestSample
RNA_ID_int_currentRun %>% ungroup() %>% summarize(earliest = min(as.Date(sample_date))) -> earliestSample

print(paste("LOG: current run ID:", unique(RNA_ID_int_currentRun$BSF_run)))
print(paste("LOG: earliest sample in current run:", earliestSample$earliest))
print(paste("LOG: latest sample in current run:", latestSample$latest))

### kill if no mutations are found
if (length(unique(sewage_samps.dt$ID)) < 1){
  print(paste("PROGRESS: start to loop over sample location"))
  quit(save="no")
}

### FROM HERE LOOP OVER EACH SAMPLE LOCATION
print(paste("PROGRESS: start to loop over sample location"))
timestamp()
for (s in 1:length(unique(sewage_samps.dt$ID))) {
    sample_identifier = unique(sewage_samps.dt$ID)[s]
    sewage_samps.dt %>% filter(ID == sample_identifier) -> sewage_samps.sdt

    for (r in 1:length(unique(sewage_samps.sdt$LocationID))) {
        roi = unique(sewage_samps.sdt$LocationID)[r]
        roiname = unique(sewage_samps.sdt$LocationName)[r]
        print(paste("PROGRESS:", roiname, roi))

        ### define data collector
        plantFittedData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_id = character(), sample_date = character(), value = numeric(), norm.value = numeric() )
        plantFullData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_date = character(), value = numeric(), norm.value = numeric(), marker = character(), singlevalue = numeric() )

        ### filter data set for current location
        sewage_samps.sdt %>% filter(LocationID == roi) %>% dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC") %>% filter(NUC %in% allmut) -> sdt

        ## count how many sample from this plant were in the last BSF run
        metaDT %>% filter(BSF_sample_name %in% sdt$ID) %>% filter(BSF_run %in% last_BSF_run_id) %>% dplyr::select(BSF_run) %>% group_by(BSF_run) %>% summarize(n = n()) -> count_last_BSF_run_id

        ### FROM HERE LOOP OVER TIME POINTS
        sdt %>% dplyr::select(sample_date_decimal, sample_date) %>% distinct() -> timePointsCombinations
        timePoints <- timePointsCombinations$sample_date_decimal[order(timePointsCombinations$sample_date_decimal)]
        timePoints_classic <- timePointsCombinations$sample_date[order(timePointsCombinations$sample_date_decimal)]

        for (t in 1:length(timePoints)) {
                timepoint <- timePoints[t]
                timepoint_classic <- timePoints_classic[t]
                timepoint_day <- decimalDate(timepoint_classic,0)
                T <- which(timePoints_classic == timepoint_classic)
                timepoint <- timePoints[T]

                sdt %>% filter(sample_date %in% timepoint_classic ) -> ssdt
                sampleID <- paste0(unique(ssdt$ID))

                print(paste("PROGRESS:", sampleID, "@", roiname, "@", timepoint_classic, "(", signif(timepoint, digits = 10), ")"))

                ## use different modes to detect variants to be tested
                if(opt$detectmode == "umm"){
                    ## define variants for which at least minMarker and at least minMarkerRatio are found
                    ## define variants which could be detected by unique markers
                    ## define variants with at least minUniqMarkerRatio and minUniqMarker of all uniq markers are found
                    ## remove all variants not detected based on minMarker and all which could be, but were not detected, based on uniq markers
                    ## repeat above step until no change (or 10 times); by removing variants, marker can become unique for other variants
                    ## define variants for which at least minUniqMarker and minUniqMarkerRatio are found in the reduced marker definition
                    ssdt %>% filter(sample_date_decimal == timepoint) %>%  mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) %>% group_by(Variants, ID, sample_date_decimal) %>% left_join(y = moi_marker_count, by = "Variants") %>% ungroup() %>% filter(value.freq>zeroo) %>% filter(value.depth > min.depth) %>% group_by(Variants, ID, sample_date_decimal) %>% mutate(N=n()) %>% mutate(r=N/n) %>% filter(r >= minMarkerRatio) %>% filter(N >= minMarker) %>% summarize(.groups = "drop_last") -> ssdt3
                    markerPerVariant <- unique(ssdt3$Variants)

                    if(TRUE){
                        ssdt %>% filter(sample_date_decimal %in% timepoint) %>% filter(!grepl(";",Variants)) %>% group_by(Variants, ID, sample_date_decimal) %>% left_join(y = moi_uniq_marker_count, by = "Variants") %>% ungroup() %>% filter(value.freq>zeroo) %>% filter(value.depth > min.depth) %>% group_by(Variants, ID, sample_date_decimal) %>% mutate(N=n()) %>% mutate(r=N/n) %>% filter(r > minUniqMarkerRatio) %>% filter(N >= minUniqMarker) %>% summarize(.groups = "drop_last") %>% ungroup() %>% dplyr::select(Variants) %>% distinct() -> uniqMarkerPerVariant_detected
                        ssdt %>% filter(sample_date_decimal %in% timepoint) %>% filter(!grepl(";",Variants)) %>% filter(value.depth > min.depth) %>% ungroup() %>% dplyr::select(Variants) %>% distinct()  -> uniqMarkerPerVariant_could_be_detected
                        uniqMarkerPerVariant_not_detected <- uniqMarkerPerVariant_could_be_detected$Variants[uniqMarkerPerVariant_could_be_detected$Variants %notin% uniqMarkerPerVariant_detected$Variants]
                        ssdt %>% filter(sample_date_decimal %in% timepoint) %>%  mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) %>% filter(Variants %in% markerPerVariant) %>% filter(Variants %notin% uniqMarkerPerVariant_not_detected) %>% group_by(ID, sample_date, sample_date_decimal, value.freq, value.depth, LocationID, LocationName, NUC) %>% summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>% ungroup() %>% dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC")-> ssdt_reduced
                        ssdt_reduced %>% filter(sample_date_decimal %in% timepoint) %>% filter(!grepl(";",Variants)) %>% group_by(Variants, ID, sample_date_decimal) %>% left_join(y = moi_uniq_marker_count, by = "Variants") %>% ungroup() %>% filter(value.freq>zeroo) %>% filter(value.depth > min.depth) %>% group_by(Variants, ID, sample_date_decimal) %>% mutate(N=n()) %>% mutate(r=N/n) %>% filter(r > minUniqMarkerRatio) %>% filter(N >= minUniqMarker) %>% summarize(.groups = "drop_last") -> ssdt2
                        uniqMarkerPerVariant <- unique(ssdt2$Variants)
                    }

                    # repeat above bloc until no more change
                    c=1
                    list() -> uniqMarkerPerVariant.iterative
                    uniqMarkerPerVariant.iterative[[c]] <- uniqMarkerPerVariant
                    list() -> uniqMarkerPerVariant_could_be_detected.iterative
                    uniqMarkerPerVariant_could_be_detected.iterative[[c]] <- uniqMarkerPerVariant_could_be_detected$Variants
                    while(c < 20){
                        c <- c + 1
                        ssdt_reduced %>% filter(sample_date_decimal %in% timepoint) %>% filter(!grepl(";",Variants)) %>% group_by(Variants, ID, sample_date_decimal) %>% left_join(y = moi_uniq_marker_count, by = "Variants") %>% ungroup() %>% filter(value.freq>zeroo) %>% filter(value.depth > min.depth) %>% group_by(Variants, ID, sample_date_decimal) %>% mutate(N=n()) %>% mutate(r=N/n) %>% filter(r > minUniqMarkerRatio) %>% filter(N >= minUniqMarker) %>% summarize(.groups = "drop_last") %>% ungroup() %>% dplyr::select(Variants) %>% distinct() -> uniqMarkerPerVariant_detected
                        ssdt_reduced %>% filter(sample_date_decimal %in% timepoint) %>% filter(!grepl(";",Variants)) %>% filter(value.depth > min.depth) %>% ungroup() %>% dplyr::select(Variants) %>% distinct()  -> uniqMarkerPerVariant_could_be_detected
                        uniqMarkerPerVariant_not_detected <- uniqMarkerPerVariant_could_be_detected$Variants[uniqMarkerPerVariant_could_be_detected$Variants %notin% uniqMarkerPerVariant_detected$Variants]
                        ssdt_reduced %>% filter(sample_date_decimal %in% timepoint) %>%  mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) %>% filter(Variants %in% markerPerVariant) %>% filter(Variants %notin% uniqMarkerPerVariant_not_detected) %>% group_by(ID, sample_date, sample_date_decimal, value.freq, value.depth, LocationID, LocationName, NUC) %>% summarize(Variants = paste(Variants, sep = ";", collapse = ";"), .groups = "drop_last") %>% ungroup() %>% dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC")-> ssdt_reduced
                        ssdt_reduced %>% filter(sample_date_decimal %in% timepoint) %>% filter(!grepl(";",Variants)) %>% group_by(Variants, ID, sample_date_decimal) %>% left_join(y = moi_uniq_marker_count, by = "Variants") %>% ungroup() %>% filter(value.freq>zeroo) %>% filter(value.depth > min.depth) %>% group_by(Variants, ID, sample_date_decimal) %>% mutate(N=n()) %>% mutate(r=N/n) %>% filter(r > minUniqMarkerRatio) %>% filter(N >= minUniqMarker) %>% summarize(.groups = "drop_last") -> ssdt2.2
                        uniqMarkerPerVariant.iterative[[c]] <- unique(ssdt2.2$Variants)
                        uniqMarkerPerVariant_could_be_detected.iterative[[c]] <- uniqMarkerPerVariant_could_be_detected$Variants
                        if(all(uniqMarkerPerVariant_could_be_detected.iterative[[c]] %in% uniqMarkerPerVariant_could_be_detected.iterative[[c-1]])){
                            c = 100
                        }
                    }
                    uniqMarkerPerVariant <- uniqMarkerPerVariant.iterative[[length(uniqMarkerPerVariant.iterative)]]

                    ## accept all markerPerVariant if ancestors to accepted by uniqMarkerPerVariant
                    uniqMarkerPerVariant_dealiased <- unlist(lapply(as.list(uniqMarkerPerVariant), dealias))
                    markerPerVariant_dealiased <- unlist(lapply(as.list(markerPerVariant), dealias))

                    specifiedLineages <- names(table(c(uniqMarkerPerVariant, markerPerVariant))[table(c(uniqMarkerPerVariant, markerPerVariant)) > 1])

                    if(length(unique(c(markerPerVariant, uniqMarkerPerVariant))) > 1){
                      for (ii in 1:length(markerPerVariant)){
                        if(markerPerVariant[ii] %in% specifiedLineages){
                          next
                        } else if (any(grepl(markerPerVariant_dealiased[ii], uniqMarkerPerVariant_dealiased, fixed = TRUE))){
                          specifiedLineages <- unique(c(specifiedLineages, markerPerVariant[ii]))
                        }
                      }
                      specifiedLineages <- unique(specifiedLineages)
                      specifiedLineages[!is.na(specifiedLineages)] -> specifiedLineages
                      rm(ssdt2, ssdt3, uniqMarkerPerVariant, markerPerVariant, ssdt_reduced, uniqMarkerPerVariant.iterative, uniqMarkerPerVariant_could_be_detected.iterative)

                    }
                } else if(opt$detectmode == "cons"){
                    ## count variant with most marker mutations
                    ## weight markers by number variants per marker
                    ssdt %>% rowwise() %>% mutate(v = -1 * (1/length(unlist(strsplit(Variants, split = ";")))) * log(1/length(unlist(strsplit(Variants, split = ";"))))) %>%  mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) %>% group_by(Variants) %>% summarize(n = n(), score = sum(v)) -> ssdt3
                    ssdt3 %>% arrange(desc(score)) %>% top_n(10) -> ssdt4
                    pnorm(abs(diff(ssdt4$score))[1], mean = mean(abs(diff(ssdt4$score))), sd = sd(abs(diff(ssdt4$score))), lower.tail = FALSE, log.p = FALSE) -> siggap
                    if( siggap > 0.05 ){
                      print(paste("WARNING: score difference between top rank (", abs(diff(ssdt4$score))[1],") and top 10 differences (", mean(abs(diff(ssdt4$score))),") not significant (", signif(siggap, digits = 3), ")"))
                    } else{
                      print(paste("WARNING: score difference between top rank (", abs(diff(ssdt4$score))[1],") and top 10 differences (", mean(abs(diff(ssdt4$score))),") looks fine (", signif(siggap, digits = 3), ")"))
                    }
                    specifiedLineages <- ssdt4$Variants[1]
                    rm(ssdt3, ssdt4, siggap)
                }

                print(paste("LOG: (", t, ")", timepoint_classic, roiname, paste("(", length(specifiedLineages), "detected lineages)")))
                print(paste0("LOG: detected lineages (", length(specifiedLineages), "): ", paste(specifiedLineages, sep = ", ", collapse = ", ")))

                if( length(specifiedLineages) > 0){
                  print(paste("LOG:", "TRY REGRESSION WITH", length(specifiedLineages), "LINEAGES"))
                  # join model matrix columns to measured variables
                  smmat <- as.data.frame(mmat)[,which(colnames(mmat) %in% c("NUC", specifiedLineages))]
                  if (length(specifiedLineages) > 1){
                    smmat$groupflag <- apply( as.data.frame(mmat)[,which(colnames(mmat) %in% c(specifiedLineages))], 1 , paste , collapse = "" )
                  } else{
                    smmat$groupflag <-  as.data.frame(mmat)[,which(colnames(mmat) %in% c(specifiedLineages))]
                  }
                  left_join(x = ssdt, y = smmat, by = c("NUC")) -> ssdt
                  ssdt %>% ungroup() %>% mutate(groupflag = paste(groupflag, ifelse(grepl(";", Variants), "M", "U"), sample_date_decimal)) -> ssdt

                  ## remove mutations which are not marker of any of the specifiedLineages
                  ## remove mutations which are marker of all of the specifiedLineages
                  if( sum(colnames(ssdt) %in% specifiedLineages) > 1){
                    ssdt[( rowSums(as.data.frame(ssdt)[,colnames(ssdt) %in% specifiedLineages ]) > 0 & rowSums(as.data.frame(ssdt)[,colnames(ssdt) %in% specifiedLineages ]) < length(specifiedLineages)),] -> ssdt
                  } else{
                    ssdt[as.data.frame(ssdt)[,colnames(ssdt) %in% specifiedLineages ] > 0,] -> ssdt
                  }

                  ## remove zeros and low depth
                  ssdt %>% filter(value.freq > 0)  %>% filter(value.depth > min.depth) -> ssdt

                  ## transform to avoid 1
                  ssdt %>% mutate(value.freq = (value.freq * (value.depth-1) + 0.5)/value.depth) -> ssdt

                  # remove "OTHERS"
                  ssdt %>% mutate(Variants = gsub("other;*", "", Variants)) -> ssdt

                  ## remove outlier per group flag
                  dim(ssdt)[1] -> mutationsBeforeOutlierRemoval
                  ssdt %>% group_by(groupflag) %>% mutate(iqr = IQR(value.freq)) %>% mutate(upperbond = quantile(value.freq, 0.75) + 1.5 * iqr, lowerbond = quantile(value.freq, 0.25) - 1.5 * iqr) %>% filter(value.freq <= upperbond & value.freq >= lowerbond) %>% ungroup() %>% dplyr::select(-"groupflag", -"iqr", -"upperbond", -"lowerbond") -> ssdt
                  dim(ssdt)[1] -> mutationsAfterOutlierRemoval
                  if(mutationsAfterOutlierRemoval < mutationsBeforeOutlierRemoval){
                    print(paste("LOG: ", mutationsBeforeOutlierRemoval-mutationsAfterOutlierRemoval, "mutations ignored since classified as outlier", "(", mutationsAfterOutlierRemoval, " remaining from", mutationsBeforeOutlierRemoval, ")"))
                  }

                  ## generate regression formula
                  if ( length(specifiedLineages) > 1 ){
                    formula <- as.formula(paste("value.freq", paste(specifiedLineages, collapse='+'), sep = " ~" ))
                  } else{
                    ssdt %>% filter(grepl(specifiedLineages, Variants)) -> ssdt
                    formula <- as.formula(paste("value.freq", "1", sep = " ~" ))
                  }

                  ## check if the current timepoint is still represented in the data
                  ssdt %>% filter(sample_date_decimal %in% timepoint) -> ssdt2
                  if(dim(ssdt2)[1] == 0){
                    print(paste("LOG: since no mutations found in timepoint of interest (", timepoint, ") regression will be skipped"))
                    next;
                  }

                  ## add weight (1/(dayDiff+1)*log10(sequencingdepth)
                  ssdt %>% rowwise() %>%  mutate(timeweight = 1/(abs(timePoints[t] - sample_date_decimal)*(leapYear(floor(timePoints[t]))) + 1)) %>% mutate(weight = log10(value.depth)*timeweight) -> ssdt

                  ## make regression
                  method = "SIMPLEX"
                  fit1 <- tryCatch(gamlss(formula, data = ssdt, family = SIMPLEX, trace = FALSE, weights = weight),error=function(e) e, warning=function(w) w)
                  if(any(grepl("warning|error", class(fit1)))){
                    method = "BE"
                    fit1 <- tryCatch(gamlss(formula, data = ssdt, family = BE, trace = FALSE, weights = weight),error=function(e) e, warning=function(w) w)
                    print(paste("LOG: fall back to BE"))
                    if(any(grepl("warning|error", class(fit1)))){
                      print(paste("LOG: BE did not converge either at", timePoints[t], " .. skipped"))
                      next; ## remove if unconverged BE results should be used
                    }
                  }
                  ssdt$fit1 <- predict(fit1, type="response", what="mu")

                  ssdt %>% dplyr::select(value.freq, all_of(specifiedLineages), fit1) -> ssdt_toOp
                  startValues <- starterV(ssdt_toOp)

                  O = 5
                  matrix(rep(0, length(specifiedLineages)*O), O) -> optimizedN
                  for (o in 1:O){
                    optim(par=startValues, fn=objfnct, data=ssdt_toOp, method = "L-BFGS-B", lower = 0, upper = 1) -> optimized
                    optimizedN[o,] <- optimized$par
                  }
                  optimized$par <- apply(optimizedN, 2, median)

                  ## normalize fitted value if sum is > 1
                  ssdt %>% filter(!grepl(";", Variants)) %>% filter(Variants %in% specifiedLineages) %>% mutate(fit1 = signif(fit1, 5))%>% dplyr::select("Variants", "fit1") %>% arrange(Variants) %>%  distinct()  %>% ungroup()  %>% mutate(T = sum(fit1)) %>% mutate(fit2 = ifelse(T>1, fit1/(T+0.00001), fit1)) -> ssdtFit
                  data.table(Variants = specifiedLineages, fit1 = signif(optimized$par, 5))  %>% arrange(Variants) %>%  distinct()  %>% ungroup()  %>% mutate(T = sum(fit1)) %>% mutate(fit2 = ifelse(T>1, fit1/(T+0.00001), fit1)) -> ssdtFit

                  ## add sample ID to output
                  if (any(ssdt$sample_date_decimal == timePoints[t])){
                    ssdt %>% filter(sample_date_decimal == timePoints[t]) %>% dplyr::select(ID) %>% distinct() -> sample_ID
                  } else{
                    data.frame(ID = ".") -> sample_ID
                    print(paste("Warning:", "no sample ID for", roi, "@", timePoints[t]))
                  }
                  ssdtFit$ID = rep(unique(sample_ID$ID)[length(unique(sample_ID$ID))], n = length(ssdtFit$Variants))

                  if(unique(ssdtFit$T)>1){
                    print(paste("LOG: fitted value corrected for >1: T =", unique(ssdtFit$T), "; method used =", method))
                  }

                  ssdt %>% ungroup() %>% mutate(T = unique(ssdtFit$T)) %>% mutate(fit2 = ifelse(T > 1, fit1/T, fit1)) %>% filter(sample_date_decimal == timePoints[t]) -> ssdt2
                  plantFullData <- rbind(plantFullData, data.table(variant = ssdt2$Variants, LocationID =  rep(roi, length(ssdt2$fit1)), LocationName  =  rep(roiname, length(ssdt2$fit1)), sample_date = as.character(ssdt2$sample_date), value = ssdt2$fit1, norm.value = ssdt2$fit2, marker = ssdt2$NUC, singlevalue = ssdt2$value.freq ))

                  ## save norm.fitted values into global DF; set all untested variants to 0
                  sampleFittedData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_id = character(), sample_date = character(), value = numeric(), norm.value = numeric() )
                  for (j in specifiedLineages){
                    ssdtFit %>% filter(Variants == j) -> extractedFt
                    if (dim(extractedFt)[1] > 0){
                      plantFittedData <- rbind(plantFittedData, data.table(variant = rep(j, length(extractedFt$fit2)), LocationID =  rep(roi, length(extractedFt$fit2)), LocationName  =  rep(roiname, length(extractedFt$fit2)), sample_id = extractedFt$ID, sample_date =  rep(as.character(unique(ssdt2$sample_date)), length(extractedFt$fit2)), value = extractedFt$fit1, norm.value = extractedFt$fit2 ))
                      sampleFittedData <- rbind(sampleFittedData, data.table(variant = rep(j, length(extractedFt$fit2)), LocationID =  rep(roi, length(extractedFt$fit2)), LocationName  =  rep(roiname, length(extractedFt$fit2)), sample_id = extractedFt$ID, sample_date =  rep(as.character(unique(ssdt2$sample_date)), length(extractedFt$fit2)), value = extractedFt$fit1, norm.value = extractedFt$fit2 ))
                    } else{
                      plantFittedData <- rbind(plantFittedData, data.table(variant = j, LocationID = roi, LocationName  = roiname, sample_id = extractedFt$ID, sample_date = as.character(unique(ssdt2$sample_date)), value = 0, norm.value = 0 ))
                      sampleFittedData <- rbind(sampleFittedData, data.table(variant = j, LocationID = roi, LocationName  = roiname, sample_id = extractedFt$ID, sample_date = as.character(unique(ssdt2$sample_date)), value = 0, norm.value = 0 ))
                    }
                  }
                  rm(formula, fit1, ssdt, ssdt2)
                } else {
                  print(paste("LOG:", "NOTHIN FOR REGRESSION AT", timePoints[t], "since number of lineages detected ==", length(specifiedLineages), "; ... ignore sample"))
                }


                if(exists("sampleFittedData")){
                    if(dim(sampleFittedData)[1] > 0){
                        left_join(x = variants_total, y = sampleFittedData %>% dplyr::select(variant, norm.value), by = c("Variants" = "variant")) -> sampleVarData
                        sampleVarData %>% mutate(norm.value = ifelse(is.na(norm.value), 0, norm.value)) -> sampleVarData
                        sampleVarData$sampleID <- sampleID
                        sampleVarData$LocationID <- roi
                        sampleVarData$LocationName <- roiname
                        sampleVarData$sample_date <- as.character(timepoint_classic)
                        sampleVarData %>% rename(variant = Variants, deduced.freq = norm.value) %>% dplyr::select(sampleID, LocationID, LocationName, sample_date, variant, deduced.freq) -> sampleVarData
                        globalVarData <- rbind(globalVarData, sampleVarData)


                        print(paste("LOG: plotting sankey", sampleID, roiname, timePoints_classic[t]))


                        sampleFittedData %>% mutate(latest = max(sample_date, na.rm = TRUE)) %>% filter(sample_date == latest) %>% filter(norm.value > 0) %>% summarize(variant = variant, freq = norm.value) -> sankey.dt

                        sampleFittedData %>% summarize(latest = max(sample_date, na.rm = TRUE)) -> sankey_date
                        sankey.dt %>% mutate(freq = freq/sum(freq)) -> sankey.dt
                        sankey.dt %>% group_by(variant = "B") %>% summarize(freq = 1-sum(freq)) -> sankey.dt0
                        rbind(sankey.dt, sankey.dt0) -> sankey.dt
                        sankey.dt %>% filter(freq > zeroo/5 | variant == "B") -> sankey.dt

                        sankey.dt$variant_dealias <- unlist(lapply(as.list(sankey.dt$variant), dealias))
                        for ( baseForColor in baseForColorsVariants ){
                          if(any(grepl(baseForColor, sankey.dt$variant_dealias, fixed = TRUE))){
                              if(baseForColor %notin% sankey.dt$variant_dealias){
                                  rbind(sankey.dt, data.frame(variant = realias(baseForColor), freq = 0, variant_dealias = baseForColor)) -> sankey.dt
                              }
                          }
                        }

                        sankey.dt %>% rowwise() %>% mutate(levels = length(unlist(strsplit(variant_dealias, split = "\\.")))) %>% group_by(1) %>% summarize( maxLev = max(levels)) -> maxLev
                        sankey.dt %>% rowwise() %>% mutate(long_variant = expandVar(variant_dealias, maxLev$maxLev)) -> sankey.dt
                        ancestor <- getTree(sankey.dt)
                        Nlevels <- max(unlist(lapply(ancestor, length)))
                        sankey.dt <- makeObs(sankey.dt, Nlevels, ancestor, sankeyPrecision)

                        sankey.dt %>% make_long(names(sankey.dt)) -> sankey.dtt

                        sankey.dtt %>% filter(!is.na(node)) -> sankey.dtt
                        max(as.numeric(gsub("level", "", sankey.dtt$next_x)), na.rm = TRUE) -> maxLev
                        sankey.dtt %>% group_by(x, node) %>% mutate(n = paste0("[", 100*n()/sankeyPrecision, "%]")) %>% rowwise() %>% mutate(level = as.numeric(gsub("level", "", x))) %>% mutate(label = ifelse(level == maxLev, paste(node, n), node)) -> sankey.dtt

                        sankey.dtt %>% rowwise() %>% mutate(node = ifelse(is.na(node), node, dealias(node))) %>% mutate(next_node = ifelse(is.na(next_node), next_node, dealias(next_node))) -> sankey.dtt

                        ggplot(sankey.dtt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = label)) + geom_sankey(flow.alpha = .6, node.color = "gray30", type ='alluvial') + geom_sankey_label(size = 3, color = "white", fill = "gray40", position = position_nudge(x = 0.05, y = 0), na.rm = TRUE, type ='alluvial', hjust = 0) + theme_sankey(base_size = 16) + labs(x = NULL) + ggtitle(roiname, subtitle = sankey_date) + theme(legend.position = "none", axis.text.x = element_blank(), plot.title = element_text(hjust = 0), plot.subtitle=element_text(hjust = 0)) + scale_fill_viridis_d(alpha = 1, begin = 0.025, end = .975, direction = 1, option = "D") -> pp

                        filename <- paste0(outdir, "/figs/sankey/",  paste('/sankeyPlot', sampleID, timePoints_classic[t], roi, sep="_"), ".pdf")
                        ggsave(filename = filename, plot = pp, width = 8, height = 4.5)
                        fwrite(as.list(c("sankey", "plot", sampleID, timePoints_classic[t], roiname, filename)), file = summaryDataFile, append = TRUE, sep = "\t")
                        rm(pp, filename, sankey.dt, sankey.dtt)
                    }


                    print(paste("LOG: infer expected mutation freqs", sampleID, roiname, timePoints_classic[t]))


                    ## infer expected AF based on AF per variant (as seen in pango) and deduced variant frequency
                    full_join(x = sampleFittedData, y = mstatgroup, by = c("variant" = "groupName"), suffix = c(".var", ".mut")) %>% dplyr::select(variant, value.var, sensitivity, AA_change, nucc) %>% mutate(sensitivity = ifelse(sensitivity > 0.8, sensitivity, 0)) %>% rowwise() %>% mutate(varAFcontribution = value.var*sensitivity) %>% group_by(nucc, AA_change) %>% summarize(AF = sum(varAFcontribution, na.rm = TRUE), .groups = "keep") %>% mutate(AF = ifelse(AF > 1, 1, AF)) -> expectedAF

                    ## extract observed AF
                    sewage_samps.sdt %>% filter(ID == sampleID) -> observedAF
                    full_join(x= expectedAF, y = observedAF, by = c("nucc" = "NUC")) %>% dplyr::select(nucc, AA_change, AF, value.freq, value.depth, Variants) %>% rename(expected.freq = AF, observed.freq = value.freq, observed.depth = value.depth, marker4Variants = Variants) -> expected_vs_observed_AF

                    ## contrast observed and expected AF
                    expected_vs_observed_AF %>% mutate(expected.freq = ifelse(is.na(expected.freq), 0, expected.freq)) %>% mutate(observed.freq = ifelse(is.na(observed.freq), 0, observed.freq)) -> expected_vs_observed_AF

                    ## filter sig. excess mutations based on beta distriubtion

                    beta_shape1 <- 2.2
                    pval_th <- 0.01
                    expected_vs_observed_AF  %>% filter(expected.freq > 0 | observed.freq > 0) %>% mutate(expected.freq.trans = ifelse(expected.freq < zeroo, zeroo, expected.freq)) %>% mutate(expected.freq.trans = (expected.freq.trans * (10000 - 1 + 0.5)/10000) ) %>% mutate(observed.freq.trans = ifelse(observed.freq < zeroo, zeroo, observed.freq)) %>% mutate(observed.freq.trans = (observed.freq.trans * (10000 - 1 + 0.5)/10000) ) -> expected_vs_observed_AF_transformed
                    expected_vs_observed_AF_transformed %>% mutate(pval = pbeta(observed.freq.trans, beta_shape1, betaParamFromMean(expected.freq.trans, beta_shape1), ncp = 0, lower.tail = FALSE, log.p = FALSE)) -> expected_vs_observed_AF_transformed
                    expected_vs_observed_AF_transformed$qval <- p.adjust(expected_vs_observed_AF_transformed$pval, method = "fdr")
                    expected_vs_observed_AF_transformed %>% filter(pval < pval_th) -> expected_vs_observed_AF_transformed_filtered
                    ggplot(expected_vs_observed_AF_transformed, aes(x = expected.freq, y = observed.freq)) + geom_smooth(method = "lm") + geom_abline(intercept = 0, slope = 1, color = "black") + theme_bw() + geom_point(size = 3, alpha = .5, aes(color = -1*log(qval))) + scale_color_viridis(option = "D", direction = -1, name = "pLog10(qval)") -> expected_vs_observed_AF_scatterPlot
                    filename <- paste0(outdir, "/figs/scatter/",  paste('/scatterPlot', sampleID, timePoints_classic[t], roi, sep="_"), ".pdf")
                    ggsave(filename = filename, plot = expected_vs_observed_AF_scatterPlot, width = 8, height = 4.5)
                    fwrite(as.list(c("scatter", "plot", sampleID, timePoints_classic[t], roiname, filename)), file = summaryDataFile, append = TRUE, sep = "\t")


                    globalAFdata   <- rbind(globalAFdata,
                                            data.table(
                                              sampleID      = sampleID,
                                              LocationID    = roi,
                                              LocationName  = roiname,
                                              sample_date   = as.character(timepoint_classic),
                                              nuc_mutation  = expected_vs_observed_AF_transformed$nucc,
                                              aa_mutation   = expected_vs_observed_AF_transformed$AA_change,
                                              observed      = expected_vs_observed_AF_transformed$observed.freq,
                                              expected      = expected_vs_observed_AF_transformed$expected.freq,
                                              excess        = expected_vs_observed_AF_transformed$observed.freq - expected_vs_observed_AF_transformed$expected.freq,
                                              pvalue        = expected_vs_observed_AF_transformed$qval
                                            )
                                      );

                    ## generate outbreak.info like heatmap
                    ###  col: enriched mutations
                    ###  row: variants for which at least on of the mutations is > 90%
                    expected_vs_observed_AF_transformed_filtered %>% group_by(nucc) %>% summarize(AA_change = AA_change, label = paste(nucc, AA_change), expected.freq = expected.freq, observed.freq = observed.freq, excess.freq = observed.freq-expected.freq, qval = qval) %>% arrange(desc(excess.freq)) -> outbreak.selection
                    data.table::melt(outbreak.selection, id.vars = c("label"), measure.vars = c("expected.freq", "observed.freq")) -> outbreak.freqs
                    if(any(dim(outbreak.selection) == 0)){
                      next;
                    }
                    mstat %>% filter(nucc %in% outbreak.selection$nucc) -> outbreak.dt
                    data.table::dcast(outbreak.dt, formula = ID ~ paste(nucc, AA_change), value.var = "sensitivity", fill = 0) -> outbreak.dt
                    data.table::melt(outbreak.dt, variable.name = "label", value.name = "AF.freq") -> outbreak.dt
                    outbreak.dt %>% group_by(ID) %>% mutate(max = max(AF.freq, na.rm = TRUE)) %>% filter(max >= min(max(outbreak.dt$AF.freq), .9)) %>% dplyr::select(-"max") -> outbreak.dt
                    outbreak.dt %>% mutate(dealiasID = dealias(ID)) -> outbreak.dt
                    # remove variants which are represented by a ancestor with same fingerprint
                    outbreak.dt %>% group_by(ID) %>%arrange(label) %>% mutate(fingerprint = paste(ifelse(AF.freq > 0.5, 1, 0), collapse = "")) %>% group_by(fingerprint) %>% mutate(groupIds = paste(unique(dealiasID), sep = ";", collapse = ";")) %>% rowwise() %>% mutate(keep = collapse2mother(dealiasID, groupIds)) %>% filter(keep) -> outbreak.dt

                    # remove mutations which are not seen in all three plots
                    table(c(unique(as.vector(outbreak.dt$label)), unique(outbreak.freqs$label), unique(outbreak.selection$label))) -> setCompletness
                    names(setCompletness[setCompletness > 2]) -> labelsToUse

                    outbreak.dt %>% filter(label %in%  labelsToUse) -> outbreak.dt
                    outbreak.freqs %>% filter(label %in%  labelsToUse) -> outbreak.freqs
                    outbreak.selection %>% filter(label %in%  labelsToUse) -> outbreak.selection




                    outbreak.dt$label <- factor(outbreak.dt$label, levels = outbreak.selection$label)
                    outbreak.freqs$label <- factor(outbreak.freqs$label, levels = outbreak.selection$label)
                    #outbreak.freqs %>% filter(!is.na(label) ) -> outbreak.freqs
                    outbreak.selection$label <- factor(outbreak.selection$label, levels = outbreak.selection$label)
                    #outbreak.selection %>% filter(!is.na(label)) -> outbreak.selection


                    ggplot(data = outbreak.dt, aes(x = label, y = ID, fill = AF.freq)) + geom_tile(color = "white", alpha = 0.8) + scale_fill_viridis(name = "Allele\nfrequency", trans = "sqrt", option = "B", begin = 0.1, end = 0.9) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("Excess Mutations") + ylab("Variants") -> outbreakInfoPlot1

                    ggplot(data = outbreak.selection, aes(x = label, y = excess.freq, fill = -1*log10(qval))) + geom_col(width = 0.66) + theme_bw() + xlab("") + ylab("Excess AF") + scale_fill_fermenter(name = "pLog10(p-value)", palette = "Oranges", direction = 1, na.value = "grey33") -> outbreakInfoPlot2

                    ggplot(data = outbreak.freqs, aes(x = label, y = value, fill = variable)) + geom_col(position = "dodge", width = .6) + theme_bw() + xlab("") + ylab("AF") + scale_fill_manual(label = c("Expected", "Observed"), breaks = c("expected.freq", "observed.freq"), values = c("#81b29a", "#e07a5f"), name = "Allele\nfrequency") -> outbreakInfoPlot3

                    outbreakInfoPlot3 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + outbreakInfoPlot2 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + outbreakInfoPlot1 + plot_layout(nrow = 3, byrow = FALSE, heights = c(1, 3,10)) -> outbreakInfoPlot
                    filename <- paste0(outdir, "/figs/excessmutation/",  paste('/excessmutationPlot', sampleID, timePoints_classic[t], roi, sep="_"), ".pdf")
                    ggsave(filename = filename, plot = outbreakInfoPlot, width = 8, height = 16)
                    fwrite(as.list(c("excessmutation", "plot", sampleID, timePoints_classic[t], roiname, filename)), file = summaryDataFile, append = TRUE, sep = "\t")

                    filename <- paste0(outdir, "/figs/excessmutation/",  paste('/excessmutationData', sampleID, timePoints_classic[t], roi, sep="_"), ".Rdata")
                    save(list = c("outbreak.dt", "outbreak.selection", "outbreak.freqs"), file = filename)
                    fwrite(as.list(c("excessmutation", "data", sampleID, timePoints_classic[t], roiname, filename)), file = summaryDataFile, append = TRUE, sep = "\t")

                    #paste0("[", ceiling(length(unique(globalAFdata$nuc_mutation))/2), "-of:", paste(unique(globalAFdata$nuc_mutation), collapse = ", "), "]")

                } else{
                    print(paste("LOG: no variant deduced", sampleID, roiname, timePoints_classic[t]))
                }

        }

        globalFittedData <- rbind(plantFittedData, globalFittedData)
        globalFullData <- rbind(plantFullData, globalFullData)

        rm(plantFittedData, plantFullData, colorAssignment)
    }
}

print(paste("PROGRESS: loop over WWTP finished"))
## dump global data into output file
print(paste("PROGRESS: writing result tables"))
globalFittedData %>% distinct() -> globalFittedData
globalFullData %>% distinct() -> globalFullData
globalAFdata %>% distinct() -> globalAFdata
globalVarData %>% distinct() -> globalVarData

fwrite(globalFittedData , file = paste0(outdir, "/globalFittedData.csv"), sep = "\t")
#fwrite(globalFullData, file = paste0(outdir, "/globalFullData.csv"), sep = "\t")
fwrite(globalAFdata, file = paste0(outdir, "/globalAFdata.csv"), sep = "\t")
fwrite(globalVarData, file = paste0(outdir, "/globalVarData.csv"), sep = "\t")


## generate over view plots
if(1){
  ## filter mutations which are significant (p<0.1) and seen more >1 time
  globalAFdata %>% mutate(excess = as.numeric(excess), N = length(unique(sampleID))) %>% filter(pvalue < 0.1) %>% group_by(nuc_mutation) %>% mutate(n = length(unique(sampleID))) %>% mutate(r = n/N) -> globalAFdata.plotting
  globalAFdata.plotting %>% filter(n > 1) -> globalAFdata.plotting

   ## cluster mutation with binary = Jacchard distance
  globalAFdata.plotting %>% mutate(label = paste0(nuc_mutation, " [", aa_mutation, "]")) -> globalAFdata.plotting
  reshape2::dcast(globalAFdata.plotting, label~sampleID, value.var = "excess", fill = 0) -> mat
  rownames(mat) <- mat$label
  mat[-1] -> mat
  hclust(dist(mat, method = "binary")) -> distancing
  cutree(distancing, h = .5) -> clustering
  left_join(x = globalAFdata.plotting, y = data.table(label = names(clustering), cluster = clustering), by = "label") -> globalAFdata.plotting  

  ## sort dt according clustering and date
  globalAFdata.plotting$label <- factor(globalAFdata.plotting$label, levels = distancing$labels[distancing$order])
  globalAFdata.plotting %>% rowwise() %>% mutate(sample_label = paste(sampleID, sample_date, sep = " @ ")) -> globalAFdata.plotting   
  globalAFdata.plotting$sample_label <- factor(globalAFdata.plotting$sample_label, levels=unique(globalAFdata.plotting$sample_label[order(globalAFdata.plotting$sample_date)]))

  ggplot(data = globalAFdata.plotting, aes(x = label, y = sample_label)) + geom_tile(aes(fill = excess)) + facet_grid(LocationName~cluster, scale="free", space = "free") + theme_bw() + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_fermenter(palette = 7, direction = 1, name = "excess\nallele freq.") + xlab("") + ylab("") -> pp
  plot.width  <- 4+length(unique(globalAFdata.plotting$label))/5+length(unique(globalAFdata.plotting$cluster))/20
  plot.height <- 4+length(unique(globalAFdata.plotting$sample_label))/5+length(unique(globalAFdata.plotting$LocationID))/10


  filename <- paste0(outdir, "/figs/overview/",  paste('/heatmap', "excessMutations", sep="_"), ".pdf")
  ggsave(filename = filename, plot = pp, width = plot.width, height = plot.height)
  fwrite(as.list(c("overview", "heatmap", "all", "any", "everywhere", filename)), file = summaryDataFile, append = TRUE, sep = "\t")


  ## add detected variants per sample    
  globalAFdata.plotting -> globalAFdata.plotting.2
  globalVarData %>% filter(deduced.freq > 0) -> globalVarData.plotting   
  left_join(x = globalAFdata.plotting.2, y = globalVarData.plotting, by = c("sampleID", "LocationID", "LocationName", "sample_date"), multiple = "all") -> globalAFdata.plotting.2
  globalAFdata.plotting.2 %>% group_by(sampleID) %>% summarize(n = n()) -> samples

  # test for each cluster - variant combination if they significantly co-occure (fisher exact test)
  data.table(clusterIdx = character(), clusterMember = character(), variant = character(), coocurence = numeric(), oddsRatio = numeric(), pvalue = numeric()) -> collectEnrichment

  for (c in unique(globalAFdata.plotting.2$cluster)){
    globalAFdata.plotting.2 %>% filter(cluster == c) -> globalAFdata.plotting.2_
    clustmembers <- length(unique(globalAFdata.plotting.2_$label))

    globalAFdata.plotting.2_ %>% group_by(sampleID) %>% summarize(n = length(unique(label))) %>% filter(n > clustmembers/2) -> clustdetected

    for (v in unique(globalAFdata.plotting.2_$variant)){
      globalAFdata.plotting.2 %>% filter(variant == v) -> globalAFdata.plotting.2__
      globalAFdata.plotting.2__ %>% group_by(sampleID) %>% filter(variant == v) %>% dplyr::select(sampleID) -> vardetected

      clust_and_var <- sum(clustdetected$sampleID %in% vardetected$sampleID)
      samples_no_var <- sum(samples$sampleID %notin% vardetected$sampleID)
      samples_no_clust <- sum(samples$sampleID %notin% clustdetected$sampleID)
      no_clust_no_var <- sum(samples$sampleID[samples$sampleID %notin% vardetected$sampleID] %in% samples$sampleID[samples$sampleID %notin% clustdetected$sampleID])
      confusionmat <- matrix(c(no_clust_no_var,samples_no_clust,samples_no_var,clust_and_var), nrow = 2)
      fisher.test(confusionmat) ->fishStats
      rbind(collectEnrichment, data.table(clusterIdx = c, clusterMember = paste(unique(globalAFdata.plotting.2_$label), collapse = ';'), variant = v, coocurence = clust_and_var, oddsRatio = fishStats$estimate, pvalue = fishStats$p.value)) -> collectEnrichment

    }
  }
  p.adjust(collectEnrichment$pvalue, method = "fdr") -> collectEnrichment$qvalue
  collectEnrichment %>% mutate(shape = ifelse(oddsRatio>1,">1","<1")) -> collectEnrichment
  ggplot(data = collectEnrichment, aes(x = variant, y = clusterIdx)) + geom_point(aes(size = oddsRatio, fill = round(-1*log10(pvalue), digits = 3), shape = shape), color = "grey80") + coord_fixed() + theme_bw() + theme(legend.position="left", axis.text.x = element_text(angle = 90, hjust = 1), legend.direction = "vertical", legend.box = "vertical") + xlab("Variant") + ylab("Mutation Cluster") + scale_fill_distiller(name = "pLog10(p-value)", direction = 1) + scale_shape_manual(breaks = c("<1", ">1"), values = c(25, 24), name = "direction") -> pq

  plot.width   <- 2 + length(unique(collectEnrichment$variant))/2
  plot.height  <- 2 + length(unique(collectEnrichment$clusterIdx))/3

  filename <- paste0(outdir, "/figs/overview/",  paste('/bubble', "clusterVariantCooccurence", sep="_"), ".pdf")
  ggsave(filename = filename, plot = pq, width = plot.width, height = plot.height)
  fwrite(as.list(c("overview", "bubble", "all", "any", "everywhere", filename)), file = summaryDataFile, append = TRUE, sep = "\t")


  collectEnrichment  %>% mutate(clusterMembers = strsplit(as.character(clusterMember), ";")) %>% unnest() %>% separate(col = clusterMembers, into = c("nuc", "aa"), sep = " ") %>% group_by(clusterIdx) %>% summarize(members = covspectrumLink(nuc)) -> collectEnrichment_covspectrumLink
  left_join(x = collectEnrichment, y = collectEnrichment_covspectrumLink, by = "clusterIdx") -> collectEnrichment_covspectrumLink
  filename <- paste0(outdir, "/figs/overview/",  paste('/table', "clusterVariantCooccurence", sep="_"), ".csv")
  fwrite(collectEnrichment_covspectrumLink, file=filename, sep = "\t")

  ## add detected variants per sample    
  globalAFdata.plotting -> globalAFdata.plotting.3
  globalVarData %>% filter(deduced.freq > 0) -> globalVarData.plotting
  left_join(x = globalAFdata.plotting.3, y = globalVarData.plotting, by = c("sampleID", "LocationID", "LocationName", "sample_date"), multiple = "all") -> globalAFdata.plotting.3
  globalAFdata.plotting.3 %>% group_by(sampleID) %>% summarize(n = n()) -> samples

  # test for each mutation - variant combination if they significantly co-occure (fisher exact test)
  data.table(mutation = character(), variant = character(), coocurence = numeric(), oddsRatio = numeric(), pvalue = numeric()) -> collectEnrichment

  for (c in unique(globalAFdata.plotting.3$nuc_mutation)){
    globalAFdata.plotting.3 %>% filter(nuc_mutation == c) -> globalAFdata.plotting.3_

    globalAFdata.plotting.3_ %>% group_by(sampleID) %>% summarize(n = length(unique(label))) %>% filter(n >= 1) -> mutdetected

    for (v in unique(globalAFdata.plotting.3_$variant)){
      globalAFdata.plotting.3 %>% filter(variant == v) -> globalAFdata.plotting.3__
      globalAFdata.plotting.3__ %>% group_by(sampleID) %>% filter(variant == v) %>% dplyr::select(sampleID) -> vardetected

      mut_and_var <- sum(mutdetected$sampleID %in% vardetected$sampleID)
      samples_no_var <- sum(samples$sampleID %notin% vardetected$sampleID)
      samples_no_mut <- sum(samples$sampleID %notin% mutdetected$sampleID)
      no_mut_no_var <- sum(samples$sampleID[samples$sampleID %notin% vardetected$sampleID] %in% samples$sampleID[samples$sampleID %notin% mutdetected$sampleID])
      confusionmat <- matrix(c(no_mut_no_var,samples_no_mut,samples_no_var,mut_and_var), nrow = 2)
      fisher.test(confusionmat) ->fishStats
      rbind(collectEnrichment, data.table(mutation = paste(unique(globalAFdata.plotting.3_$label), collapse = ';'), variant = v, coocurence = mut_and_var, oddsRatio = fishStats$estimate, pvalue = fishStats$p.value)) -> collectEnrichment

    }
  }
  p.adjust(collectEnrichment$pvalue, method = "fdr") -> collectEnrichment$qvalue
  collectEnrichment %>% mutate(shape = ifelse(oddsRatio>1,">1","<1")) -> collectEnrichment
  ggplot(data = collectEnrichment, aes(x = variant, y = mutation)) + geom_point(aes(size = oddsRatio, fill = round(-1*log10(pvalue), digits = 3), shape = shape), color = "grey80") + coord_fixed() + theme_bw() + theme(legend.position="left", axis.text.x = element_text(angle = 90, hjust = 1), legend.direction = "vertical", legend.box = "vertical") + xlab("Variant") + ylab("Mutation mutation") + scale_fill_distiller(name = "pLog10(p-value)", direction = 1) + scale_shape_manual(breaks = c("<1", ">1"), values = c(25, 24), name = "direction") -> pq

  plot.width   <- 2 + length(unique(collectEnrichment$variant))/2
  plot.height  <- 2 + length(unique(collectEnrichment$mutation))/3

  filename <- paste0(outdir, "/figs/overview/",  paste('/bubble', "mutationVariantCooccurence", sep="_"), ".pdf")
  ggsave(filename = filename, plot = pq, width = plot.width, height = plot.height)
  fwrite(as.list(c("overview", "bubble", "all", "any", "everywhere", filename)), file = summaryDataFile, append = TRUE, sep = "\t")


}

timestamp()
