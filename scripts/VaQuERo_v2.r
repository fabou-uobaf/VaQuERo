## load libraries
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("gamlss"))
suppressPackageStartupMessages(library("ggmap"))
suppressPackageStartupMessages(library("tmaptools"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("betareg"))
suppressPackageStartupMessages(library("ggspatial"))
suppressPackageStartupMessages(library("sf"))
suppressPackageStartupMessages(library("rnaturalearth"))
suppressPackageStartupMessages(library("rnaturalearthdata"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("rjson"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("ggsankey"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("viridis"))


timestamp()

# get Options
option_list = list(
  make_option(c("--dir"), type="character", default="ExampleOutput",
              help="Directory to write results [default %default]", metavar="character"),
  make_option(c("--country"), type="character", default="Austria",
              help="Name of country used to produce map [default %default]", metavar="character"),
  make_option(c("--bbsouth"), type="double", default=46.38,
              help="Bounding box most south point [default %default]", metavar="character"),
  make_option(c("--bbnorth"), type="double", default=49.01,
              help="Bounding box most norther point [default %default]", metavar="character"),
  make_option(c("--bbwest"), type="double", default=9.53,
              help="Bounding box most western point [default %default]", metavar="character"),
  make_option(c("--bbeast"), type="double", default=17.15,
              help="Bounding box most easter point [default %default]", metavar="character"),
  make_option(c("--metadata"), type="character", default="data/metaDataSub.tsv",
              help="Path to meta data input file [default %default]", metavar="character"),
  make_option(c("--marker"), type="character", default="resources/mutations_list.csv",
              help="Path to marker mutation input file [default %default]", metavar="character"),
  make_option(c("--smarker"), type="character", default="resources/mutations_special.csv",
              help="Path to special mutation input file [default %default]", metavar="character"),
  make_option(c("--pmarker"), type="character", default="resources/mutations_problematic_all.csv",
              help="Path to problematic mutation input file, which will be ignored throughout the analysis [default %default]", metavar="character"),
  make_option(c("--data2"), type="character", default="data/mutationDataSub_sparse.tsv",
              help="Path to data input file in deprecated sparse matrix format [default %default]", metavar="character"),
  make_option(c("--data"), type="character", default="data/mutationDataSub.tsv.g",
              help="Path to data input file in tidy table format [default %default]", metavar="character"),
  make_option(c("--inputformat"), type="character", default="tidy",
              help="Input data format. Should be 'tidy' or 'sparse' [default %default]", metavar="character"),
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
  make_option(c("--recent"), type="integer", default=99,
              help="How old (in days) most recent sample might be to be still considered in overview maps [default %default]", metavar="character"),
  make_option(c("--plottp"), type="integer", default=3,
              help="Produce timecourse plots only if more than this timepoints are available [default %default]", metavar="character"),
  make_option(c("--minuniqmark"), type="integer", default=1,
              help="Minimal absolute number of uniq markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--minuniqmarkfrac"), type="double", default=0.4,
              help="Minimal fraction of uniq markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--minqmark"), type="integer", default=3,
              help="Minimal absolute number of markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--minmarkfrac"), type="double", default=0.4,
              help="Minimal fraction of markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--smoothingsamples"), type="integer", default=2,
              help="Number of previous timepoints use for smoothing [default %default]", metavar="character"),
  make_option(c("--smoothingtime"), type="double", default=2,
              help="Previous timepoints for smoothing are ignored if more days than this days apart [default %default]", metavar="character"),
  make_option(c("--voi"), type="character", default="BA.4,BA.5,BA.1,BA.2,B.1.617.2,B.1.1.7",
              help="List of variants which should be plotted in more detail. List separated by semicolon [default %default]", metavar="character"),
  make_option(c("--highlight"), type="character", default="BA.4,BA.5,BA.1,BA.2,B.1.617.2,B.1.1.7",
              help="List of variants which should be plotted at the bottom axis. List separated by semicolon [default %default]", metavar="character"),
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
  opt$metadata = "data/metaData_subset.csv"
  opt$data2="data/sewage_samples_merge_samp_lofreq_dp_AF001.snpeff_afs.dp.pq.tab.gz"
  opt$inputformat = "sparse"
  opt$marker="VaQuERo/resources/mutations_list_grouped_pango_2022-11-29_Austria.csv"
  opt$smarker="VaQuERo/resources/mutations_special_2022-09-13.csv"
  opt$pmarker="VaQuERo/resources/mutations_problematic_vss1_v3.csv"
  opt$zero=0.02
  opt$depth=250
  opt$minuniqmark=1
  opt$minuniqmarkfrac=0.5
  opt$minqmark=4
  opt$minmarkfrac=0.5
  opt$smoothingsamples=2
  opt$smoothingtime=2
  opt$voi="B.1.1.7,BA.1,BA.2,BA.4,BA.5,B.1.617.2,XBB,BF.7,BM.1.1,BA.2.12.1,BN.1,BA.2.75.2"
  opt$highlight="XBB,BF.7,BM.1.1,BA.2.12.1,BN.1,BA.2.75.2,BQ.1.1"
  opt$colorBase="B.1.617.2,BA.1,BA.2,BA.4,BA.5,BQ.1.1"
  opt$recent <- 9999
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




## print parameter to Log

print("##~LOG~PARAMETERS~####################")
print(opt)
print("##~LOG~PARAMETERS~####################")
writeLines("\n\n\n")

## read in config file to overwrite all para below

## set variable parameters
summaryDataFile <- paste0(opt$dir, "/summary.csv")
markermutationFile  <- opt$marker
specialmutationFile <- opt$smarker
problematicmutationFile <- opt$pmarker
sankeyPrecision <- 1000

mapCountry  <- opt$country
mapMargines <- c(opt$bbsouth,opt$bbwest,opt$bbnorth,opt$bbeast)
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
if( ! dir.exists(paste0(outdir, "/figs/specialMutations"))){
  dir.create(paste0(outdir, "/figs/specialMutations"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/fullview"))){
  dir.create(paste0(outdir, "/figs/fullview"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/detail"))){
  dir.create(paste0(outdir, "/figs/detail"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/stackview"))){
  dir.create(paste0(outdir, "/figs/stackview"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/maps"))){
  dir.create(paste0(outdir, "/figs/maps"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/sankey"))){
  dir.create(paste0(outdir, "/figs/sankey"), showWarnings = FALSE)
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
recentEnought = opt$recent # sample with less than this days in the past are not used for the summary map
tpLimitToPlot = opt$plottp # produce plot only if at least than n time points
minUniqMarker <- opt$minuniqmark # minmal absolute number of uniq markers that variant is considered detected
minUniqMarkerRatio <- opt$minuniqmarkfrac # minmal fraction of uniq markers that variant is considered detected (0.25)
minMarker <- opt$minqmark # minmal absolute number of sensitive markers that variant is considered detected
minMarkerRatio <- opt$minmarkfrac # minmal fraction of sensitive markers that variant is considered detected (0.25)
timeLag <- opt$smoothingsamples # number of previous timepoints use for smoothing
timeStart <- 1 # set to 1 if all time points should be considered; 2 if first should not be considered
timeLagDay <- opt$smoothingtime # previous timepoints for smoothing are ignored if more days before
VoI                 <- unlist(str_split(opt$voi, pattern=c(";",",")))
highlightedVariants <- unlist(str_split(opt$highlight, pattern=c(";",","))[[2]])
baseForColorsVariants <- unlist(str_split(opt$colorBase, pattern=c(";",","))[[2]])
baseForColorsVariants <- unlist(lapply(as.list(baseForColorsVariants), dealias))
colorSets <- c("Blues", "Greens", "Purples", "Oranges", "Reds", "Greys")


## define global variables to fill
globalFittedData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_id = character(), sample_date = character(), value = numeric() )
globalFullData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_date = character(), value = numeric(), marker = character(), singlevalue = numeric() )
globalFullSoiData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_date = character(), marker = character(), value = numeric())


## read mutations of interest from file
print(paste0("PROGRESS: read mutation files "))
moi <- fread(file = markermutationFile)
unite(moi, NUC, c(4,3,5), ALT, sep = "", remove = FALSE) -> moi

## read special mutations of interest from file
soi <- fread(file = specialmutationFile)
unite(soi, NUC, c(4,3,5), ALT, sep = "", remove = FALSE) -> soi
soi %>% dplyr::select("Variants","NUC", "AA") -> special

## read special mutations of interest from file
poi <- fread(file = problematicmutationFile)
unite(poi, NUC, c(4,3,5), ALT, sep = "", remove = FALSE) -> poi
poi %>% dplyr::select("PrimerSet","NUC", "AA") -> problematic

## remove all problematic sites from soi and moi
moi %>% filter(NUC %notin% poi$NUC) -> moi
soi %>% filter(NUC %notin% poi$NUC) -> soi

## define subset of marker and special mutations
allmut   <- unique(c(moi$NUC, soi$NUC))
exsoimut <- allmut[allmut %notin% moi$NUC]
soimut   <- soi$NUC

## read in meta data
print(paste0("PROGRESS: read and process meta data "))
metaDT       <- fread(file = opt$metadata)
unique(metaDT$BSF_sample_name) -> sampleoi

## check if all voi can be detected in principle given used parameters
moi %>% dplyr::select("Variants") %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest() %>% distinct() %>% filter(Variants != "other")-> variants_total
moi %>% filter(!grepl(";", Variants)) %>% group_by(Variants) %>% summarise(n = n()) %>% filter(n >= minUniqMarker) -> variants_w_enough_uniq

variants_total$Variants[variants_total$Variants %notin% variants_w_enough_uniq$Variants] -> variants_futile
if(length(variants_futile) > 0){
  for (c in seq_along(variants_futile)){
    message(paste("Warning(s):", variants_futile[c], "cannot be detected directly since less than", minUniqMarker, "unique markers defined."))
  }
}

# check for same date, same location issue
# introduce artifical decimal date
metaDT %>% group_by(LocationID, sample_date) %>% mutate(n = rank(BSF_sample_name)-1)  %>% ungroup() %>% rowwise() %>% mutate(sample_date_decimal = decimalDate(sample_date, n)) %>% dplyr::select(-n) -> metaDT

# define which locations are used for the report
metaDT %>% dplyr::select("LocationID") %>% distinct() -> locationsReportedOn


# Generate map of all locations sequenced in current run

## determine last sequencing batch ID
as.character(metaDT %>% filter(!is.na(BSF_run)) %>% filter(!is.na(BSF_start_date)) %>% group_by(BSF_run, BSF_start_date) %>% summarize(.groups = "drop") %>% ungroup() %>% filter(as.Date(BSF_start_date) == max(as.Date(BSF_start_date))) %>% dplyr::select("BSF_run")) -> last_BSF_run_id
print(paste("LOG: last_BSF_run_id", last_BSF_run_id))

## modify status of samples in last run based on number N in consensus
metaDT %>% filter(BSF_run == last_BSF_run_id) %>% filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date) -> mapSeqDT
mapSeqDT %>% mutate(status = ifelse(is.na(status), "fail", status)) -> mapSeqDT
mapSeqDT %>% mutate(status = ifelse(N_in_Consensus > N_in_Consensus_filter & status == "pass", "fail", status)) -> mapSeqDT
mapSeqDT %>% mutate(status = ifelse(status == "fail" & N_in_Consensus < N_in_Consensus_detection_filter, "detected", status)) -> mapSeqDT

## remove samples with "in_run" flag set in status
if(any(grepl("in_run", mapSeqDT$status))){
  in_run_count <- sum(grepl("in_run", mapSeqDT$status))
  print(paste0("WARNING: status 'in_run' found ", in_run_count," times; set to 'fail'"))
  mapSeqDT %>% mutate(status = ifelse(status == "in_run", "fail", status)) -> mapSeqDT
}

# print log how many samples from last run passed filter
print(paste0("LOG: print STATUS counts: "))
print(table(mapSeqDT$status))

## generate empty map
print(paste0("PROGRESS: print STATUS map"))
timestamp()

World <- ne_countries(scale = "medium", returnclass = "sf")
Country <- subset(World, name_sort == mapCountry)

s <- ggplot()
s <- s + geom_sf(data = World, fill = "grey95")
s <- s + geom_sf(data = Country, fill = "antiquewhite")
s <- s + theme_minimal()
s <- s + coord_sf(ylim = mapMargines[c(1,3)], xlim = mapMargines[c(2,4)], expand = FALSE)
s <- s + annotation_scale(location = "bl")
s <- s + annotation_north_arrow(location = "tl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering)
s <- s + theme(axis.text = element_blank(), legend.direction = "vertical", legend.box = "horizontal", legend.position = "bottom")
s <- s + geom_point(data=mapSeqDT, aes(y=dcpLatitude, x=dcpLongitude, shape = status, fill = status, size = as.numeric(connected_people)), alpha = 0.25)
s <- s + scale_fill_manual(values = c("detected" = "deepskyblue", "fail" = "firebrick", "pass" = "limegreen"), name = "Seq. Status")
s <- s + guides(size = guide_legend(title = "Population", nrow = 2))
s <- s + scale_size(range = c(2, 6), labels = scales::comma, breaks = c(50000, 100000, 500000, 1000000))
s <- s + scale_shape_manual(values = c("detected" = 24, "fail" = 25, "pass" = 21), name = "Seq. Status") # set to c(24,25,21)
s <- s + xlab("") + ylab("")
filename <- paste0(outdir, "/figs/maps/STATUS.pdf")
ggsave(filename = filename, plot = s)
fwrite(as.list(c("statusmapplot", "status", filename)), file = summaryDataFile, append = TRUE, sep = "\t")
rm(s, filename)


## exchange VoI with merged VoI
VoIm <- unique(moi$Variants[grep(";", moi$Variants, invert = TRUE)])
for (i in 1:length(VoI)){
    if (VoI[i] %in% VoIm){
      VoI[i] <- VoIm[VoI[i] == VoIm]
    }
}
for (i in 1:length(highlightedVariants)){
  if (highlightedVariants[i] %in% VoIm){
    highlightedVariants[i] <- VoIm[highlightedVariants[i] == VoIm]
  }
}


## generate general model matrix
moi %>% filter(NUC %notin% exsoimut) %>% dplyr::select("Variants","NUC") -> marker
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
  sewage_samps %>% filter(POS %in% unique(sort(c(moi$Postion, soi$Postion)))) -> sewage_samps

  ## get sample names
  sample_names = grep("\\.AF$", names(sewage_samps), value = TRUE)
  sample_names = gsub("\\.AF$","",sample_names)

  ## remove all positions which are not mutations of interest
  unite(sewage_samps, NUC, c(3,2,4), ALT, sep = "", remove = FALSE) -> sewage_samps
  sewage_samps[sewage_samps$NUC %in% allmut,] -> sewage_samps

  ## remove all positions which are problematic
  sewage_samps[sewage_samps$NUC %notin% poi$NUC,] -> sewage_samps

  ## set missing mutation frequencies with 0
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
  sewage_samps.dt$value.depth[is.na(sewage_samps.dt$value.depth)] <- 333  ###########!!! CHANGE TO 0
  sewage_samps.dt$value.qual[is.na(sewage_samps.dt$value.qual)] <- 0

  ## remove all positions which are not mutations of interest
  sewage_samps.dt %>% filter(NUC %in% allmut) -> sewage_samps.dt

  ## remove all samples which are not specified in metadata
  sewage_samps.dt %>% filter(ID %in% sampleoi) -> sewage_samps.dt

  rm(sewage_samps)
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

### FROM HERE LOOP OVER EACH SEWAGE PLANTS
print(paste("PROGRESS: start to loop over WWTP"))
timestamp()
for (r in 1:length(unique(sewage_samps.dt$LocationID))) {
    roi = unique(sewage_samps.dt$LocationID)[r]
    roiname = unique(sewage_samps.dt$LocationName)[r]
    print(paste("PROGRESS:", roiname, roi))

    ### define data collector
    plantFittedData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_id = character(), sample_date = character(), value = numeric() )
    plantFullData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_date = character(), value = numeric(), marker = character(), singlevalue = numeric() )
    plantFullSoiData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_date = character(), marker = character(), value = numeric())

    ### filter data set for current location
    sewage_samps.dt %>% filter(LocationID == roi) %>% dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC") %>% filter(NUC %notin% exsoimut) -> sdt
    sewage_samps.dt %>% filter(LocationID == roi) %>% dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC") %>% filter(NUC %in% soimut) -> spemut_sdt

    ### collect frequency of special mutations
    left_join(x = spemut_sdt, y = soi, by = "NUC") %>% rowwise() %>% mutate(marker = paste(sep = "|", NUC, paste(Gene, AA, sep = ":"), paste0("[", Variants.y, "]")) ) %>% dplyr::select(Variants.y, LocationID, LocationName, sample_date, marker, value.freq) -> plantFullSoiData
    colnames(plantFullSoiData) <- c("variant", "LocationID", "LocationName", "sample_date", "marker", "value")
    rbind(plantFullSoiData, globalFullSoiData) -> globalFullSoiData

    ## count how many sample from this plant were in the last BSF run
    metaDT %>% filter(BSF_sample_name %in% sdt$ID) %>% filter(BSF_run %in% last_BSF_run_id) %>% dplyr::select(BSF_run) %>% group_by(BSF_run) %>% summarize(n = n()) -> count_last_BSF_run_id

    ## define variants for which at least minUniqMarker are found in any timepoint
    ## at least minUniqMarkerRatio of all uniq markers are found in any timepoint
    ## define marker for which at least minMarker are found in any timepoint
    ## at least minMarkerRatio of all uniq markers are found in any timepoint
    sdt %>% filter(!grepl(";",Variants)) %>% group_by(Variants, ID) %>% mutate(n=n()) %>% ungroup() %>% filter(value.freq>zeroo) %>% filter(value.depth > min.depth) %>% group_by(Variants, ID) %>% mutate(N=n()) %>% mutate(r=N/n) %>% filter(r > minUniqMarkerRatio) %>% filter(N >= minUniqMarker) %>% summarize(.groups = "drop_last") -> sdt2
    sdt %>%  mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) %>% group_by(Variants, ID) %>% mutate(n=n()) %>% ungroup() %>% filter(value.freq>zeroo) %>% filter(value.depth > min.depth) %>% group_by(Variants, ID) %>% mutate(N=n()) %>% mutate(r=N/n) %>% filter(r >= minMarkerRatio) %>% filter(N >= minMarker) %>% summarize(.groups = "drop_last") -> sdt3
    uniqMarkerPerVariant <- unique(sdt2$Variants)
    markerPerVariant <- unique(sdt3$Variants)

    ## accept all uniqMarkerPerVariant (with more than minUniqMarker and minUniqMarkerRatio)
    ## accept all markerPerVariant if ancestors to accepted by uniqMarkerPerVariant
    uniqMarkerPerVariant_dealiased <- unlist(lapply(as.list(uniqMarkerPerVariant), dealias))
    markerPerVariant_dealiased <- unlist(lapply(as.list(markerPerVariant), dealias))

    specifiedLineagesTimecourse <- uniqMarkerPerVariant
    if(length(markerPerVariant) > 0){
      for (ii in 1:length(markerPerVariant)){
        if(any(grepl(markerPerVariant_dealiased[ii], uniqMarkerPerVariant_dealiased, fixed = TRUE))){
          specifiedLineagesTimecourse <- c(specifiedLineagesTimecourse, markerPerVariant[ii])
          }
        }
    }
    specifiedLineagesTimecourse <- unique(specifiedLineagesTimecourse)
    specifiedLineagesTimecourse[!(is.na(specifiedLineagesTimecourse))] -> specifiedLineagesTimecourse

    ### FROM HERE LOOP OVER TIME POINTS
    sdt %>% dplyr::select(sample_date_decimal, sample_date) %>% distinct() -> timePointsCombinations
    timePoints <- timePointsCombinations$sample_date_decimal[order(timePointsCombinations$sample_date_decimal)]
    timePoints_classic <- timePointsCombinations$sample_date[order(timePointsCombinations$sample_date_decimal)]


    if (length(timePoints) >= timeStart){
    for (t in timeStart:length(timePoints)){

      timepoint <- timePoints[t]
      timepoint_classic <- timePoints_classic[t]
      timepoint_day <- decimalDate(timepoint_classic,0)
      T <- which(timePoints_classic == timepoint_classic)
      timepoint <- timePoints[T]
      print(paste("PROGRESS:", roiname, "@", timepoint_classic, "(", signif(timepoint, digits = 10), ")"))

      lowerDiff <- ifelse(min(T) > timeLag, timeLag, min(T)-1)
      upperDiff <- ifelse((max(T)+timeLag) <= length(timePoints), timeLag, length(timePoints)-max(T))
      diffs <- (-lowerDiff:upperDiff)
      allTimepoints <- timePoints[min(T+diffs):max(T+diffs)]
      allTimepoints_classic <- timePoints_classic[min(T+diffs):max(T+diffs)]
      if(any(abs(as.numeric(timepoint_day - allTimepoints)) <= (timeLagDay/(leapYear(floor(timepoint_day)))))){
        allTimepoints_classic <- allTimepoints_classic[abs(as.numeric(timepoint_day - allTimepoints)) <= timeLagDay/(leapYear(floor(timepoint_day)))]
        allTimepoints <- allTimepoints[abs(as.numeric(timepoint_day - allTimepoints)) <= timeLagDay/(leapYear(floor(timepoint_day)))]

      } else{
        allTimepoints <- timepoint
        allTimepoints_classic <- allTimepoints_classic
      }

      sdt %>% filter(sample_date %in% c(timepoint_classic, allTimepoints_classic) ) -> ssdt

      # find lineages for which at least minUniqMarker are found in current timepoint(s)
      # and minUniqMarkerRatio of all uniq markers are found in current timepoint
      ssdt %>% filter(sample_date_decimal %in% timepoint) %>% filter(!grepl(";",Variants)) %>% group_by(Variants, ID) %>% mutate(n=n()) %>% ungroup() %>% filter(value.freq>zeroo) %>% filter(value.depth > min.depth) %>% group_by(Variants, ID) %>% mutate(N=n()) %>% mutate(r=N/n) %>% filter(r > minUniqMarkerRatio) %>% filter(N >= minUniqMarker) %>% summarize(.groups = "drop_last") -> ssdt2
      ssdt %>% filter(sample_date_decimal == timepoint) %>%  mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) %>% group_by(Variants, sample_date, sample_date_decimal) %>% mutate(n=n()) %>% ungroup() %>% filter(value.freq>zeroo) %>% filter(value.depth > min.depth) %>% group_by(Variants, sample_date, sample_date_decimal) %>% mutate(N=n()) %>% mutate(r=N/n) %>% filter(r >= minMarkerRatio) %>% filter(N >= minMarker) %>% summarize(.groups = "drop_last") -> ssdt3

      uniqMarkerPerVariant <- unique(ssdt2$Variants)
      markerPerVariant <- unique(ssdt3$Variants)
      ## accept all uniqMarkerPerVariant (with more than minUniqMarker and minUniqMarkerRatio)
      ## accept all markerPerVariant if ancestors to accepted by uniqMarkerPerVariant
      uniqMarkerPerVariant_dealiased <- unlist(lapply(as.list(uniqMarkerPerVariant), dealias))
      markerPerVariant_dealiased <- unlist(lapply(as.list(markerPerVariant), dealias))

      specifiedLineages <- uniqMarkerPerVariant

      if(length(unique(markerPerVariant, uniqMarkerPerVariant)) > 1){
        for (ii in 1:length(markerPerVariant)){
          if(any(grepl(markerPerVariant_dealiased[ii], uniqMarkerPerVariant_dealiased, fixed = TRUE))){
            specifiedLineages <- c(specifiedLineages, markerPerVariant[ii])
          }
        }
        specifiedLineages <- unique(specifiedLineages)
        specifiedLineages[!is.na(specifiedLineages)] -> specifiedLineages
      }
      rm(ssdt2, ssdt3, uniqMarkerPerVariant, markerPerVariant)

      print(paste("LOG:", timepoint_classic, roiname, paste("(",length(allTimepoints[allTimepoints %in% timepoint]), "same day;", "+", length(allTimepoints[allTimepoints %notin% timepoint]), "neighboring TP;", length(specifiedLineages), "detected lineages)")))


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
          ssdt[( rowSums(as.data.frame(ssdt)[,colnames(ssdt) %in% unique(specifiedLineages, highlightedVariants) ]) > 0 & rowSums(as.data.frame(ssdt)[,colnames(ssdt) %in% unique(specifiedLineages, highlightedVariants) ]) < length(unique(specifiedLineages, highlightedVariants))),] -> detour4log
        } else{
          ssdt[as.data.frame(ssdt)[,colnames(ssdt) %in% specifiedLineages ] > 0,] -> ssdt
          ssdt[as.data.frame(ssdt)[,colnames(ssdt) %in% unique(specifiedLineages, highlightedVariants) ] > 0,] -> detour4log
        }

        ## remove zeros and low depth
        ssdt %>% filter(value.freq > zeroo)  %>% filter(value.depth > min.depth) -> ssdt

        ## transform to avoid 1
        ssdt %>% mutate(value.freq = (value.freq * (value.depth-1) + 0.5)/value.depth) -> ssdt


        # remove "OTHERS"
        ssdt %>% filter(!grepl("other", Variants)) -> ssdt

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
          ssdt %>% filter(!grepl(";", Variants)) -> ssdt
          if ( length(unique(ssdt$Variants)) > 1 ){
            formula <- as.formula(paste("value.freq", "Variants", sep = " ~" ))
          } else{
            formula <- as.formula(paste("value.freq", "1", sep = " ~" ))
          }
        }

        ## check if the current timepoint is still represented in the data
        ssdt %>% filter(sample_date_decimal %in% timepoint) -> ssdt2
        if(dim(ssdt2)[1] == 0){
          print(paste("LOG: since no mutations found in timepoint of interest (", timepoint, ") regression will be skipped"))
          next;
        }

        ## add weight (1/(dayDiff+1)*log10(sequencingdepth)
        ssdt %>% rowwise() %>%  mutate(timeweight = 1/(abs(timePoints[t] - sample_date_decimal)*(leapYear(floor(timePoints[t]))) + 1)) %>% mutate(weight = log10(value.depth)*timeweight) -> ssdt
        detour4log %>% rowwise() %>%  mutate(timeweight = 1/(abs(timePoints[t] - sample_date_decimal)*(leapYear(floor(timePoints[t]))) + 1)) %>% mutate(weight = log10(value.depth)*timeweight) %>% filter(!grepl(";.+;", Variants)) -> detour4log

        print(data.frame(date = detour4log$sample_date, decimal = detour4log$sample_date_decimal, Tweight = round(detour4log$timeweight, digits = 2), depth = detour4log$value.depth, weight = round(detour4log$weight, digits = 2), value = round(detour4log$value.freq, digits = 2), variants = detour4log$Variants, mutation = detour4log$NUC))
        rm(detour4log)

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

        O = 1
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

        ssdt %>% ungroup() %>% filter(Variants %in% specifiedLineagesTimecourse) %>% mutate(T = unique(ssdtFit$T)) %>% mutate(fit2 = ifelse(T > 1, fit1/T, fit1)) %>% filter(sample_date_decimal == timePoints[t]) -> ssdt2
        plantFullData <- rbind(plantFullData, data.table(variant = ssdt2$Variants, LocationID =  rep(roi, length(ssdt2$fit1)), LocationName  =  rep(roiname, length(ssdt2$fit1)), sample_date = as.character(ssdt2$sample_date), value = ssdt2$fit2, marker = ssdt2$NUC, singlevalue = ssdt2$value.freq ))

        ## save norm.fitted values into global DF; set all untested variants to 0
        for (j in specifiedLineagesTimecourse){
          ssdtFit %>% filter(Variants == j) -> extractedFt
          if (dim(extractedFt)[1] > 0){
            plantFittedData <- rbind(plantFittedData, data.table(variant = rep(j, length(extractedFt$fit2)), LocationID =  rep(roi, length(extractedFt$fit2)), LocationName  =  rep(roiname, length(extractedFt$fit2)), sample_id = extractedFt$ID, sample_date =  rep(as.character(unique(ssdt2$sample_date)), length(extractedFt$fit2)), value = extractedFt$fit2 ))
          } else{
            plantFittedData <- rbind(plantFittedData, data.table(variant = j, LocationID = roi, LocationName  = roiname, sample_id = extractedFt$ID, sample_date = as.character(unique(ssdt2$sample_date)), value = 0 ))
          }
        }
        rm(formula, fit1, ssdt, ssdt2)
      } else {
        print(paste("LOG:", "NOTHIN FOR REGRESSION AT", timePoints[t], "since number of lineages detected ==", length(specifiedLineages), "; ... ignore TIMEPOINT"))


      }
    }
    }

    print(paste("     PROGRESS: start plotting ", roiname))

    if(dim(plantFittedData)[1] > 0){
        plantFittedData %>% mutate(latest = max(sample_date, na.rm = TRUE)) %>% filter(sample_date == latest) %>% filter(value > 0) %>% summarize(variant = variant, freq = value) -> sankey.dt

        plantFittedData %>% summarize(latest = max(sample_date, na.rm = TRUE)) -> sankey_date
        sankey.dt %>% mutate(freq = freq/sum(freq)) -> sankey.dt
        sankey.dt %>% group_by(variant = "B") %>% summarize(freq = 1-sum(freq)) -> sankey.dt0
        rbind(sankey.dt, sankey.dt0) -> sankey.dt
        sankey.dt %>% filter(freq > zeroo/10 | variant == "B") -> sankey.dt

        sankey.dt$variant_dealias <- unlist(lapply(as.list(sankey.dt$variant), dealias))
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

        ggplot(sankey.dtt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = label)) + geom_sankey(flow.alpha = .6, node.color = "gray30", type ='alluvial') + geom_sankey_label(size = 3, color = "white", fill = "gray40", position = position_nudge(x = 0.05, y = 0), na.rm = TRUE, type ='alluvial', hjust = 0) + theme_sankey(base_size = 16) + labs(x = NULL) + ggtitle(roiname, subtitle = sankey_date) + theme(legend.position = "none", axis.text.x = element_blank(), plot.title = element_text(hjust = 0), plot.subtitle=element_text(hjust = 0)) + scale_fill_manual(values = col2var, breaks = var2col) -> pp

        filename <- paste0(outdir, "/figs/sankey/",  paste('/klaerwerk', roi, sep="_"), ".pdf")
        ggsave(filename = filename, plot = pp, width = 8, height = 4.5)
        fwrite(as.list(c("sankey", "WWTP", roiname, filename)), file = summaryDataFile, append = TRUE, sep = "\t")
        rm(pp, filename, sankey.dt, sankey.dtt)
    }

    if(length(unique(plantFittedData$sample_date)) >= tpLimitToPlot){

      if(dim(plantFittedData)[1] >= 1){
        print(paste("     PROGRESS: plotting stackedview", roiname))

        ## print stacked area overview of all detected lineages
        plantFittedData2 <- plantFittedData
        plantFittedData2$variant <- factor(plantFittedData2$variant, levels = rev(c(highlightedVariants, unique(plantFittedData2$variant)[unique(plantFittedData2$variant) %notin% highlightedVariants])))

        plantFittedData2 %>% group_by(LocationID) %>% mutate(n = length(unique(sample_date))) %>% ungroup() %>% group_by(LocationID, LocationName, sample_date, variant) %>% summarize(value = mean(value), .groups = "drop") %>% ungroup() %>% group_by(LocationID, LocationName, sample_date) %>% mutate(T = sum(value)) %>% rowwise() %>% summarize(variant = variant, value = ifelse(T>1, value/(T+0.00001), value), .groups = "drop") -> plottng_data

        ## get color by lineage base
        #@@@@
        ColorBaseData <- data.table(variant = character(), base = character())
        for (varr in unique(plottng_data$variant)){
          foundbase = "other"
          varrr <- dealias(varr)
          for (base in baseForColorsVariants){
            if(grepl(base, varrr, fixed = TRUE)){
              foundbase = base
            }
          }
          rbind(data.table(variant = varr, base = foundbase), ColorBaseData) -> ColorBaseData
        }
        ColorBaseData$variant_dealiased <- unlist(lapply(as.list(as.character(ColorBaseData$variant)), dealias))
        ColorBaseData[order(ColorBaseData$variant_dealiased, decreasing = TRUE),] -> ColorBaseData

        ColorBaseData %>% group_by(base) %>%  mutate(n = n()) %>% arrange(base, variant_dealiased) %>% dplyr::mutate(id = cur_group_id()) %>% mutate(id = ifelse(id > 6, 6, id)) -> ColorBaseData
        ColorBaseData %>% group_by(id) %>% mutate(i = row_number()) %>% rowwise() %>% mutate(col = getColor(n, id, i)) -> ColorBaseData

        ggplot(data = plottng_data, aes(x = as.Date(sample_date), y = value, fill = variant)) -> q3
        q3 <- q3 + geom_area(position = "stack", alpha = 0.6)
        q3 <- q3 + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5)
        q3 <- q3 + facet_wrap(~LocationName, ncol = 4)
        #q3 <- q3 + scale_fill_viridis_d(alpha = 0.6, begin = .05, end = .95, option = "H", direction = +1, name="")
        q3 <- q3  + scale_fill_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "")
        q3 <- q3 + scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks = c(0,0.5,1))
        q3 <- q3 + theme_bw() + theme(legend.position="bottom", strip.text.x = element_text(size = 10), panel.grid.minor = element_blank(), panel.spacing.y = unit(0, "lines"),  strip.background = element_rect(fill="grey97"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        q3 <- q3 + scale_x_date(date_breaks = "2 month", date_labels =  "%b %y", limits = c(as.Date(NA), as.Date(latestSample$latest)))
        q3 <- q3 + ylab(paste0("Variantenanteil [1/1]") ) + xlab("")
        filename <- paste0(outdir, '/figs/stackview', paste('/klaerwerk', roi, "all", sep="_"), ".pdf")
        ggsave(filename = filename, plot = q3, width = plotWidth, height = plotHeight)
        fwrite(as.list(c("stackOverview", c( ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "all", filename))), file = summaryDataFile, append = TRUE, sep = "\t")


        ## extract color assigment in this first plot to be reused in later plots of same WW plant
        g <- ggplot_build(q3)
        GrpIdx <- g$data[[1]]["group"]
        LabelTxt <- g$plot$scales$scales[[1]]$range$range
        FillIdx <- g$data[[1]]["fill"]
        data.table(category = LabelTxt[GrpIdx$group], fill = FillIdx$fill) %>% ungroup() %>% distinct() -> colorAssignment
        rm(q3, filename, plantFittedData2, FillIdx, LabelTxt, GrpIdx, g)

        ## print faceted line plot of fitted values for all detected lineages
        plantFittedData %>% filter(sample_date > as.Date(latestSample$latest)-(2*recentEnought)) %>% group_by(variant) %>% mutate(maxValue = pmax(value)) %>% filter(maxValue > 0) %>% ggplot(aes(x = as.Date(sample_date), y = value, col = variant)) + geom_line(alpha = 0.6) + geom_point(alpha = 0.66, shape = 13) + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5) + scale_color_manual(values = colorAssignment$fill, breaks = colorAssignment$category, name = "Varianten") + scale_y_continuous(labels=scales::percent, trans = "log10") + theme_bw() + theme(legend.position="bottom", strip.background = element_rect(fill="grey97")) + scale_x_date(date_breaks = "2 weeks", date_labels =  "%b %d", limits = c(as.Date(latestSample$latest)-(2*recentEnought), as.Date(latestSample$latest))) + ylab(paste0("Variantenanteil [1/1]") ) + xlab("") -> q2

        spemut_sdt %>% filter(value.freq > 0) %>% dplyr::select("ID", "sample_date", "value.freq", "LocationID", "LocationName", "NUC") %>% filter(sample_date > as.Date(latestSample$latest)-(2*recentEnought)) -> spemut_draw2
        if(dim(spemut_draw2)[1] > 0){
          print(paste("     PROGRESS: plotting special mutations", roiname))
          left_join(x = spemut_draw2, y = soi, by = "NUC") -> spemut_draw2
          spemut_draw2 %>% rowwise() %>% mutate(marker =  paste(sep = "|", NUC, paste(Gene, AA, sep = ":"), paste0("[", Variants, "]")) ) -> spemut_draw2
          colnames(spemut_draw2)[colnames(spemut_draw2) == "Variants"] <- "variant"
          q2 + geom_point(data = spemut_draw2, aes(x = as.Date(sample_date), y = value.freq, shape = marker, fill = marker), color = "black", size = 2, alpha = .45) -> q2
          q2 + scale_shape_manual(values = c(0:25, 33:50)) -> q2
          q2 + guides(shape = guide_legend(title = "Spezial-Mutationen", ncol = 2, title.position = "top"), fill = guide_legend(title = "Spezial-Mutationen", ncol = 2, title.position = "top"), col = guide_legend(title = "Varianten", ncol = 1, title.position = "top")) -> q2

          filename <- paste0(outdir, '/figs/specialMutations', paste('/klaerwerk', roi, "all", sep="_"), ".pdf")
          #ggsave(filename = filename, plot = q2, width = plotWidth, height = plotHeight)
          fwrite(as.list(c("specialMutations", c( ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "all",filename))), file = summaryDataFile, append = TRUE, sep = "\t")
        }
        rm(q2,filename,spemut_draw2)

      }

      if(dim(filter(plantFittedData, variant %in% VoI))[1] >= 1){

        ## print faceted line plot of fitted values for all VoI
        colorAssignment %>% filter(category %in% VoI) -> colorAssignmentVoI
        plantFittedData %>% filter(variant %in% VoI) %>% filter(sample_date > as.Date(latestSample$latest)-(2*recentEnought)) %>% group_by(variant) %>% mutate(maxValue = pmax(value)) %>% filter(maxValue > 0) %>% ggplot(aes(x = as.Date(sample_date), y = value, col = variant)) + geom_line(alpha = 0.6) + geom_point(alpha = 0.66, shape = 13) + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5)  + scale_color_manual(values = colorAssignmentVoI$fill, breaks = colorAssignmentVoI$category, name = "Varianten") + scale_y_continuous(labels=scales::percent, trans = "log10") + theme_bw() + theme(legend.position="bottom", strip.background = element_rect(fill="grey97")) + scale_x_date(date_breaks = "2 weeks", date_labels =  "%b %d", limits = c(as.Date(latestSample$latest)-(2*recentEnought), as.Date(latestSample$latest))) + ylab(paste0("Variantenanteil [1/1]") ) + xlab("") -> q1

        spemut_sdt %>% filter(value.freq > 0) %>% dplyr::select("ID", "sample_date", "value.freq", "LocationID", "LocationName", "NUC")  %>% filter(sample_date > as.Date(latestSample$latest)-(2*recentEnought)) -> spemut_draw1
        if(dim(spemut_draw1)[1] > 0){
          print(paste("     PROGRESS: plotting special mutations VoI", roiname))
          left_join(x = spemut_draw1, y = soi, by = "NUC") -> spemut_draw1
          spemut_draw1 %>% rowwise() %>% mutate(marker = paste(sep = "|", NUC, paste(Gene, AA, sep = ":"), paste0("[", Variants, "]")) ) -> spemut_draw1
          colnames(spemut_draw1)[colnames(spemut_draw1) == "Variants"] <- "variant"
          q1 + geom_point(data = spemut_draw1, aes(x = as.Date(sample_date), y = value.freq, shape = marker, fill = marker), color = "black", size = 2, alpha = .45) -> q1
          q1 + scale_shape_manual(values = c(0:25, 33:50)) -> q1
          q1 + guides(shape = guide_legend(title = "Spezial-Mutationen", ncol = 2, title.position = "top"), fill = guide_legend(title = "Spezial-Mutationen", ncol = 2, title.position = "top"), col = guide_legend(title = "Varianten", ncol = 1, title.position = "top")) -> q1


          filename <- paste0(outdir, '/figs/specialMutations', paste('/klaerwerk', roi, "VoI", sep="_"), ".pdf")
          ggsave(filename = filename, plot = q1, width = plotWidth, height = plotHeight*1.5)
          fwrite(as.list(c("specialMutations", c( ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "VoI", filename))), file = summaryDataFile, append = TRUE, sep = "\t")
        }
        rm(q1, filename, colorAssignmentVoI, spemut_draw1)
      }

      if(dim(filter(plantFullData, variant %in% VoI))[1] >= 1){
        print(paste("     PROGRESS: plotting variantDetails VoI", roiname))
        ## print faceted line plot of fitted values plus point plot of measured AF for all VoIs
        colorAssignment %>% filter(category %in% VoI) -> colorAssignmentVoI
        plantFullData %>% filter(variant %in% VoI) %>% ggplot()  + geom_line(data = subset(plantFittedData, variant %in% VoI), alpha = 0.6, size = 2, aes(x = as.Date(sample_date), y = value, col = variant)) + geom_jitter(aes(x = as.Date(sample_date), y = singlevalue), alpha = 0.33, size = 1.5, width = 0.33, height = 0) + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5) + facet_wrap(~variant)  + scale_color_manual(values = colorAssignmentVoI$fill, breaks = colorAssignmentVoI$category, name = "") + scale_y_continuous(labels=scales::percent, trans = "log10") + theme_bw() + theme(legend.position="bottom", strip.background = element_rect(fill="grey97"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_x_date(date_breaks = "2 month", date_labels =  "%b %y", limits = c(as.Date(NA), as.Date(latestSample$latest))) + ylab(paste0("Variantenanteil [1/1]") ) + xlab("") -> p1

        filename <- paste0(outdir, '/figs/detail', paste('/klaerwerk', roi, "VoI", sep="_"), ".pdf")
        ggsave(filename = filename, plot = p1, width = plotWidth, height = plotHeight)
        fwrite(as.list(c("variantDetail", c( ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "VoI", filename))), file = summaryDataFile, append = TRUE, sep = "\t")

        rm(p1, filename, colorAssignmentVoI)
      }

      if(dim(plantFullData)[1] >= 1){
          print(paste("     PROGRESS: plotting variant Details", roiname))
        ## print faceted line plot of fitted values plus point plot of measured AF for all lineages
        plantFullData %>% ggplot() +  geom_line(data = plantFittedData, alpha = 0.6, size = 2, aes(x = as.Date(sample_date), y = value, col = variant)) + geom_jitter(aes(x = as.Date(sample_date), y = singlevalue), alpha = 0.33, size = 1.5, width = 0.33, height = 0) + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5) + facet_wrap(~variant) + scale_color_manual(values = colorAssignment$fill, breaks = colorAssignment$category, name = "") + scale_y_continuous(labels=scales::percent, trans = "log10") + theme_bw() + theme(legend.position="bottom", strip.background = element_rect(fill="grey97"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_x_date(date_breaks = "3 month", date_labels =  "%b", limits = c(as.Date(NA), as.Date(latestSample$latest))) + ylab(paste0("Variantenanteil [1/1]") ) + xlab("") -> p2

        filename <- paste0(outdir, '/figs/detail', paste('/klaerwerk', roi, "all", sep="_"), ".pdf")
        ggsave(filename = filename, plot = p2, width = plotWidth, height = plotHeight*1.5)
        fwrite(as.list(c("variantDetail", c(ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "all", filename))), file = summaryDataFile, append = TRUE, sep = "\t")
        rm(p2, filename)
      }
    }

    globalFittedData <- rbind(plantFittedData, globalFittedData)
    globalFullData <- rbind(plantFullData, globalFullData)

    rm(plantFittedData, plantFullData, colorAssignment)
}

print(paste("PROGRESS: loop over WWTP finished"))
## dump global data into output file
print(paste("PROGRESS: writing result tables"))
fwrite(globalFittedData, file = paste0(outdir, "/globalFittedData.csv"), sep = "\t")
fwrite(globalFullData, file = paste0(outdir, "/globalFullData.csv"), sep = "\t")
fwrite(globalFullSoiData, file = paste0(outdir, "/globalSpecialmutData.csv"), sep = "\t")

## print overview of all VoI and all plants
print(paste("PROGRESS: start plotting overviews"))




##

if(dim(globalFittedData)[1] > 0){
    print(paste("     PROGRESS: plotting overview Sankey plot"))

    ## make complete plot
    globalFittedData %>% filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date) %>% group_by(LocationID) %>% mutate(latest = max(sample_date, na.rm = TRUE)) %>% filter(sample_date == latest) %>% summarize(variant = variant, freq = value) -> sankey_all.dt
    unique(sankey_all.dt$variant) -> var2col
    unlist(lapply(as.list(var2col), dealias)) -> dealiasvar2col
    sort(dealiasvar2col) -> var2col
    viridis_pal(alpha = 1, begin = 0.05, end = .95, direction = 1, option = "D")(length(var2col)) -> col2var

    if(dim(sankey_all.dt)[1] > 0){

        globalFittedData %>% filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date) %>% group_by(LocationID) %>% mutate(latest = max(sample_date, na.rm = TRUE)) %>% filter(sample_date == latest) %>% ungroup() %>% summarize(earliest = min(sample_date, na.rm = TRUE), latest = max(sample_date, na.rm = TRUE)) -> sankey_date

        left_join(x = sankey_all.dt, y = (metaDT %>% dplyr::select(LocationID, connected_people) %>% distinct())) %>% group_by(variant) %>% summarize(freq = weighted.mean(freq, w = connected_people, na.rm = TRUE)) %>% filter(freq > 0) -> sankey.dt

        sankey.dt %>% mutate(freq = freq/sum(freq)) -> sankey.dt
        sankey.dt %>% group_by(variant = "B") %>% summarize(freq = 1-sum(freq)) -> sankey.dt0

        rbind(sankey.dt, sankey.dt0) -> sankey.dt
        sankey.dt %>% filter(freq > zeroo/10 | variant == "B") -> sankey.dt

        sankey.dt$variant_dealias <- unlist(lapply(as.list(sankey.dt$variant), dealias))
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

        ggplot(sankey.dtt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = label)) + geom_sankey(flow.alpha = .6, node.color = "gray30", type ='alluvial') + geom_sankey_label(size = 3, color = "white", fill = "gray40", position = position_nudge(x = 0.05, y = 0), na.rm = TRUE, type ='alluvial', hjust = 0) + theme_sankey(base_size = 16) + labs(x = NULL) + ggtitle("Gewichtetes Mittel: Österreich", subtitle = paste( sankey_date$earliest , "bis", sankey_date$latest)) + theme(legend.position = "none", axis.text.x = element_blank(), plot.title = element_text(hjust = 0), plot.subtitle=element_text(hjust = 0)) + scale_fill_manual(values = col2var, breaks = var2col) -> pp

        filename <- paste0(outdir, "/figs/sankey/Overview_",  "Austria", ".pdf")
        ggsave(filename = filename, plot = pp)
        fwrite(as.list(c("sankey", "Overview", "Austria", filename)), file = summaryDataFile, append = TRUE, sep = "\t")
        rm(pp, filename, sankey.dt, sankey.dtt)
    }

    ## make plot per state
    left_join(x = sankey_all.dt, y = metaDT, by = "LocationID") %>% dplyr::select(LocationID, variant, freq, state) %>% distinct() -> sankey_allState.dt
    sankey_state_plot_list <- list()
    for (stateoi in unique(sankey_allState.dt$state)){
      sankey_allState.dt %>% filter(state == stateoi) -> sankey_state.dt
      if(dim(sankey_state.dt)[1] > 0){

          globalFittedData %>% filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date) %>% filter(LocationID %in% unique(sankey_state.dt$LocationID)) %>% group_by(LocationID) %>% mutate(latest = max(sample_date, na.rm = TRUE)) %>% filter(sample_date == latest) %>% ungroup() %>% summarize(earliest = min(sample_date, na.rm = TRUE), latest = max(sample_date, na.rm = TRUE)) -> sankey_date

          left_join(x = sankey_state.dt, y = (metaDT %>% dplyr::select(LocationID, connected_people) %>% distinct())) %>% group_by(variant) %>% summarize(freq = weighted.mean(freq, w = connected_people, na.rm = TRUE)) %>% filter(freq > 0) -> sankey.dt

          sankey.dt %>% mutate(freq = freq/sum(freq)) -> sankey.dt
          sankey.dt %>% group_by(variant = "B") %>% summarize(freq = 1-sum(freq)) -> sankey.dt0

          rbind(sankey.dt, sankey.dt0) -> sankey.dt
          sankey.dt %>% filter(freq > zeroo/10 | variant == "B") -> sankey.dt

          sankey.dt$variant_dealias <- unlist(lapply(as.list(sankey.dt$variant), dealias))
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

          ggplot(sankey.dtt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = label)) + geom_sankey(flow.alpha = .6, node.color = "gray30", type ='alluvial') + geom_sankey_label(size = 3, color = "white", fill = "gray40", position = position_nudge(x = 0.05, y = 0), na.rm = TRUE, type ='alluvial', hjust = 0) + theme_sankey(base_size = 16) + labs(x = NULL) + ggtitle(paste0("Gewichtetes Mittel: ", stateoi), subtitle = paste( sankey_date$earliest , "bis", sankey_date$latest)) + theme(legend.position = "none", axis.text.x = element_blank(), plot.title = element_text(hjust = 0), plot.subtitle=element_text(hjust = 0)) + scale_fill_manual(values = col2var, breaks = var2col) -> pp
          sankey_state_plot_list[[length(sankey_state_plot_list)+1]] <- pp

          filename <- paste0(outdir, "/figs/sankey/Overview_",  gsub(" ", "_", stateoi), ".pdf")
          ggsave(filename = filename, plot = pp)
          fwrite(as.list(c("sankey", "Overview", stateoi, filename)), file = summaryDataFile, append = TRUE, sep = "\t")
          rm(pp, filename, sankey_state.dt, sankey.dt, sankey.dtt)
      }
      plot_grid(plotlist = sankey_state_plot_list, ncol = 3) -> ppp
      filename <- paste0(outdir, "/figs/sankey/Overview_byState.pdf")
      ggsave(filename = filename, plot = ppp, width = 20, height = 20)
      fwrite(as.list(c("sankey", "Overview", "allStates", filename)), file = summaryDataFile, append = TRUE, sep = "\t")
      rm(ppp, filename, sankey_state.dt, sankey.dt, sankey.dtt)


    }
}

##


globalFittedData %>% group_by(LocationID) %>% mutate(n = length(unique(sample_date))) %>% filter(n > tpLimitToPlot*2) %>% filter(variant %in% VoI) %>% ungroup() %>% group_by(LocationID, LocationName, sample_date, variant) %>% summarize(value = mean(value), .groups = "drop") -> globalFittedData2
if (dim(globalFittedData2)[1] >= tpLimitToPlot){
  print(paste("     PROGRESS: plotting overview VoI"))
  ggplot(data = globalFittedData2, aes(x = as.Date(sample_date), y = value, fill = variant)) + geom_area(position = "stack", alpha = 0.7)  + facet_wrap(~LocationName) + scale_fill_viridis_d(alpha = 0.6, begin = .05, end = .95, option = "H", direction = +1, name="") + scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks = c(0,0.5,1)) + theme_minimal() + theme(legend.position="bottom", strip.text.x = element_text(size = 4.5), panel.grid.minor = element_blank(), panel.spacing.y = unit(0, "lines"), legend.direction="horizontal") + guides(fill = guide_legend(title = "", ncol = 10)) + scale_x_date(date_breaks = "2 month", date_labels =  "%b") + ylab(paste0("Variantenanteil [1/1]") ) + xlab("") -> r2
  filename <- paste0(outdir, '/figs/fullview', paste('/klaerwerke', "VoI", sep="_"), ".pdf")
  ggsave(filename = filename, plot = r2, width = plotWidth*1.5, height = plotHeight*1.5)
  fwrite(as.list(c("fullOverview", c(roiname, "VoI", filename))), file = summaryDataFile, append = TRUE, sep = "\t")
  rm(r2,filename)
}
rm(globalFittedData2)

## print overview of all variants and all plants
globalFittedData$variant <- factor(globalFittedData$variant, levels = rev(c(highlightedVariants, unique(globalFittedData$variant)[unique(globalFittedData$variant) %notin% highlightedVariants])))
globalFittedData %>% group_by(LocationID) %>% mutate(n = length(unique(sample_date))) %>% filter(n > tpLimitToPlot*2) %>% ungroup() %>% group_by(LocationID, LocationName, sample_date, variant) %>% summarize(value = mean(value), .groups = "drop") -> globalFittedData2
if (dim(globalFittedData2)[1] >= tpLimitToPlot){
  print(paste("     PROGRESS: plotting overview"))
  ggplot(data = globalFittedData2, aes(x = as.Date(sample_date), y = value, fill = variant)) + geom_area(position = "stack", alpha = 0.7)  + facet_wrap(~LocationName) + scale_fill_viridis_d(alpha = 0.6, begin = .05, end = .95, option = "H", direction = +1, name="") + scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks = c(0,0.5,1)) + theme_minimal() + theme(legend.position="bottom", strip.text.x = element_text(size = 4.5), panel.grid.minor = element_blank(), panel.spacing.y = unit(0, "lines"), legend.direction="horizontal") + guides(fill = guide_legend(title = "", ncol = 10)) + scale_x_date(date_breaks = "2 month", date_labels =  "%b") + ylab(paste0("Variantenanteil [1/1]") ) + xlab("") -> r1
  filename <- paste0(outdir, '/figs/fullview', paste('/klaerwerke', "all", sep="_"), ".pdf")
  ggsave(filename = filename, plot = r1, width = plotWidth*1.5, height = plotHeight*1.5)
  fwrite(as.list(c("fullOverview", c(roiname, "all", filename))), file = summaryDataFile, append = TRUE, sep = "\t")
  rm(r1,filename)
}
rm(globalFittedData2)


# Generate map of each VoI detected recently

print(paste("PROGRESS: start plotting maps"))
for (voi in VoI){
  print(paste("     PROGRESS: considering", voi))

  globalFittedData %>% filter(variant == voi) %>% filter(! is.na(sample_id) ) -> mapdata
  metaDT$sample_date <- as.character(metaDT$sample_date)
  left_join(x = mapdata, y = metaDT, by = c("sample_id" = "BSF_sample_name", "LocationName", "LocationID", "sample_date")) -> mapdata

  if (dim(mapdata)[1] > 0){
    mapdata %>% group_by(LocationID) %>% mutate(latest = max(as.Date(sample_date))) %>% filter(latest == sample_date) %>% filter(BSF_run == last_BSF_run_id) %>% filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date) -> mapdata
    mapdata %>% filter(value > 0) -> mapdata

    if (length(mapdata$value) > 0 ){
      print(paste("     PROGRESS: plotting map for", voi))
      print(data.frame(location=mapdata$LocationName, dates = mapdata$sample_date, values = mapdata$value))

      s <- ggplot()
      s <- s + geom_sf(data = World, fill = "grey95")
      s <- s + geom_sf(data = Country, fill = "antiquewhite")
      s <- s + theme_minimal()
      s <- s + coord_sf(ylim = mapMargines[c(1,3)], xlim = mapMargines[c(2,4)], expand = FALSE)
      s <- s + annotation_scale(location = "bl")
      s <- s + annotation_north_arrow(location = "tl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering)
      s <- s + geom_point(data=subset(mapdata, status == "pass"), aes(y=dcpLatitude, x=dcpLongitude, size = as.numeric(connected_people)), alpha = 1, shape = 3)
      s <- s + geom_point(data=mapdata, aes(y=dcpLatitude, x=dcpLongitude, size = as.numeric(connected_people), col = as.numeric(value)), alpha = 0.8)
      s <- s + facet_wrap(~variant, nrow = 2)
      s <- s + theme(axis.text = element_blank(), legend.direction = "vertical", legend.box = "horizontal", legend.position = "bottom")
      s <- s + guides(size = guide_legend(title = "Population", nrow = 2), col = guide_colourbar(title = paste0("Anteil ", voi), direction = "horizontal", title.position = "top"))
      s <- s + scale_color_viridis_c(limits = c(0,1), direction = 1, begin = .05, end = .95, option = "C")
      s <- s + scale_size(range = c(2, 6), labels = scales::comma, breaks = c(50000, 100000, 500000, 1000000))

      filename <- paste0(outdir, '/figs/maps/', voi, ".pdf")
      ggsave(filename = filename, plot = s)
      fwrite(as.list(c("mapplot", voi, filename)), file = summaryDataFile, append = TRUE, sep = "\t")
      rm(s, filename, mapdata)
    }
  }
}

timestamp()
