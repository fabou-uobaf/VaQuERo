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
  make_option(c("--alias_fh"), type="character", default="https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json",
              help="Json file linking pango lineage nomenclature to predecessor. Default value requires internet connection [default %default]", metavar="character"),
  make_option(c("--plotwidth"), type="double", default=8,
              help="Base size of plot width [default %default]", metavar="character"),
  make_option(c("--plotheight"), type="double", default=4.5,
              help="Base size of plot height [default %default]", metavar="character"),
  make_option(c("--ninconsens"), type="double", default="0.4",
              help="Minimal fraction of genome covered by reads to be considered (0-1) [default %default]", metavar="character"),
  make_option(c("--zero"), type="double", default=0.02,
              help="Minimal allele frequency to be considered [default %default]", metavar="double"),
  make_option(c("--alphaprime"), type="double", default=3.2,
              help="One of the precomputed beta distribution parameter used to model the mutation frequency dinstribution [default %default]", metavar="double"),
  make_option(c("--depth"), type="integer", default=75,
              help="Minimal depth at mutation locus to be considered [default %default]", metavar="character"),
  make_option(c("--recent"), type="integer", default=999,
              help="How old (in days) most recent sample might be to be still considered in overview maps [default %default]", metavar="character"),
  make_option(c("--plottp"), type="integer", default=3,
              help="Produce timecourse plots only if more than this timepoints are available [default %default]", metavar="character"),
  make_option(c("--minuniqmark"), type="integer", default=1,
              help="Minimal absolute number of uniq markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--minuniqmarkfrac"), type="double", default=0.4,
              help="Minimal fraction of uniq markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--mininfofrac"), type="double", default=0.6,
              help="Minimal fraction of markers that variant is considered detected [default %default]", metavar="character"),
  make_option(c("--addUniqZeros"), type="logical", default="FALSE", action = "store_true",
            help="Should unobserved uniq marker added with AF=0 [default %default]", metavar="character"),
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
  make_option(c("--debug"), type="logical", default="FALSE", action = "store_true",
              help="Toggle to use run with provided example data [default %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



#####################################################
####  parameter setting for interactive debugging
if(opt$debug){
    opt$dir = "output-variants"
    opt$metadata = "data/metaData_general.csv"
    opt$data="data/mutationData_DB_NationMonitoringSites.tsv.gz"
    opt$marker="VaQuERo/resources/mutations_list_grouped_pango_codonPhased_2025-05-05.csv"
    opt$pmarker="VaQuERo/resources/mutations_problematic_2024-10-04.csv"
    opt$smarker="VaQuERo/resources/mutations_special_2022-12-21.csv"
    opt$zero=0.01
    opt$ninconsens=0.4
    opt$depth=50
    opt$minuniqmark=1
    opt$minuniqmarkfrac=0.4
    opt$mininfofrac=0.6
    opt$addUniqZeros=TRUE
    opt$alphaprime=2.2
    opt$smoothingsamples=2
    opt$smoothingtime=15
    opt$voi="XEC,LP.8.1,NB.1.8.1,XFG"
    opt$highlight="XEC,KP.3,KP.2,BA.2.86"
    opt$colorBase="XEC,KP.3,BA.2.86"
    opt$recent <- 23
    print("Warning: command line option overwritten")
}

#####################################################

## define functions
options(warn=-1)


## loading external function library, expected to be in same dir as execution script
sourceFileBase = "^VaQueR_functions.r$"
sourceDir <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
        } else {
                # 'source'd via R console
                return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
        }
}
if(!interactive()){
  sourceFile <- list.files(
    sourceDir(),
    pattern = sourceFileBase,
    recursive = TRUE,
    full.names = TRUE
  )
} else{
  sourceFile <- c()
}

if(length(sourceFile) == 1 && file.exists(sourceFile)){
  print(paste("LOG: loading function source file", sourceFile))
  opt$source = sourceFile
  source(sourceFile)
} else{
  sourceFile <- list.files(
    ".",
    pattern = sourceFileBase,
    recursive = TRUE,
    full.names = TRUE
  )
  if(length(sourceFile) == 1 && file.exists(sourceFile)){
    print(paste("LOG: loading function source file", sourceFile))
    opt$source = sourceFile
    source(sourceFile)
  } else{
    print(paste("ERROR: source file", sourceFileBase, "not found. Please double check if it is in the same directory as the analysis script VaQuERo_v2.R"))
    exit()
  }
}

## print parameter to Log

writeLines("##~LOG~PARAMETERS~####################")
str(opt)
writeLines("##~LOG~PARAMETERS~####################")
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
alphaprime <- opt$alphaprime



## create directory to write plots
writeLines(paste0("PROGRESS: create directory "))
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
aliases <- fromJSON(file =opt$alias_fh)
if(!exists("aliases")){
  writeLines(paste("Error(s): lineage alias json file", aliases_fh ,"not available. Please check internet connection or set it manually local copy via --alias_fh"))
  stop()
}
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
minInfoRatio <- opt$mininfofrac # minmal fraction of information from sensitive markers that variant is considered detected (0.6)
timeLag <- opt$smoothingsamples # number of previous timepoints use for smoothing
timeStart <- 1 # set to 1 if all time points should be considered; 2 if first should not be considered
timeLagDay <- opt$smoothingtime # previous timepoints for smoothing are ignored if more days before
VoI                 <- grep(";|,", unlist(str_split(opt$voi, pattern=c(";",","))), value = TRUE, invert = TRUE)
highlightedVariants <- unlist(str_split(opt$highlight, pattern=c(";",","))[[2]])
baseForColorsVariants <- unlist(str_split(opt$colorBase, pattern=c(";",","))[[2]])
baseForColorsVariants <- unlist(lapply(as.list(baseForColorsVariants), dealias))
colorSets <- c("Blues", "Greens", "Oranges", "Purples", "Reds", "Greys")


## define global variables to fill
globalFittedDataList  <- list()
globalFullDataList  <- list()
globalFullSoiDataList  <- list()

globalFittedData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_id = character(), sample_date = character(), value = numeric() )
globalFullData <- data.table(variant = character(), all_variants = character(), LocationID = character(), LocationName  = character(), sample_date = character(), value = numeric(), marker = character(), singlevalue = numeric() )
globalFullSoiData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_date = character(), marker = character(), value = numeric())


## read mutations of interest from file
writeLines(paste0("PROGRESS: read mutation files "))
moi <- fread(file = markermutationFile)
unite(moi, NUC, c(4,3,5), ALT, sep = "", remove = FALSE) -> moi
moi %>% mutate(Variants = gsub("other;?", "", Variants)) %>% mutate(Variants = gsub(";;", ";", Variants)) -> moi   ### just there to fix GISAID issue with BQ.1.1 vs BQ.1
moi %>% dplyr::select("Variants") %>% group_by(Variants) %>% summarize(n = n(), .groups = "keep")  -> moi_markerCombination_count
moi %>% filter(!grepl(";", Variants)) %>% dplyr::select("Variants") %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(cols = c(Variants)) %>% group_by(Variants) %>% summarize(n = n(), .groups = "keep") -> moi_uniq_marker_count

## counter number of mutations for each variant and calculate information content for each variant
variant_counter <- moi %>%
                      dplyr::select("Variants") %>%
                      mutate(Variants = strsplit(as.character(Variants), ";")) %>%
                      unnest(cols = c(Variants)) %>%
                      pull(Variants) %>%
                      unique() %>%
                      length()

moi_marker_information <- moi %>%
                      dplyr::select("Variants", "NUC") %>%
                      mutate(Variants = strsplit(as.character(Variants), ";")) %>%
                      unnest(cols = c(Variants)) %>%
                      group_by(NUC) %>%
                      mutate(i = -1*log(n()/variant_counter)) %>%
                      dplyr::select("NUC", "i") %>%
                      distinct()

moi_variant_information <- moi %>%
                      dplyr::select("Variants", "NUC") %>%
                      mutate(Variants = strsplit(as.character(Variants), ";")) %>%
                      unnest(cols = c(Variants)) %>%
                      group_by(NUC) %>%
                      mutate(i = -1*log(n()/variant_counter)) %>%
                      group_by(Variants) %>%
                      summarize(n = n(), i = sum(i), .groups = "keep")


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
writeLines(paste0("PROGRESS: read and process meta data "))
metaDT       <- fread(file = opt$metadata)
unique(metaDT$BSF_sample_name) -> sampleoi

## check if all voi can be detected in principle given used parameters
moi %>% dplyr::select("Variants") %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest() %>% distinct() %>% filter(Variants != "other")-> variants_total
moi %>% filter(!grepl(";", Variants)) %>% group_by(Variants) %>% summarize(n = n(),  .groups = "keep") %>% filter(n >= minUniqMarker) -> variants_w_enough_uniq

variants_total$Variants[variants_total$Variants %notin% variants_w_enough_uniq$Variants] -> variants_futile
if(length(variants_futile) > 0){
  for (c in seq_along(variants_futile)){
    writeLines(paste("Warning(s):", variants_futile[c], "has less than", minUniqMarker, "unique markers defined."))
  }
}

# check for same date, same location issue
# introduce artifical decimal date
metaDT %>% group_by(LocationID, sample_date) %>% mutate(n = rank(BSF_sample_name)-1)  %>% ungroup() %>% rowwise() %>% mutate(sample_date_decimal = decimalDate(sample_date, n)) %>% dplyr::select(-n) -> metaDT

# define which locations are used for the report
metaDT %>% dplyr::select("LocationID") %>% distinct() -> locationsReportedOn


# Generate map of all locations sequenced in current run

## determine last sequencing batch ID
metaDT %>% filter(!is.na(BSF_run)) %>%
    filter(!is.na(BSF_start_date)) %>%
    group_by(BSF_run, BSF_start_date) %>%
    summarize(.groups = "keep") %>%
    ungroup() %>%
    mutate(BSF_start_date_max = max(as.Date(BSF_start_date), na.rm = TRUE)) %>%
    filter(as.Date(BSF_start_date) == max(as.Date(BSF_start_date))) %>%
    dplyr::select("BSF_run") %>%
    as.character() -> last_BSF_run_id
writeLines(paste("LOG: last_BSF_run_id", last_BSF_run_id))

## modify status of samples in last run based on number N in consensus
metaDT %>% filter(BSF_run == last_BSF_run_id) %>% filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date) -> mapSeqDT
mapSeqDT %>% mutate(status = ifelse(is.na(status), "fail", status)) -> mapSeqDT
mapSeqDT %>% mutate(status = ifelse(N_in_Consensus > N_in_Consensus_filter & status == "pass", "fail", status)) -> mapSeqDT
mapSeqDT %>% mutate(status = ifelse(status == "fail" & N_in_Consensus < N_in_Consensus_detection_filter, "detected", status)) -> mapSeqDT

## remove samples with "in_run" flag set in status
if(any(grepl("in_run", mapSeqDT$status))){
  in_run_count <- sum(grepl("in_run", mapSeqDT$status))
  writeLines(paste0("WARNING: status 'in_run' found ", in_run_count," times; set to 'fail'"))
  mapSeqDT %>% mutate(status = ifelse(status == "in_run", "fail", status)) -> mapSeqDT
}

# print log how many samples from last run passed filter
writeLines(paste0("LOG: print STATUS counts: "))
print(table(mapSeqDT$status))

## generate empty map
writeLines(paste0("PROGRESS: print STATUS map"))
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
ggsave(filename = filename, plot = s, width = 7, height = 7)
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
  writeLines(paste0("PROGRESS: read AF data (deprecated file format!) "))
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
  writeLines(paste0("PROGRESS: read AF data "))
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

  ## remove all positions which are not mutations of interest
  sewage_samps.dt %>% filter(NUC %in% allmut) -> sewage_samps.dt

  ## remove all samples which are not specified in metadata
  sewage_samps.dt %>% filter(ID %in% sampleoi) -> sewage_samps.dt

  rm(sewage_samps)
}


## add location to sewage_samps.dt
sewage_samps.dt %>% mutate(RNA_ID_int = gsub("_S\\d+$","", ID)) -> sewage_samps.dt
sewage_samps.dt %>% mutate(RNA_ID_int = gsub("_\\d$","", RNA_ID_int))  -> sewage_samps.dt

if(any(is.na(sewage_samps.dt$LocationID))){
  writeLines(paste("WARNING: samples without assigned location were found and removed."))
  sewage_samps.dt %>% filter(!is.na(LocationID)) -> sewage_samps.dt
}

metaDT %>% filter(! is.na(LocationID) ) %>% dplyr::select("RNA_ID_int", "LocationID", "LocationName") -> sample_location
left_join(x = sewage_samps.dt, y = sample_location, by = "RNA_ID_int", multiple = "all") -> sewage_samps.dt

## remove samples which are not include_in_report == TRUE
metaDT %>% filter(is.na(include_in_report) | include_in_report == TRUE ) %>% dplyr::select("RNA_ID_int", "BSF_sample_name") -> sample_includeInReport
sewage_samps.dt %>% filter(ID %in% sample_includeInReport$BSF_sample_name) -> sewage_samps.dt

## remove samples with >N_in_Consensus_filter N
metaDT  %>% filter(N_in_Consensus < N_in_Consensus_filter) %>% dplyr::select("RNA_ID_int", "BSF_sample_name") -> passed_samples
sewage_samps.dt <- sewage_samps.dt[sewage_samps.dt$ID %in% passed_samples$BSF_sample_name,]

## add variant of interests to table
left_join(x=sewage_samps.dt, y=moi, by = c("POS"="Postion", "REF", "ALT", "NUC"), multiple = "all") -> sewage_samps.dt

## add sampling date
metaDT %>% dplyr::select("RNA_ID_int", "sample_date", "sample_date_decimal") -> sample_dates
left_join(x=sewage_samps.dt, y=sample_dates, by = "RNA_ID_int", multiple = "all") -> sewage_samps.dt

## get most recent sampling date from last Run
metaDT %>% filter(BSF_start_date == sort(metaDT$BSF_start_date)[length(sort(metaDT$BSF_start_date))]) %>% dplyr::select("BSF_run", "RNA_ID_int", "BSF_sample_name", "sample_date") -> RNA_ID_int_currentRun

metaDT %>% ungroup() %>% summarize(latest = max(as.Date(sample_date)), .groups = "keep") -> latestSample
RNA_ID_int_currentRun %>% ungroup() %>% summarize(earliest = min(as.Date(sample_date)), .groups = "keep") -> earliestSample

writeLines(paste("LOG: current run ID:", unique(RNA_ID_int_currentRun$BSF_run)))
writeLines(paste("LOG: earliest sample in current run:", earliestSample$earliest))
writeLines(paste("LOG: latest sample in current run:", latestSample$latest))

### FROM HERE LOOP OVER EACH SEWAGE PLANTS
writeLines(paste("PROGRESS: start to loop over WWTP"))
timestamp()
# r <-  grep("ATTP_9", unique(sewage_samps.dt$LocationID))

for (r in 1:length(unique(sewage_samps.dt$LocationID))) {
    roi = unique(sewage_samps.dt$LocationID)[r]
    roiname = sewage_samps.dt %>% filter(LocationID == roi ) %>% pull(LocationName) %>% unique()
    writeLines(paste("PROGRESS: start processing", r, roiname, roi))

    ### define data collector
    plantFittedData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_id = character(), sample_date = character(), value = numeric() )
    plantFullData <- data.table(variant = character(), all_variants = character(), LocationID = character(), LocationName  = character(), sample_date = character(), value = numeric(), marker = character(), singlevalue = numeric() )
    plantFullSoiData <- data.table(variant = character(), LocationID = character(), LocationName  = character(), sample_date = character(), marker = character(), value = numeric())

    ### filter data set for current location
    sewage_samps.dt %>% filter(LocationID == roi) %>% dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC") %>% filter(NUC %notin% exsoimut) -> sdt
    sewage_samps.dt %>% filter(LocationID == roi) %>% dplyr::select("Variants", "ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC") %>% filter(NUC %in% soimut) -> spemut_sdt

    ### collect frequency of special mutations
    left_join(x = spemut_sdt, y = soi, by = "NUC", multiple = "all") %>% rowwise() %>% mutate(marker = paste(sep = "|", NUC, paste(Gene, AA, sep = ":"), paste0("[", Variants.y, "]")) ) %>% dplyr::select(Variants.y, LocationID, LocationName, sample_date, marker, value.freq) -> plantFullSoiData
    colnames(plantFullSoiData) <- c("variant", "LocationID", "LocationName", "sample_date", "marker", "value")
    rbind(plantFullSoiData, globalFullSoiData) -> globalFullSoiData

    ## count how many sample from this plant were in the last BSF run
    metaDT %>% filter(BSF_sample_name %in% sdt$ID) %>% filter(BSF_run %in% last_BSF_run_id) %>% dplyr::select(BSF_run) %>% group_by(BSF_run) %>% summarize(n = n(), .groups = "keep") -> count_last_BSF_run_id

    specifiedLineagesTimecourse <- moi_variant_information$Variants

    ### FROM HERE LOOP OVER TIME POINTS
    sdt %>% dplyr::select(sample_date_decimal, sample_date) %>% distinct() -> timePointsCombinations
    timePoints <- timePointsCombinations$sample_date_decimal[order(timePointsCombinations$sample_date_decimal)]
    timePoints_classic <- timePointsCombinations$sample_date[order(timePointsCombinations$sample_date_decimal)]

    if (length(timePoints) >= timeStart){
          t <- timeStart
          catchCounter = 0
          while(t <= length(timePoints)){
            checkSuccess<-tryCatch( {
                timepoint <- timePoints[t]
                timepoint_classic <- timePoints_classic[t]
                ref_timepoint_classic <- timepoint_classic
                ref_timepoint         <- timepoint
                timepoint_day <- decimalDate(timepoint_classic,0)
                T <- which(timePoints_classic == timepoint_classic)
                timepoint <- timePoints[T]
                writeLines(paste("PROGRESS:", roiname, paste0("(tp: ", t, ")"), "@", paste(ref_timepoint_classic, collapse=", "), "(", paste(signif(timepoint, digits = 10), collapse=", "), ")"))

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
                  allTimepoints_classic <- rep(timepoint_classic, times = length(timepoint))
                }

                sdt %>% filter(sample_date %in% c(timepoint_classic, allTimepoints_classic) ) -> ssdt

                # detect lineage in current center timepoint
                adaptivMinInfoRatio <- 1-((1-minInfoRatio)^(1.1))
                specifiedLineages <- detect_lineages(DT_ = ssdt, timepoint_ = timepoint, minInfoRatio_ = adaptivMinInfoRatio)

                # deprecated; use neigboring timepoints individially and merge resulting list
                if(0){
                  # check if current timepoint is squeezed between two timepoint
                  # detect lineage in neighboring timepoints
                  # add lineages detected in any to specifiedLineages
                  # or if last time point, check which lineages are detecten in any of the previous two timePoints
                  # add them
                  prev_timepoints <- allTimepoints[allTimepoints < timepoint]
                  aft_timepoints <- allTimepoints[allTimepoints > timepoint]
                  if(length(prev_timepoints)>0 & length(aft_timepoints)>0){
                    writeLines("LOG: augment variant detection with previous and following time points")
                    if(abs(max(prev_timepoints) - timepoint)*365 < opt$smoothingtime){
                      prev_specifiedLineages <- detect_lineages(DT_ = ssdt, timepoint_ = max(prev_timepoints), minInfoRatio_ = minInfoRatio)
                      specifiedLineages <- unique(c(specifiedLineages, prev_specifiedLineages))
                    }
                    if(abs(min(aft_timepoints) - timepoint)*365 < opt$smoothingtime){
                      aft_specifiedLineages <- detect_lineages(DT_ = ssdt, timepoint_ = min(prev_timepoints), minInfoRatio_ = minInfoRatio)
                      specifiedLineages <- unique(c(specifiedLineages, aft_specifiedLineages))
                    }
                  }
                  if(timepoint == max(timePoints) & length(prev_timepoints)>=1){
                    writeLines("LOG: augment variant detection be the two previous time points")

                    prev_timepoint <- sort(prev_timepoints, decreasing=TRUE)[1]
                    if(abs(prev_timepoint - timepoint)*365 < opt$smoothingtime){
                      prev_specifiedLineages <- detect_lineages(DT_ = ssdt, timepoint_ = prev_timepoint, minInfoRatio_ = minInfoRatio)
                      specifiedLineages <- unique(c(specifiedLineages, prev_specifiedLineages))
                    }
                  }
                  if(timepoint == max(timePoints) & length(prev_timepoints)>=2){
                    prevprev_timepoint <- sort(prev_timepoints, decreasing=TRUE)[2]
                    if(abs(prevprev_timepoint - timepoint)*365 < opt$smoothingtime){
                      prevprev_specifiedLineages <- detect_lineages(DT_ = ssdt, timepoint_ = prevprev_timepoint)
                      specifiedLineages <- unique(c(specifiedLineages, prevprev_specifiedLineages))
                    }
                  }
                }

                # use neigboring timepoints merged to deduce resulting list
                if(1){
                  # detect lineage in merged data from neighboring timepoints + current timepoint
                  # add lineages detected in any to specifiedLineages
                  prev_timepoints <- allTimepoints[allTimepoints < timepoint]
                  aft_timepoints <- allTimepoints[allTimepoints > timepoint]

                  adaptivMinInfoRatio <- 1-((1-minInfoRatio)^(1.1^length(allTimepoints)))
                  merged_specifiedLineages <- detect_lineages(DT_ = ssdt, timepoint_ = c(prev_timepoints, timepoint, aft_timepoints), minInfoRatio_ = adaptivMinInfoRatio)
                  specifiedLineages <- unique(c(specifiedLineages, merged_specifiedLineages))
                }



                writeLines(paste("LOG: (", t, ")", timepoint_classic, roiname, paste("(",length(allTimepoints[allTimepoints %in% timepoint]), "same day;", "+", length(allTimepoints[allTimepoints %notin% timepoint]), "neighboring TP;", length(specifiedLineages), "detected lineages)")))
                writeLines(paste0("LOG: detected lineages (", length(specifiedLineages), "): ", paste(specifiedLineages, sep = ", ", collapse = ", ")))

                if( length(specifiedLineages) > 0){
                  writeLines(paste("LOG:", "PERFORM REGRESSION WITH", length(specifiedLineages), "LINEAGES"))

                  # join model matrix columns to measured variables
                  smmat <- as.data.frame(mmat)[,which(colnames(mmat) %in% c("NUC", specifiedLineages))]
                  if (length(specifiedLineages) > 1){
                    smmat$groupflag <- apply( as.data.frame(mmat)[,which(colnames(mmat) %in% c(specifiedLineages))], 1 , paste , collapse = "" )
                  } else{
                    smmat$groupflag <-  as.data.frame(mmat)[,which(colnames(mmat) %in% c(specifiedLineages))]
                  }
                  left_join(x = ssdt, y = smmat, by = c("NUC"), multiple = "all") -> ssdt
                  ssdt %>% ungroup() %>% mutate(groupflag = paste(groupflag, ifelse(grepl(";", Variants), "M", "U"), sample_date_decimal)) -> ssdt

                  ## remove mutations which are not marker of any of the specifiedLineages
                  ## keep mutations which are marker of all of the specifiedLineages
                  if( sum(colnames(ssdt) %in% specifiedLineages) > 1){
                    ssdt[( rowSums(as.data.frame(ssdt)[,colnames(ssdt) %in% specifiedLineages ]) > 0 & rowSums(as.data.frame(ssdt)[,colnames(ssdt) %in% specifiedLineages ]) <= length(specifiedLineages)),] -> ssdt
                    ssdt[( rowSums(as.data.frame(ssdt)[,colnames(ssdt) %in% unique(specifiedLineages, highlightedVariants) ]) > 0 & rowSums(as.data.frame(ssdt)[,colnames(ssdt) %in% unique(specifiedLineages, highlightedVariants) ]) <= length(unique(specifiedLineages, highlightedVariants))),] -> detour4log
                  } else{
                    ssdt[as.data.frame(ssdt)[,colnames(ssdt) %in% specifiedLineages ] > 0,] -> ssdt
                    ssdt[as.data.frame(ssdt)[,colnames(ssdt) %in% unique(specifiedLineages, highlightedVariants) ] > 0,] -> detour4log
                  }
                  ### remove mutations which are marker of all of the specifiedLineages
                  #ssdt <- ssdt[rowSums(as.data.frame(ssdt)[colnames(ssdt) %in% specifiedLineages]) < length(specifiedLineages)]

                  ## remove zeros and low depth
                  ssdt %>% filter(value.freq > 0)  %>% filter(value.depth > min.depth) -> ssdt

                  ## transform to avoid 1
                  ssdt %>% mutate(value.freq = (value.freq * (value.depth-1) + 0.5)/value.depth) -> ssdt

                  # remove "OTHERS"
                  ssdt %>% mutate(Variants = gsub("other;*", "", Variants)) -> ssdt

                  ## remove outlier per group flag
                  ## remove if outside 1.5*IQR+/-0.05
                  dim(ssdt)[1] -> mutationsBeforeOutlierRemoval
                  ssdt %>% group_by(groupflag) %>% mutate(iqr = IQR(value.freq)) %>% mutate(upperbond = quantile(value.freq, 0.75) + 1.5 * iqr + 0.05, lowerbond = quantile(value.freq, 0.25) - 1.5 * iqr - 0.05) %>% filter((value.freq <= upperbond ) & (value.freq >= lowerbond)) %>% ungroup() %>% dplyr::select(-"groupflag", -"iqr", -"upperbond", -"lowerbond") -> ssdt
                  dim(ssdt)[1] -> mutationsAfterOutlierRemoval
                  if(mutationsAfterOutlierRemoval < mutationsBeforeOutlierRemoval){
                    writeLines(paste("LOG: ", mutationsBeforeOutlierRemoval-mutationsAfterOutlierRemoval, "mutations ignored since classified as outlier", "(", mutationsAfterOutlierRemoval, " remaining from", mutationsBeforeOutlierRemoval, ")"))
                  }

                  ## add unobserved uniq markers for specifiedLineages (after ignoring nonSpecifiedLineages) with AF 0
                  if(opt$addUniqZeros){
                    missedUniqMutations <- moi %>%
                            mutate(Variants = strsplit(as.character(Variants), ";")) %>%
                            unnest(Variants) %>%
                            filter(Variants %in% specifiedLineages) %>%
                            group_by(NUC) %>%
                            summarize(n = n(), Variants = paste(Variants, sep =";", collapse =";"), .groups = "drop") %>%
                            filter(n == 1) %>%
                            filter(NUC %notin% ssdt$NUC)
                    if(dim(missedUniqMutations)[1] > 0){

                      depthweight <- 10*min.depth
                      writeLines(paste("LOG: add", length(unique(missedUniqMutations$NUC)), " mutation to AF =", 0.5/depthweight, " since unobserved marker for", length(unique(missedUniqMutations$Variants)), "variants"))

                      for (j in 1:dim(missedUniqMutations)[1]){

                        mockssdt <- c()
                        mockssdt  <- c(mockssdt, as.character(missedUniqMutations[j,3]))
                        mockssdt  <- c(mockssdt, as.character(sdt %>% filter(sample_date %in% c(ref_timepoint_classic) ) %>% dplyr::select(ID) %>% distinct() %>% pull(ID) %>% paste(collapse = ",", sep = ",")))
                        mockssdt  <- c(mockssdt, as.character(ref_timepoint_classic))
                        mockssdt  <- c(mockssdt, as.character(ref_timepoint))
                        mockssdt  <- c(mockssdt, as.character(0.5/depthweight))
                        mockssdt  <- c(mockssdt, as.character(depthweight))
                        mockssdt  <- c(mockssdt, as.character(roi))
                        mockssdt  <- c(mockssdt, as.character(roiname))
                        mockssdt  <- c(mockssdt, as.character(missedUniqMutations[j,1]))
                        mockssdt <- as.data.table(matrix(mockssdt, nrow = 1))
                        colnames(mockssdt) <- colnames(ssdt)[1:9]
                        mockssdt$sample_date <- as.IDate(mockssdt$sample_date)
                        mockssdt$sample_date_decimal <- as.numeric(mockssdt$sample_date_decimal)
                        mockssdt$value.freq <- as.numeric(mockssdt$value.freq)
                        mockssdt$value.depth <- as.numeric(mockssdt$value.depth)


                        mockmmat <- c()
                        for (iij in colnames(ssdt)[colnames(ssdt) %in% specifiedLineages]){
                            Mm <- ifelse(iij == missedUniqMutations[j,3], 1, 0)
                            mockmmat  <- c(mockmmat, as.character(Mm))
                        }
                        mockmmat <- as.data.table(matrix(mockmmat, nrow = 1))
                        colnames(mockmmat) <- colnames(ssdt)[colnames(ssdt) %in% specifiedLineages]
                        mockmmat <- as.data.table(lapply(mockmmat, as.numeric))

                        mock <- cbind(mockssdt, mockmmat)
                        if(any(colnames(ssdt) == "groupflag")){
                          mock$groupflag <- "mockMissedUniq"
                        }
                        ssdt <- rbind(ssdt, mock)
                      }
                    }
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
                  writeLines(paste("  WARNING: since no mutations found in timepoint of interest (", timepoint, ") regression will be skipped"))
                    t <- t + 1
                    next;
                  }

                  ## add weight (1/(dayDiff+1))  *  log10(sequencingdepth)  *  (1/number_of_variants)^2  *  1/mutationLength
                  ssdt %>% rowwise() %>% mutate(varweight = (1/(sum(unlist(str_split(Variants, pattern = ";")) %in% specifiedLineages)))^2 ) %>%  mutate(timeweight = 1/(abs(timePoints[t] - sample_date_decimal)*(leapYear(floor(timePoints[t]))) + 1)) %>% mutate(mutationweight = 1/nchar(gsub("\\d", "", NUC))/2) %>% mutate(weight = 1*timeweight*varweight*mutationweight) -> ssdt
                  detour4log %>% rowwise() %>% mutate(varweight = (1/(sum(unlist(str_split(Variants, pattern = ";")) %in% specifiedLineages)))^2 ) %>%  mutate(timeweight = 1/(abs(timePoints[t] - sample_date_decimal)*(leapYear(floor(timePoints[t]))) + 1)) %>% mutate(weight = timeweight*varweight) %>% filter(!grepl(";.+;", Variants)) -> detour4log
                  detour4log %>% rowwise() %>% mutate(varweight = (1/(sum(unlist(str_split(Variants, pattern = ";")) %in% specifiedLineages)))^2 ) %>%  mutate(timeweight = 1/(abs(timePoints[t] - sample_date_decimal)*(leapYear(floor(timePoints[t]))) + 1)) %>% mutate(weight = timeweight*varweight) %>% dplyr::select(sample_date, NUC, value.freq, weight, all_of(specifiedLineages)) -> detour4log_printer
                  detour4log_printer %>% rowwise() %>% mutate(value.freq = signif(value.freq, digits = 2), weight = signif(weight, digits = 2)) -> detour4log_printer
                  detour4log %>% arrange(Variants) -> detour4log
                  print(data.frame(date = detour4log$sample_date, Tweight = round(detour4log$timeweight, digits = 2), depth = detour4log$value.depth, weight = round(detour4log$weight, digits = 2), value = round(detour4log$value.freq, digits = 2), variants = detour4log$Variants, mutation = detour4log$NUC))
                  rm(detour4log, detour4log_printer)

                  ## make regression
                  if(0){
                    method = "SIMPLEX"
                    fit1 <- tryCatch(gamlss(formula, data = ssdt, family = SIMPLEX, trace = FALSE, weights = weight),error=function(e) e, warning=function(w) w)
                    if(any(grepl("warning|error", class(fit1)))){
                      method = "BE"
                      fit1 <- tryCatch(gamlss(formula, data = ssdt, family = BE, trace = FALSE, weights = weight),error=function(e) e, warning=function(w) w)
                      writeLines(paste("LOG: fall back to BE"))
                      if(any(grepl("warning|error", class(fit1)))){
                        writeLines(paste("  WARNING: BE did not converge at", timePoints[t], " .. skipped"))
                        t <- t + 1
                        next; ## remove if unconverged BE results should be used
                      }
                    }
                  }
                  if(1){
                    method = "BE"
                    fit1 <- tryCatch(gamlss(formula, data = ssdt, family = BE, trace = FALSE, weights = weight),error=function(e) e, warning=function(w) w)
                    writeLines(paste("LOG: fall back to BE"))
                    if(any(grepl("warning|error", class(fit1)))){
                      writeLines(paste("LOG: BE did not converge either at", timePoints[t], " .. skipped"))
                      t <- t + 1
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
                  writeLines(paste("WARNING:", "no sample ID for", roi, "@", timePoints[t]))
                  }
                  ssdtFit$ID = rep(unique(sample_ID$ID)[length(unique(sample_ID$ID))], n = length(ssdtFit$Variants))

                  if(unique(ssdtFit$T)>1){
                    writeLines(paste("  LOG: fitted value corrected for >1: T =", unique(ssdtFit$T), "; method used =", method))
                  }

                  ssdt %>% ungroup() %>% mutate(all_Variants = Variants) %>% mutate(Variants = strsplit(as.character(Variants), ";")) %>% unnest(Variants) %>% filter(Variants %in% specifiedLineages) %>% mutate(T = unique(ssdtFit$T)) %>% mutate(fit2 = ifelse(T > 1, fit1/T, fit1)) %>% filter(sample_date_decimal == timePoints[t]) -> ssdt2
                  plantFullData <- rbind(plantFullData, data.table(variant = ssdt2$Variants, all_variants = ssdt2$all_Variants, LocationID =  rep(roi, length(ssdt2$fit1)), LocationName  =  rep(roiname, length(ssdt2$fit1)), sample_date = as.character(ssdt2$sample_date), value = ssdt2$fit2, marker = ssdt2$NUC, singlevalue = ssdt2$value.freq ))

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
                  writeLines(paste("LOG:", "NOTHIN FOR REGRESSION AT", timePoints[t], "since number of lineages detected ==", length(specifiedLineages), "; ... ignore TIMEPOINT"))
                }
            }, error=function(e) {print("WARNING: Iteration failed somewhere")})
            if(is.null(checkSuccess)){
              t <- t + 1
              catchCounter = 0
              rm(checkSuccess)
              next
            } else if(catchCounter < 3){
              writeLines(paste("LOG: repeat loop for timepoint", t, "for", catchCounter, "times"))
              catchCounter <- catchCounter + 1
              rm(checkSuccess)
              next
            } else{
              writeLines(paste("LOG: stop repeating loop for", t))
              t <- t + 1
              catchCounter = 0
              rm(checkSuccess)
              next
            }
          }
    }
    writeLines(paste0("PROGRESS: finished variant quantification for all timePoints "))

    # average per sample_date to take care of more than one sample a day
    plantFittedData <- plantFittedData %>% group_by(variant, LocationID, LocationName, sample_date) %>% summarize(value = mean(value), .groups = "keep")
    sampleID2sample_date <- metaDT %>% filter(LocationID == roi) %>% group_by(sample_date) %>% summarize(BSF_sample_name = paste(BSF_sample_name, collapse = ","), .groups = "keep")
    sampleID2sample_date$sample_date <- as.character(sampleID2sample_date$sample_date)
    plantFittedData <- left_join(x = plantFittedData, y = sampleID2sample_date, by = "sample_date")
    plantFittedData <- plantFittedData %>% dplyr::select(variant, LocationID, LocationName, sample_id=BSF_sample_name, sample_date, value)

    # remove variants which are never observed
    if(dim(plantFittedData)[1] == 0){
      next;
    }
    plantFittedData %>% group_by(variant) %>% mutate(maxValue = max(value)) %>% filter(maxValue>0) %>% dplyr::select(-"maxValue") -> plantFittedData


    writeLines(paste("PROGRESS: start plotting ", roiname))

    if(dim(plantFittedData)[1] > 0){
        plantFittedData %>% mutate(latest = max(sample_date, na.rm = TRUE)) %>% filter(sample_date == latest) %>% filter(value > 0) %>% summarize(variant = variant, freq = value, .groups = "keep") %>% distinct() -> sankey.dt

        plantFittedData  %>% group_by() %>% summarize(latest = max(sample_date, na.rm = TRUE), .groups = "keep") -> sankey_date

        sankey.dt  %>% rowwise() %>% mutate(freq = ifelse(freq<0, 0, freq)) -> sankey.dt
        sankey.dt  %>% rowwise() %>% mutate(freq = ifelse(freq>1, 1, freq)) -> sankey.dt
        if(sum(sankey.dt$freq)>1){
          sankey.dt  %>% ungroup() %>% mutate(freq = freq/sum(freq)) -> sankey.dt
        }
        sankey.dt %>% group_by(variant = "0") %>% summarize(freq = 1-sum(freq), .groups = "keep") -> sankey.dt0
        rbind(sankey.dt, sankey.dt0) -> sankey.dt
        sankey.dt %>% rowwise() %>% filter(freq > zeroo/10) -> sankey.dt

        sankey.dt$variant_dealias <- unlist(lapply(as.list(sankey.dt$variant), dealias))
        for ( baseForColor in baseForColorsVariants ){
          if(any(grepl(baseForColor, sankey.dt$variant_dealias, fixed = TRUE))){
              if(baseForColor %notin% sankey.dt$variant_dealias){
                  rbind(sankey.dt, data.frame(variant = realias(baseForColor), freq = 0, variant_dealias = baseForColor)) -> sankey.dt
              }
          }
        }

        sankey.dt %>% rowwise() %>% mutate(levels = length(unlist(strsplit(variant_dealias, split = "\\.")))) %>% group_by(1) %>% summarize( maxLev = max(levels), .groups = "keep") -> maxLev
        sankey.dt %>% rowwise() %>% mutate(long_variant = expandVar(variant_dealias, maxLev$maxLev)) -> sankey.dt
        ancestor <- getTree(sankey.dt)
        Nlevels <- max(unlist(lapply(ancestor, length)))

        sankey.dt <- makeObs(sankey.dt, Nlevels, ancestor, sankeyPrecision)

        sankey.dt %>% make_long(names(sankey.dt)) -> sankey.dtt

        sankey.dtt %>% filter(!is.na(node)) -> sankey.dtt
        max(as.numeric(gsub("level", "", sankey.dtt$next_x)), na.rm = TRUE) -> maxLev

        #sankey.dtt %>% group_by(x, node) %>% mutate(n = paste0("[", 100*n()/sankeyPrecision, "%]")) %>% rowwise() %>% mutate(level = as.numeric(gsub("level", "", x))) %>% mutate(label = ifelse(level == maxLev, paste(node, n), node)) -> sankey.dtt
        sankey.dtt %>% group_by(x, node) %>% mutate(n = paste0("[", 100*n()/sankeyPrecision, "%]")) %>% rowwise() %>% mutate(level = as.numeric(gsub("level", "", x))) %>% mutate(label = paste(node, n)) -> sankey.dtt

        sankey.dtt %>% rowwise() %>% mutate(node = ifelse(is.na(node), node, dealias(node))) %>% mutate(next_node = ifelse(is.na(next_node), next_node, dealias(next_node))) -> sankey.dtt

        sankey.dtt %>% rowwise() %>% mutate(label = ifelse((node == "0" & level>0), gsub("^0", "unclassified", label), label)) -> sankey.dtt
        sankey.dtt <- sankey.dtt %>% mutate(label = ifelse((grepl("unclassified", label) & level != maxLev), NA, label))
        sankey.dtt %>% group_by(node, x) %>% mutate(label = ifelse(length(unique(next_node)) == 1 & node == next_node & level != maxLev, NA, label)) -> sankey.dtt

        ggplot(sankey.dtt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = label)) +
            geom_sankey(flow.alpha = .6, type ='alluvial') +
            geom_sankey_label(size = 2, color = "white", fill = "gray40", position = position_nudge(x = 0.05, y = 0), na.rm = TRUE, type ='alluvial', hjust = 0) +
            theme_sankey(base_size = 16) +
            labs(x = NULL) +
            ggtitle(roiname, subtitle = sankey_date) +
            theme(legend.position = "none", axis.text.x = element_blank(), plot.title = element_text(hjust = 0), plot.subtitle=element_text(hjust = 0)) +
            scale_fill_viridis_d(alpha = 1, begin = 0.025, end = .975, direction = 1, option = "D") +
            scale_x_discrete(expand = expansion(mult = c(0, .1), add = c(.1, .8))) -> pp

        filename <- paste0(outdir, "/figs/sankey/",  paste('/wwtp', roi, sep="_"), ".pdf")
        ggsave(filename = filename, plot = pp, width = 6, height = 4.4)
        fwrite(as.list(c("sankey", "WWTP", roiname, filename)), file = summaryDataFile, append = TRUE, sep = "\t")
        rm(pp, filename, sankey.dt, sankey.dtt)
    }

    if(length(unique(plantFittedData$sample_date)) >= tpLimitToPlot){

      if(dim(plantFittedData)[1] >= 1){
        writeLines(paste("PROGRESS: plotting stackedview", roiname))

        ## print stacked area overview of all detected lineages
        plantFittedData2 <- plantFittedData

        ## sort variants according dealiased names, priotizing highlightedVariants
        plantFittedData2 %>%
          dplyr::select(variant) %>%
          distinct() %>%
          rowwise() %>%
          mutate(variant_dealias = dealias(variant)) %>%
          arrange(variant_dealias) %>%
          pull(variant) -> variant_order

        plantFittedData2$variant <- factor(plantFittedData2$variant, levels = variant_order)

        ## avoid >1 sums per time point
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


        plottng_data$variant_dealiased <- unlist(lapply(as.list(as.character(plottng_data$variant)), dealias))
        stackorder <- plottng_data %>% ungroup() %>%  dplyr::select(variant, variant_dealiased) %>% distinct() %>% arrange(variant_dealiased) %>% pull(variant)
        stackorder <- as.character(stackorder)

        stackorder_dealias <- unlist(lapply(as.list(stackorder), dealias))
        highlightedVariants_dealias <- unlist(lapply(as.list(highlightedVariants), dealias))
        stackorder_highlight <- c()
        for (v in seq_along(highlightedVariants)){
          vv <- highlightedVariants[v]
          vd <- highlightedVariants_dealias[v]
          nv <- grep(vd, stackorder_dealias)
          stackorder_highlight <- c(stackorder_highlight, nv[nv %notin% stackorder_highlight])
        }
        stackorder_highlight <- c(stackorder_highlight, seq_along(stackorder)[seq_along(stackorder) %notin% stackorder_highlight])
        stackorder_highlight <- stackorder[stackorder_highlight]


        plottng_data$variant <- factor(plottng_data$variant, levels = rev(stackorder_highlight))

        plottng_data %>% filter(sample_date == max(sample_date, na.rm = TRUE)) %>% filter(value > 0) -> plottng_labels
        plottng_labels$variant <- factor(plottng_labels$variant, levels = rev(stackorder_highlight[stackorder_highlight %in% plottng_labels$variant]))

        ggplot(data = plottng_data, aes(x = as.Date(sample_date), y = value, fill = variant, color = variant)) -> q3
        q3 <- q3 + geom_area(position = "stack", alpha = 0.6)
        q3 <- q3 + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5)
        q3 <- q3 + scale_fill_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "")
        q3 <- q3 + scale_color_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "")
	      q3 <- q3 + scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks = c(0,0.5,1))
        q3 <- q3 + theme_minimal()
        q3 <- q3 + theme(legend.position="none", strip.text.x = element_text(size = 4.5), panel.grid.minor = element_blank(), panel.spacing.y = unit(0, "lines"), legend.direction="horizontal")
        q3 <- q3 + guides(fill = guide_legend(title = "", ncol = 7, override.aes = aes(label = "")), color = guide_legend(title = "", ncol = 7, override.aes = aes(label = "")))
        q3 <- q3 + scale_x_date(date_breaks = "2 month", date_labels =  "%b %y", limits = c(as.Date(NA), as.Date(latestSample$latest)+14))
        q3 <- q3 + ylab(paste0("Variantenanteil [1/1]") )
        q3 <- q3 + xlab("Kalender Woche")
        q3 <- q3 + geom_text_repel(data=plottng_labels, aes(x=as.Date(sample_date), y=value, label=variant), position = position_stack(), min.segment.length = .01, direction = "y", alpha = 0.6, xlim = c(max(as.Date(plottng_labels$sample_date)), Inf), size = 2.5)

        ggplot(data = plottng_data, aes(x = as.Date(sample_date), y = value, fill = variant, color = variant)) -> q3s
        q3s <- q3s + geom_area(alpha = 0.8)
        q3s <- q3s + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5)
        q3s <- q3s + facet_wrap(~variant, ncol = 5, scales = "free_y")
        q3s <- q3s + scale_fill_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "")
        q3s <- q3s + scale_color_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "")
        q3s <- q3s + scale_y_continuous(labels=scales::percent, n.breaks = 3)
        q3s <- q3s + theme_minimal()
        q3s <- q3s + theme(legend.position="none", axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 6), strip.text.x = element_text(size = 5), panel.grid.minor = element_blank(), panel.spacing.y = unit(0, "lines"), legend.direction="horizontal")
        q3s <- q3s + guides(fill = guide_legend(title = "", ncol = 7, override.aes = aes(label = "")), color = guide_legend(title = "", ncol = 7, override.aes = aes(label = "")))
        q3s <- q3s + scale_x_date(date_breaks = "4 month", date_labels =  "%b %y", limits = c(as.Date(NA), as.Date(latestSample$latest)+14))
        q3s <- q3s + ylab(paste0("Variantenanteil [1/1]") )
        q3s <- q3s + xlab("")


        filename <- paste0(outdir, '/figs/stackview', paste('/wwtp', roi, "all", sep="_"), ".pdf")
        ggsave(filename = filename, plot = plot_grid(q3, q3s, ncol = 1, rel_heights = c(2,3)), width = plotWidth, height = 1.4*plotHeight)
        fwrite(as.list(c("stackOverview", c( ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "all", filename))), file = summaryDataFile, append = TRUE, sep = "\t")

        rm(q3, filename, plantFittedData2, FillIdx, LabelTxt, GrpIdx, g, plottng_labels)

        ## print faceted line plot of fitted values for all detected lineages
        plantFittedData %>% filter(sample_date > as.Date(latestSample$latest)-(2*recentEnought)) %>% group_by(variant) %>% mutate(maxValue = pmax(value)) %>% filter(maxValue > 0) %>% ggplot(aes(x = as.Date(sample_date), y = value, col = variant)) + geom_line(alpha = 0.6) + geom_point(alpha = 0.66, shape = 13) + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5) + scale_color_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "Varianten") + scale_y_continuous(labels=scales::percent, trans = "log10") + theme_bw() + theme(legend.position="bottom", strip.background = element_rect(fill="grey97")) + scale_x_date(date_breaks = "2 weeks", date_labels =  "%b %d", limits = c(as.Date(latestSample$latest)-(2*recentEnought), as.Date(latestSample$latest))) + ylab(paste0("Variantenanteil [1/1]") ) + xlab("") -> q2

        spemut_sdt %>% filter(value.freq > 0) %>% dplyr::select("ID", "sample_date", "value.freq", "LocationID", "LocationName", "NUC") %>% filter(sample_date > as.Date(latestSample$latest)-(2*recentEnought)) -> spemut_draw2
        if(dim(spemut_draw2)[1] > 0){
          writeLines(paste("PROGRESS: plotting special mutations", roiname))
          left_join(x = spemut_draw2, y = soi, by = "NUC", multiple = "all") -> spemut_draw2
          spemut_draw2 %>% rowwise() %>% mutate(marker =  paste(sep = "|", NUC, paste(Gene, AA, sep = ":"), paste0("[", Variants, "]")) ) -> spemut_draw2
          colnames(spemut_draw2)[colnames(spemut_draw2) == "Variants"] <- "variant"
          q2 + geom_point(data = spemut_draw2, aes(x = as.Date(sample_date), y = value.freq, shape = marker, fill = marker), color = "black", size = 2, alpha = .45) -> q2
          q2 + scale_shape_manual(values = c(0:25, 33:50)) -> q2
          q2 + guides(shape = guide_legend(title = "Spezial-Mutationen", ncol = 2, title.position = "top"), fill = guide_legend(title = "Spezial-Mutationen", ncol = 2, title.position = "top"), col = guide_legend(title = "Varianten", ncol = 3, title.position = "top")) -> q2

          filename <- paste0(outdir, '/figs/specialMutations', paste('/wwtp', roi, "all", sep="_"), ".pdf")
          #ggsave(filename = filename, plot = q2, width = plotWidth/1.5, height = plotHeight/1.5)
          #fwrite(as.list(c("specialMutations", c( ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "all",filename))), file = summaryDataFile, append = TRUE, sep = "\t")
        }
        rm(q2,filename,spemut_draw2)

      }

      if(dim(filter(plantFittedData, variant %in% VoI))[1] >= 1){

        ## print faceted line plot of fitted values for all VoI
        plantFittedData %>% filter(variant %in% VoI) %>% filter(sample_date > as.Date(latestSample$latest)-(2*recentEnought)) %>% group_by(variant) %>% mutate(maxValue = pmax(value)) %>% filter(maxValue > 0) %>% ggplot(aes(x = as.Date(sample_date), y = value, col = variant)) + geom_line(alpha = 0.6) + geom_point(alpha = 0.66, shape = 13) + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5)  + scale_color_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "Varianten") + scale_y_continuous(labels=scales::percent, trans = "log10") + theme_bw() + theme(legend.position="bottom", strip.background = element_rect(fill="grey97")) + scale_x_date(date_breaks = "2 weeks", date_labels =  "%b %d", limits = c(as.Date(latestSample$latest)-(2*recentEnought), as.Date(latestSample$latest))) + ylab(paste0("Variantenanteil [1/1]") ) + xlab("") -> q1

        spemut_sdt %>% filter(value.freq > 0) %>% dplyr::select("ID", "sample_date", "value.freq", "LocationID", "LocationName", "NUC")  %>% filter(sample_date > as.Date(latestSample$latest)-(2*recentEnought)) -> spemut_draw1
        if(dim(spemut_draw1)[1] > 0){
          writeLines(paste("PROGRESS: plotting special mutations VoI", roiname))
          left_join(x = spemut_draw1, y = soi, by = "NUC", multiple = "all") -> spemut_draw1
          spemut_draw1 %>% rowwise() %>% mutate(marker = paste(sep = "|", NUC, paste(Gene, AA, sep = ":"), paste0("[", Variants, "]")) ) -> spemut_draw1
          colnames(spemut_draw1)[colnames(spemut_draw1) == "Variants"] <- "variant"
          q1 + geom_point(data = spemut_draw1, aes(x = as.Date(sample_date), y = value.freq, shape = marker, fill = marker), color = "black", size = 2, alpha = .45) -> q1
          q1 + scale_shape_manual(values = c(0:25, 33:50)) -> q1
          q1 + guides(shape = guide_legend(title = "Spezial-Mutationen", ncol = 2, title.position = "top"), fill = guide_legend(title = "Spezial-Mutationen", ncol = 2, title.position = "top"), col = guide_legend(title = "Varianten", ncol = 3, title.position = "top")) -> q1


          filename <- paste0(outdir, '/figs/specialMutations', paste('/wwtp', roi, "VoI", sep="_"), ".pdf")
          ggsave(filename = filename, plot = q1, width = plotWidth/1.5, height = plotHeight/1.5)
          fwrite(as.list(c("specialMutations", c( ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "VoI", filename))), file = summaryDataFile, append = TRUE, sep = "\t")
        }
        rm(q1, filename, spemut_draw1)
      }

      if(dim(filter(plantFullData, all_variants %in% VoI))[1] >= 1){
        writeLines(paste("PROGRESS: plotting variantDetails VoI", roiname))
        ## print faceted line plot of fitted values plus point plot of measured AF for all VoIs
        plantFullData %>% filter(all_variants %in% VoI) %>%
              ggplot() +
              geom_line(data = subset(plantFittedData, variant %in% VoI), alpha = 0.6, size = 2, aes(x = as.Date(sample_date), y = value, col = variant)) +
              geom_jitter(aes(x = as.Date(sample_date), y = singlevalue), alpha = 0.33, size = 1.5, width = 0.33, height = 0) +
              geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5) +
              facet_wrap(~variant)  +
              scale_color_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "") +
              scale_y_continuous(labels=scales::percent) +
              theme_bw() +
              theme(legend.position="bottom", strip.background = element_rect(fill="grey97"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
              scale_x_date(date_breaks = "2 month", date_labels =  "%b %y", limits = c(as.Date(NA), as.Date(latestSample$latest))) +
              ylab(paste0("Variantenanteil [1/1]") ) +
              xlab("") -> p1

        filename <- paste0(outdir, '/figs/detail', paste('/wwtp', roi, "VoI", sep="_"), ".pdf")
        ggsave(filename = filename, plot = p1, width = plotWidth, height = plotHeight)
        fwrite(as.list(c("variantDetail", c( ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "VoI", filename))), file = summaryDataFile, append = TRUE, sep = "\t")

        rm(p1, filename)
      }

      if(dim(plantFullData)[1] >= 1){
        writeLines(paste("PROGRESS: plotting variant Details", roiname))
        ## print faceted line plot of fitted values plus point plot of measured AF for all lineages
        plantFullData %>% filter(!grepl(";", all_variants)) %>%
              ggplot() +
              geom_line(data = plantFittedData, alpha = 0.6, size = 2, aes(x = as.Date(sample_date), y = value, col = variant)) +
              geom_jitter(aes(x = as.Date(sample_date), y = singlevalue), alpha = 0.33, size = 1.5, width = 0.33, height = 0) +
              geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5) +
              facet_wrap(~variant) +
              scale_color_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "") +
              scale_y_continuous(labels=scales::percent) +
              theme_bw() + theme(legend.position="bottom", strip.background = element_rect(fill="grey97"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
              scale_x_date(date_breaks = "3 month", date_labels =  "%b", limits = c(as.Date(NA), as.Date(latestSample$latest))) +
              ylab(paste0("Variantenanteil [1/1]") ) +
              xlab("") +
              guides(color = guide_legend(title = "", ncol = 7)) -> p2

        filename <- paste0(outdir, '/figs/detail', paste('/wwtp', roi, "all", sep="_"), ".pdf")
        ggsave(filename = filename, plot = p2, width = plotWidth, height = plotHeight*1.6)
        fwrite(as.list(c("variantDetail", c(ifelse(dim(count_last_BSF_run_id)[1] > 0, "current", "old"), roiname, "all", filename))), file = summaryDataFile, append = TRUE, sep = "\t")
        rm(p2, filename)
      }
    }

    globalFittedDataList[[length(globalFittedDataList)+1]] <- plantFittedData
    globalFullDataList[[length(globalFullDataList)+1]] <- plantFullData


    rm(plantFittedData, plantFullData)
}

globalFittedData <- rbindlist(globalFittedDataList)
globalFullData <- rbindlist(globalFullDataList)

## complete set to add zeros for unobserved variants
reshape2::dcast(globalFittedData, sample_id~variant, value.var = "value") -> globalFittedData_completed
globalFittedData_completed[is.na(globalFittedData_completed)] <- 0
globalFittedData_completed <- reshape2::melt(globalFittedData_completed, id.vars="sample_id", variable.name = "variant", value.name = "value")
globalFittedData <- right_join(x = (globalFittedData %>% dplyr::select(LocationID, LocationName, sample_id, sample_date) %>% distinct()), y = globalFittedData_completed, by = c("sample_id"))
globalFittedData <- globalFittedData %>% dplyr::select(variant, LocationID, LocationName, sample_id, sample_date, value)


writeLines(paste("PROGRESS: loop over WWTP finished"))
## dump global data into output file
writeLines(paste("PROGRESS: writing result tables"))
globalFittedData %>% distinct() -> globalFittedData
globalFullData %>% distinct() -> globalFullData
globalFullSoiData %>% distinct() -> globalFullSoiData


fwrite(globalFittedData , file = paste0(outdir, "/globalFittedData.csv"), sep = "\t")
fwrite(globalFullData, file = paste0(outdir, "/globalFullData.csv.gz"), sep = "\t")
fwrite(globalFullSoiData, file = paste0(outdir, "/globalSpecialmutData.csv.gz"), sep = "\t")

## print overview of all VoI and all plants
writeLines(paste("PROGRESS: start plotting overviews"))




##

if (dim(globalFittedData)[1] >= tpLimitToPlot){

  writeLines(paste("PROGRESS: plotting aggregated stack plot"))

  globalFittedData$sample_date <- as.Date(globalFittedData$sample_date)
  metaDT$sample_date <- as.Date(metaDT$sample_date)

  left_join(x = globalFittedData, y = (metaDT %>% dplyr::select(LocationID, connected_people)), by = "LocationID", multiple = "all") %>% distinct()  -> stacker.dt
  if(any(stacker.dt$connected_people > 0)){
    stacker.dt %>% rowwise() %>% mutate(kw = as.Date(sample_date, tryFormats = c("%Y-%m-%d")) + 6 - as.numeric(strftime(as.Date(sample_date, tryFormats = c("%Y-%m-%d")), format = "%u"))) %>% group_by(variant, kw) %>% summarize(agg_value = weighted.mean(value, connected_people), .groups = "keep") -> stacker.dt

    stacker.dt %>% group_by(variant) %>% mutate(max = max(agg_value)) %>% filter(max > 0) -> stacker.dt
    reshape2::melt(reshape2::dcast(stacker.dt, kw~variant, value.var = "agg_value", fill = 0), id.vars = c("kw"), variable.name = "variant", value.name = "agg_value") %>% filter(!is.na(kw)) -> stacker.dtt

    ColorBaseData <- data.table(variant = character(), base = character())
    for (varr in unique(stacker.dtt$variant)){
      foundbase = "other"
      varrr <- dealias(varr)
      for (base in baseForColorsVariants){
        if(grepl(base, varrr, fixed = TRUE)){
          foundbase = base
        }
      }
      rbind(data.table(variant = varr, base = foundbase), ColorBaseData) -> ColorBaseData
    }

    ## adjust if mean sum per week is >1
    stacker.dtt %>% group_by(kw) %>% mutate(T = sum(agg_value)) %>% rowwise() %>% mutate(agg_value = ifelse(T > 1, agg_value/(T+0.00001), agg_value)) %>% dplyr::select(-"T") -> stacker.dtt

    ColorBaseData$variant_dealiased <- unlist(lapply(as.list(as.character(ColorBaseData$variant)), dealias))
    ColorBaseData[order(ColorBaseData$variant_dealiased, decreasing = TRUE),] -> ColorBaseData

    ColorBaseData %>% group_by(base) %>%  mutate(n = n()) %>% arrange(base, variant_dealiased) %>% dplyr::mutate(id = cur_group_id()) %>% mutate(id = ifelse(id > 6, 6, id)) -> ColorBaseData
    ColorBaseData %>% group_by(id) %>% mutate(i = row_number()) %>% rowwise() %>% mutate(col = getColor(n, id, i)) -> ColorBaseData

    stacker.dtt$variant_dealiased <- unlist(lapply(as.list(as.character(stacker.dtt$variant)), dealias))
    stackorder <- stacker.dtt %>% ungroup() %>%  dplyr::select(variant, variant_dealiased) %>% distinct() %>% arrange(variant_dealiased)
    stacker.dtt$variant <- factor(stacker.dtt$variant, levels = rev(stackorder$variant))

    stack_order <- stacker.dtt %>% ungroup() %>% dplyr::select(variant) %>%
      distinct() %>%
      rowwise() %>%
      mutate(variant_dealias = dealias(variant)) %>%
      arrange(variant_dealias) %>%
      pull(variant)
    stack_order_dealias <- unlist(lapply(as.list(stack_order), dealias))
    highlightedVariants_dealias <- unlist(lapply(as.list(highlightedVariants), dealias))
    stack_order_highlight <- c()
    for (v in seq_along(highlightedVariants)){
      vv <- highlightedVariants[v]
      vd <- highlightedVariants_dealias[v]
      nv <- grep(vd, stack_order_dealias)
      stack_order_highlight <- c(stack_order_highlight, nv[nv %notin% stack_order_highlight])
    }
    stack_order_highlight <- c(stack_order_highlight, seq_along(stack_order)[seq_along(stack_order) %notin% stack_order_highlight])
    stack_order <- stack_order[stack_order_highlight]
    stacker.dtt$variant <- factor(stacker.dtt$variant, levels = rev(stack_order))

    stacker.labels <- stacker.dtt %>% ungroup() %>% filter(kw == max(kw)) %>% filter(agg_value > 0)
    stacker.labels$variant <- factor(stacker.labels$variant, levels = rev(stack_order[stack_order %in% stacker.labels$variant]))


    ap <- ggplot(data = stacker.dtt, aes(x = as.Date(kw), y = agg_value, fill = variant, color = variant))
    ap <- ap + geom_area(position = "stack", alpha = 0.6)
    ap <- ap + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5)
    ap <- ap + scale_fill_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "")
    ap <- ap + scale_color_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "")
    ap <- ap + scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks = c(0,0.5,1))
    ap <- ap + theme_minimal()
    ap <- ap + theme(legend.position="none", strip.text.x = element_text(size = 4.5), panel.grid.minor = element_blank(), panel.spacing.y = unit(0, "lines"), legend.direction="horizontal")
    ap <- ap + guides(fill = guide_legend(title = "", ncol = 7, override.aes = aes(label = "")), color = guide_legend(title = "", ncol = 7, override.aes = aes(label = "")))
    ap <- ap + scale_x_date(date_breaks = "2 month", date_labels =  "%b %y", limits = c(as.Date(NA), as.Date(latestSample$latest)+14))
    ap <- ap + ylab(paste0("Variantenanteil [1/1]") )
    ap <- ap + xlab("Kalender Woche")
    ap <- ap + ggtitle("Gewichtetes Mittel: Österreich")
    ap <- ap + geom_text_repel(data=stacker.labels, aes(x=as.Date(kw), y=agg_value, label=variant), position = position_stack(), min.segment.length = .01, direction = "y", alpha = 0.6, xlim = c(max(as.Date(stacker.labels$kw)), Inf), size = 2.5)

    aps <- ggplot(data = stacker.dtt, aes(x = as.Date(kw), y = agg_value, fill = variant, color = variant))
    aps <- aps + geom_area(alpha = 0.8)
    aps <- aps + geom_vline(xintercept = as.Date(latestSample$latest), linetype = "dotdash", color = "grey", size =  0.5)
    aps <- aps + facet_wrap(~variant, ncol = 5, scales = "free_y")
    aps <- aps + scale_fill_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "")
    aps <- aps + scale_color_manual(values = ColorBaseData$col, breaks = ColorBaseData$variant, name = "")
    aps <- aps + scale_y_continuous(labels=scales::percent, n.breaks = 3)
    aps <- aps + theme_minimal()
    aps <- aps + theme(legend.position="none", axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 6), strip.text.x = element_text(size = 5), panel.grid.minor = element_blank(), panel.spacing.y = unit(0, "lines"), legend.direction="horizontal")
    aps <- aps + guides(fill = guide_legend(title = "", ncol = 7, override.aes = aes(label = "")), color = guide_legend(title = "", ncol = 7, override.aes = aes(label = "")))
    aps <- aps + scale_x_date(date_breaks = "4 month", date_labels =  "%b %y", limits = c(as.Date(NA), as.Date(latestSample$latest)+14))
    aps <- aps + ylab(paste0("Variantenanteil [1/1]") )
    aps <- aps + xlab("")


    filename <- paste0(outdir, '/figs/stackview', '/', paste(opt$country, "all", sep="_"), ".pdf")
    ggsave(filename = filename, plot = plot_grid(ap, aps, ncol = 1, rel_heights = c(1,3)), width = plotWidth, height = 1.6*plotHeight)
    fwrite(as.list(c("stacked", "Overview", opt$country, filename)), file = summaryDataFile, append = TRUE, sep = "\t")

    rm(ap,filename, ColorBaseData, stacker.dt, stacker.dtt)
  }

}

if(dim(globalFittedData)[1] > 0){
    writeLines(paste("PROGRESS: plotting overview Sankey plot + Detection Plot"))

    if(as.Date(format(Sys.time(), "%Y-%m-%d")) - max(globalFittedData$sample_date, na.rm = TRUE) > recentEnought){
      recentEnought.original <- recentEnought
      recentEnought = 1+as.numeric(as.Date(format(Sys.time(), "%Y-%m-%d")) - max(globalFittedData$sample_date, na.rm = TRUE))
      writeLines(paste("WARNING: not data points classify as 'recentEnought' for overview Sankey plot; value increased from", recentEnought.original, "to", recentEnought))
    }

    ## make complete plot
    globalFittedData %>% filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date) %>%
        group_by(LocationID) %>%
        mutate(latest = max(sample_date, na.rm = TRUE)) %>%
        filter(sample_date == latest) %>%
        summarize(variant = variant, freq = value, .groups = "keep") -> sankey_all.dt


    if(dim(sankey_all.dt)[1] > 0){

        globalFittedData %>% filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date) %>%
            group_by(LocationID) %>%
            mutate(latest = max(sample_date, na.rm = TRUE)) %>%
            filter(sample_date == latest) %>%
            ungroup() %>%
            summarize(earliest = min(sample_date, na.rm = TRUE), latest = max(sample_date, na.rm = TRUE), .groups = "keep") -> sankey_date

        reshape2::dcast(sankey_all.dt, LocationID~variant, value.var = "freq") -> sankey_completed.dt
        sankey_completed.dt[is.na(sankey_completed.dt)] <- 0
        sankey_completed.dt <- reshape2::melt(sankey_completed.dt, id.vars="LocationID", variable.name = "variant", value.name = "freq")

        left_join(x = sankey_completed.dt, y = (metaDT %>% dplyr::select(LocationID, connected_people) %>%
            distinct()), by = "LocationID", multiple = "all") %>%
            group_by(variant) %>%
            summarize(freq = weighted.mean(freq, w = connected_people, na.rm = TRUE), .groups = "keep") -> sankey.dt

        if(dim(sankey.dt)[1] > 0){
          sankey.dt  %>% rowwise() %>% mutate(freq = ifelse(freq<0, 0, freq)) -> sankey.dt
          sankey.dt  %>% rowwise() %>% mutate(freq = ifelse(freq>1, 1, freq)) -> sankey.dt
          if(sum(sankey.dt$freq)>1){
            sankey.dt  %>% ungroup() %>% mutate(freq = freq/sum(freq)) -> sankey.dt
          }
          sankey.dt  %>% ungroup() -> sankey.dt
          sankey.dt -> occurence.freq

          sankey.dt %>% group_by(variant = "0") %>% summarize(freq = 1-sum(freq), .groups = "keep") -> sankey.dt0
          rbind(sankey.dt, sankey.dt0) -> sankey.dt
          sankey.dt %>% rowwise() %>% filter(freq > zeroo/10) -> sankey.dt

          sankey.dt$variant_dealias <- unlist(lapply(as.list(sankey.dt$variant), dealias))
          for ( baseForColor in baseForColorsVariants ){
            if(any(grepl(baseForColor, sankey.dt$variant_dealias, fixed = TRUE))){
                if(baseForColor %notin% sankey.dt$variant_dealias){
                    rbind(sankey.dt, data.frame(variant = realias(baseForColor), freq = 0, variant_dealias = baseForColor)) -> sankey.dt
                }
            }
          }

          sankey.dt %>% rowwise() %>% mutate(levels = length(unlist(strsplit(variant_dealias, split = "\\.")))) %>% group_by(1) %>% summarize( maxLev = max(levels), .groups = "keep") -> maxLev
          sankey.dt %>% rowwise() %>% mutate(long_variant = expandVar(variant_dealias, maxLev$maxLev)) -> sankey.dt
          ancestor <- getTree(sankey.dt)
          Nlevels <- max(unlist(lapply(ancestor, length)))

          sankey.dt <- makeObs(sankey.dt, Nlevels, ancestor, sankeyPrecision)

          sankey.dt %>% make_long(names(sankey.dt)) -> sankey.dtt

          sankey.dtt %>% filter(!is.na(node)) -> sankey.dtt
          max(as.numeric(gsub("level", "", sankey.dtt$next_x)), na.rm = TRUE) -> maxLev
          #sankey.dtt %>% group_by(x, node) %>% mutate(n = paste0("[", 100*n()/sankeyPrecision, "%]")) %>% rowwise() %>% mutate(level = as.numeric(gsub("level", "", x))) %>% mutate(label = ifelse(level == maxLev, paste(node, n), node)) -> sankey.dtt
          sankey.dtt %>% group_by(x, node) %>% mutate(n = paste0("[", 100*n()/sankeyPrecision, "%]")) %>% rowwise() %>% mutate(level = as.numeric(gsub("level", "", x))) %>% mutate(label = paste(node, n)) -> sankey.dtt

          sankey.dtt %>% rowwise() %>% mutate(node = ifelse(is.na(node), node, dealias(node))) %>% mutate(next_node = ifelse(is.na(next_node), next_node, dealias(next_node))) -> sankey.dtt

          sankey.dtt %>% rowwise() %>% mutate(label = ifelse((node == "0" & level>0), gsub("^0", "unclassified", label), label)) -> sankey.dtt

          unique(sankey.dtt$node) -> var2col
          if("0" %notin% var2col){c("0", var2col) -> var2col}
          unlist(lapply(as.list(var2col), dealias)) -> dealiasvar2col
          sort(dealiasvar2col) -> dealiasvar2col
          viridis_pal(alpha = 1, begin = 0.025, end = .975, direction = 1, option = "D")(length(dealiasvar2col)) -> var2col
          sankey.dtt <- sankey.dtt %>% mutate(label = ifelse((grepl("unclassified", label) & level != maxLev), NA, label))
          sankey.dtt <- sankey.dtt %>% group_by(node, x) %>% mutate(label = ifelse(length(unique(next_node)) == 1 & node == next_node & level != maxLev, NA, label))

          ggplot(sankey.dtt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = label)) +
                geom_sankey(flow.alpha = .6, type ='alluvial') +
                geom_sankey_label(size = 2, color = "white", fill = "gray40", position = position_nudge(x = 0.05, y = 0), na.rm = TRUE, type ='alluvial', hjust = 0) +
                theme_sankey(base_size = 16) +
                labs(x = NULL) +
                ggtitle("Gewichtetes Mittel: Österreich", subtitle = paste( sankey_date$earliest , "bis", sankey_date$latest)) +
                theme(legend.position = "none", axis.text.x = element_blank(), plot.title = element_text(hjust = 0), plot.subtitle=element_text(hjust = 0)) +
                scale_fill_manual(values = var2col, breaks = dealiasvar2col) +
                scale_x_discrete(expand = expansion(mult = c(0, .1), add = c(.1, .8))) -> pp

          filename <- paste0(outdir, "/figs/sankey/Overview_",  opt$country, ".pdf")
          ggsave(filename = filename, plot = pp, width = 6, height = 6)
          fwrite(as.list(c("sankey", "Overview", opt$country, filename)), file = summaryDataFile, append = TRUE, sep = "\t")
          rm(pp, filename, sankey.dt, sankey.dtt)


          writeLines(paste("PROGRESS: plotting WWTP detection plot"))

          reshape2::dcast(sankey_all.dt, LocationID~variant, value.var = "freq") -> sankey_completed.dt
          sankey_completed.dt[is.na(sankey_completed.dt)] <- 0
          sankey_completed.dt <- reshape2::melt(sankey_completed.dt, id.vars="LocationID", variable.name = "variant", value.name = "freq")

          left_join(x = sankey_completed.dt, y = (metaDT %>% dplyr::select(LocationID, connected_people) %>% distinct()), multiple = "all", by = "LocationID") %>%
              group_by(variant) %>% summarize(freq = weighted.mean(freq, w = connected_people, na.rm = TRUE), .groups = "keep") %>%
              filter(freq > zeroo/10) -> depicted_in_sankey

          occurence.dt <- data.table(variant = character(), d = numeric(), nd = numeric())
          for (variant_ in depicted_in_sankey$variant){
              variant_dealias_ <- dealias(variant_)
              sankey_all.dt %>% rowwise() %>%
                  mutate(variant_dealias = dealias(variant)) %>%
                  ungroup() %>%
                  mutate(N = length(unique(LocationID))) %>%
                  filter(grepl(variant_dealias_, variant_dealias)) %>%
                  filter(freq > zeroo/10) %>%
                  summarize(variant = paste0(variant_, ".*"), N = unique(N), d = length(unique(LocationID)), .groups = "keep") %>%
                  mutate(nd = N - d) %>% dplyr::select(-"N") -> occurence.rate.perVar
              occurence.dt <- rbind(occurence.dt, occurence.rate.perVar)
          }

          occurence.dt %>% mutate(r = paste0(sprintf("%.0f", 100*d/(d+nd)), "%")) -> occurence.rate
          reshape2::melt(occurence.dt, id.vars = c("variant"))  -> occurence.dt
          occurence.dt$variant <- factor(occurence.dt$variant, levels = occurence.dt %>% filter(variable == "d") %>% arrange(value) %>% pull(variant))

          dpl <- ggplot(data = occurence.dt, aes(x = variant, y = value, fill = variable))
          dpl <- dpl + geom_col(position = position_stack(reverse = TRUE))
          dpl <- dpl + scale_fill_manual(name = "", breaks = c("d", "nd"), values = brewer_pal(palette = "Set2")(2), labels = c("detektiert", "nicht detektiert"))
          dpl <- dpl + theme_bw()
          dpl <- dpl + xlab("Varianten")
          dpl <- dpl + ylab("Anzahl an Kläranlagen")
          dpl <- dpl + theme(legend.position="bottom")
          dpl <- dpl + scale_y_continuous()
          dpl <- dpl + theme(plot.title = element_text(hjust = 0), plot.subtitle=element_text(hjust = 0))
          dpl <- dpl + ggtitle("Anteil an positiven Kläranlagen per Variante", subtitle = paste( sankey_date$earliest , "bis", sankey_date$latest))
          dpl <- dpl + geom_label(data = occurence.rate, aes(x = variant, y = d, label = r), fill = "grey90",  size = 3)
          dpl <- dpl + coord_flip()

          filename1 <- paste0(outdir, "/figs/sankey/Detection_",  opt$country, ".pdf")
          ggsave(filename = filename1, plot = dpl, width = plotWidth, height = plotHeight)
          fwrite(as.list(c("detection", "Overview", opt$country, filename1)), file = summaryDataFile, append = TRUE, sep = "\t")

          writeLines(paste("PROGRESS: generating WWTP detection table"))

          sankey_all.dt %>%
              ungroup() %>% filter(variant %in% depicted_in_sankey$variant) %>%
              mutate(N = length(unique(LocationID))) %>%
              filter(freq > zeroo/10) %>%
              group_by(variant) %>% summarize(N = unique(N), d = length(unique(LocationID)), .groups = "keep") %>%
              mutate(nd = N - d) %>% dplyr::select(-"N") -> occurence.dt
          occurence.dt %>% mutate(r = paste0(sprintf("%.0f", 100*d/(d+nd)), "%")) -> occurence.rate
          reshape2::melt(occurence.dt, id.vars = c("variant"))  -> occurence.dt
          left_join(x = occurence.rate, y = occurence.freq, by = "variant", multiple = "all") -> synopsis.dt
          synopsis.dt %>% filter(freq > zeroo/10) %>% mutate(freq = signif(freq, digits = 2)) -> synopsis.dt
          colnames(synopsis.dt) <- c("Variante", "Detektiert", "Nicht detektiert", "Prozent", "Gew. Mittel")
          filename2 <- paste0(opt$dir, "/synopsis.tex")
          legendTxt <- paste("Für jede Virus Variente die Österreich in Proben von", sankey_date$earliest , "bis", sankey_date$latest ,"detektiert wurden, jeweils das gewichtete Mittel der relativen Häufigkeit und die Anzahl der Kläranlagen in denen die Variante detektiert wurde")

          makeTexTab(filename2, synopsis.dt, legendTxt)
          rm(dpl, filename1, filename2, occurence.dt, synopsis.dt, legendTxt)

        }
    }

    ## make plot per state
    if(dim(sankey_all.dt)[1] > 0 & any(colnames(metaDT) == "state")){
      left_join(x = sankey_all.dt, y = metaDT, by = "LocationID", multiple = "all") %>% dplyr::select(LocationID, variant, freq, state) %>% distinct() -> sankey_allState.dt
      sankey_state_plot_list <- list()
      for (stateoi in unique(sankey_allState.dt$state)){
        sankey_allState.dt %>% filter(state == stateoi) -> sankey_state.dt
        if(dim(sankey_state.dt)[1] > 0){

            globalFittedData %>% filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date) %>%
                filter(LocationID %in% unique(sankey_state.dt$LocationID)) %>%
                group_by(LocationID) %>%
                mutate(latest = max(sample_date, na.rm = TRUE)) %>%
                filter(sample_date == latest) %>% ungroup() %>%
                summarize(earliest = min(sample_date, na.rm = TRUE), latest = max(sample_date, na.rm = TRUE), .groups = "keep") -> sankey_date

            left_join(x = sankey_state.dt, y = (metaDT %>% dplyr::select(LocationID, connected_people) %>%
                distinct() %>% rowwise() %>%
                mutate(connected_people = ifelse(connected_people < 1, 1, connected_people))), multiple = "all", by = "LocationID") %>%
                group_by(variant) %>%
                summarize(freq = weighted.mean(freq, w = connected_people, na.rm = TRUE), .groups = "keep") -> sankey.dt

            sankey.dt  %>% rowwise() %>% mutate(freq = ifelse(freq<0, 0, freq)) -> sankey.dt
            sankey.dt  %>% rowwise() %>% mutate(freq = ifelse(freq>1, 1, freq)) -> sankey.dt

            if(sum(sankey.dt$freq)>1){
              sankey.dt  %>% ungroup() %>% mutate(freq = freq/sum(freq)) -> sankey.dt
            }
            sankey.dt %>% ungroup() -> sankey.dt
            sankey.dt -> occurence.freq

            sankey.dt %>% group_by(variant = "0") %>% summarize(freq = 1-sum(freq), .groups = "keep") -> sankey.dt0
            rbind(sankey.dt, sankey.dt0) -> sankey.dt
            sankey.dt %>% rowwise() %>% filter(freq > zeroo/10) -> sankey.dt

            sankey.dt$variant_dealias <- unlist(lapply(as.list(sankey.dt$variant), dealias))
            for ( baseForColor in baseForColorsVariants ){
              if(any(grepl(baseForColor, sankey.dt$variant_dealias, fixed = TRUE))){
                  if(baseForColor %notin% sankey.dt$variant_dealias){
                      rbind(sankey.dt, data.frame(variant = realias(baseForColor), freq = 0, variant_dealias = baseForColor)) -> sankey.dt
                  }
              }
            }

            sankey.dt %>% rowwise() %>% mutate(levels = length(unlist(strsplit(variant_dealias, split = "\\.")))) %>% group_by(1) %>% summarize( maxLev = max(levels), .groups = "keep") -> maxLev
            sankey.dt %>% rowwise() %>% mutate(long_variant = expandVar(variant_dealias, maxLev$maxLev)) -> sankey.dt
            ancestor <- getTree(sankey.dt)
            Nlevels <- max(unlist(lapply(ancestor, length)))

            sankey.dt <- makeObs(sankey.dt, Nlevels, ancestor, sankeyPrecision)

            sankey.dt %>% make_long(names(sankey.dt)) -> sankey.dtt

            sankey.dtt %>% filter(!is.na(node)) -> sankey.dtt
            max(as.numeric(gsub("level", "", sankey.dtt$next_x)), na.rm = TRUE) -> maxLev
            #sankey.dtt %>% group_by(x, node) %>% mutate(n = paste0("[", 100*n()/sankeyPrecision, "%]")) %>% rowwise() %>% mutate(level = as.numeric(gsub("level", "", x))) %>% mutate(label = ifelse(level == maxLev, paste(node, n), node)) -> sankey.dtt
            sankey.dtt %>% group_by(x, node) %>% mutate(n = paste0("[", 100*n()/sankeyPrecision, "%]")) %>% rowwise() %>% mutate(level = as.numeric(gsub("level", "", x))) %>% mutate(label = paste(node, n)) -> sankey.dtt

            sankey.dtt %>% rowwise() %>% mutate(node = ifelse(is.na(node), node, dealias(node))) %>% mutate(next_node = ifelse(is.na(next_node), next_node, dealias(next_node))) -> sankey.dtt

            sankey.dtt %>% rowwise() %>% mutate(label = ifelse((node == "0" & level>0), gsub("^0", "unclassified", label), label)) -> sankey.dtt
            sankey.dtt <- sankey.dtt %>% mutate(label = ifelse((grepl("unclassified", label) & level != maxLev), NA, label))
            sankey.dtt <- sankey.dtt %>% group_by(node, x) %>% mutate(label = ifelse(length(unique(next_node)) == 1 & node == next_node & level != maxLev, NA, label))

            unique(sankey.dtt$node) -> var2col
            if("0" %notin% var2col){c("0", var2col) -> var2col}
            unlist(lapply(as.list(var2col), dealias)) -> dealiasvar2col
            sort(dealiasvar2col) -> dealiasvar2col
            viridis_pal(alpha = 1, begin = 0.025, end = .975, direction = 1, option = "D")(length(dealiasvar2col)) -> var2col

            ggplot(sankey.dtt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = label)) +
                geom_sankey(flow.alpha = .6, type ='alluvial') +
                geom_sankey_label(size = 2, color = "white", position = position_nudge(x = 0.05, y = 0), na.rm = TRUE, type ='alluvial', hjust = 0) +
                theme_sankey(base_size = 16) + labs(x = NULL) +
                ggtitle(paste0("Gewichtetes Mittel: ", stateoi), subtitle = paste( sankey_date$earliest , "bis", sankey_date$latest)) +
                theme(legend.position = "none", axis.text.x = element_blank(), plot.title = element_text(hjust = 0), plot.subtitle=element_text(hjust = 0)) +
                scale_fill_manual(values = var2col, breaks = dealiasvar2col) +
                scale_x_discrete(expand = expansion(mult = c(0, .1), add = c(.1, .8))) -> pp
            sankey_state_plot_list[[length(sankey_state_plot_list)+1]] <- pp

            filename <- paste0(outdir, "/figs/sankey/Overview_",  gsub(" ", "_", stateoi), ".pdf")
            ggsave(filename = filename, plot = pp, width = 7, height = 7)
            fwrite(as.list(c("sankey", "Overview", stateoi, filename)), file = summaryDataFile, append = TRUE, sep = "\t")
            rm(pp, filename, sankey_state.dt, sankey.dt, sankey.dtt)
        }
      }

      plot_grid(plotlist = sankey_state_plot_list, ncol = 3) -> ppp
      filename <- paste0(outdir, "/figs/sankey/Overview_byState.pdf")
      ggsave(filename = filename, plot = ppp, width = 18, height = 18)
      fwrite(as.list(c("sankey", "Overview", "allStates", filename)), file = summaryDataFile, append = TRUE, sep = "\t")
      rm(ppp, filename)
    }
}

##


if (dim(globalFittedData)[1] == 0){
  writeLines(paste("WARNING: no varianted fitted in any of the samples; script will end here"))
  q("no")
}
globalFittedData %>% group_by(LocationID) %>% mutate(n = length(unique(sample_date))) %>% filter(n > tpLimitToPlot*2) %>% filter(variant %in% VoI) %>% ungroup() %>% group_by(LocationID, LocationName, sample_date, variant) %>% summarize(value = mean(value), .groups = "drop") -> globalFittedData2
if (dim(globalFittedData2)[1] >= tpLimitToPlot){
  writeLines(paste("PROGRESS: plotting overview VoI"))
  ggplot(data = globalFittedData2, aes(x = as.Date(sample_date), y = value, fill = variant)) + geom_area(position = "stack", alpha = 0.7)  + facet_wrap(~LocationName) + scale_fill_viridis_d(alpha = 0.6, begin = .05, end = .95, option = "H", direction = +1, name="") + scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks = c(0,0.5,1)) + theme_minimal() + theme(legend.position="bottom", strip.text.x = element_text(size = 4.5), panel.grid.minor = element_blank(), panel.spacing.y = unit(0, "lines"), legend.direction="horizontal") + guides(fill = guide_legend(title = "", ncol = 10)) + scale_x_date(date_breaks = "2 month", date_labels =  "%b") + ylab(paste0("Variantenanteil [1/1]") ) + xlab("") -> r2
  filename <- paste0(outdir, '/figs/fullview', paste('/wwtp', "VoI", sep="_"), ".pdf")
  #ggsave(filename = filename, plot = r2, width = plotWidth*1.5, height = plotHeight*1.5)
  #fwrite(as.list(c("fullOverview", c(roiname, "VoI", filename))), file = summaryDataFile, append = TRUE, sep = "\t")
  rm(r2,filename)
}
rm(globalFittedData2)

## print overview of all variants and all plants

globalFittedData$variant <- factor(globalFittedData$variant, levels = variant_order)
globalFittedData %>% group_by(LocationID) %>% mutate(n = length(unique(sample_date))) %>% filter(n > tpLimitToPlot*2) %>% ungroup() %>% group_by(LocationID, LocationName, sample_date, variant) %>% summarize(value = mean(value), .groups = "drop") -> globalFittedData2

# Generate map of each VoI detected recently

writeLines(paste("PROGRESS: start plotting maps"))
for (voi in VoI){
  writeLines(paste("  PROGRESS: considering", voi, "and sublineages"))

  voi_dealias <- dealias(voi)
  mapdata <- globalFittedData %>%
      filter(!is.na(variant)) %>%
      rowwise() %>%
      mutate(variant_dealias = dealias(variant)) %>%
      filter(grepl(voi_dealias, variant_dealias)) %>%
      group_by(LocationID, LocationName, sample_id, sample_date) %>%
      summarize(variant = paste0(voi, ".*"), value = sum(value), .groups = "keep")


  if(all(dim(mapdata) > 0)){
    metaDT$sample_date <- as.Date(metaDT$sample_date)
    mapdata$sample_date <- as.Date(mapdata$sample_date)
    left_join(x = mapdata, y = metaDT, by = c("sample_id" = "BSF_sample_name", "LocationName", "LocationID", "sample_date"), multiple = "all") -> mapdata

    if (dim(mapdata)[1] > 0){
      mapdata <- mapdata %>%
          group_by(LocationID) %>%
          mutate(latest = max(as.Date(sample_date))) %>%
          filter(latest == sample_date) %>%
          filter(BSF_run == last_BSF_run_id) %>%
          filter( (as.Date(format(Sys.time(), "%Y-%m-%d")) - recentEnought) < sample_date)
      mapdata <- mapdata %>% filter(value > 0)

      if (length(mapdata$value) > 0 ){
        writeLines(paste("  PROGRESS: plotting map for", paste0(voi,"*")))
        print(data.frame(location=mapdata$LocationName, dates = mapdata$sample_date, values = mapdata$value))

        s <- ggplot()
        s <- s + geom_sf(data = World, fill = "grey95")
        s <- s + geom_sf(data = Country, fill = "antiquewhite")
        s <- s + theme_minimal()
        s <- s + coord_sf(ylim = mapMargines[c(1,3)], xlim = mapMargines[c(2,4)], expand = FALSE)
        s <- s + annotation_scale(location = "bl")
        s <- s + annotation_north_arrow(location = "tl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering)
        s <- s + geom_point(data=(mapSeqDT %>% filter(status == "pass") ), aes(y=dcpLatitude, x=dcpLongitude, size = as.numeric(connected_people)), alpha = .4, fill = "grey80", shape = 21)
        s <- s + geom_point(data=(mapSeqDT %>% filter(status == "pass") ), aes(y=dcpLatitude, x=dcpLongitude, size = .7*as.numeric(connected_people)), alpha = .5, shape = 3)
        s <- s + geom_point(data=mapdata, aes(y=dcpLatitude, x=dcpLongitude, size = as.numeric(connected_people), col = as.numeric(value)), alpha = 0.8)
        s <- s + facet_wrap(~variant, nrow = 2)
        s <- s + theme(axis.text = element_blank(), legend.direction = "vertical", legend.box = "horizontal", legend.position = "bottom")
        s <- s + guides(size = guide_legend(title = "Population", nrow = 2), col = guide_colourbar(title = paste0("Anteil ", paste0(voi,"*")), direction = "horizontal", title.position = "top"))
        s <- s + scale_color_viridis_c(direction = 1, begin = .05, end = .95, option = "C")
        s <- s + scale_size(range = c(2, 6), labels = scales::comma, breaks = c(50000, 100000, 500000, 1000000))

        filename <- paste0(outdir, '/figs/maps/', voi, ".pdf")
        ggsave(filename = filename, plot = s, width = 7, height = 7)
        fwrite(as.list(c("mapplot", voi, filename)), file = summaryDataFile, append = TRUE, sep = "\t")
        rm(s, filename, mapdata)
      }
    }
  }
}

timestamp()
