## starting time stamp
timestamp()

## load libraries
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("betareg"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("rjson"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("NbClust"))
suppressPackageStartupMessages(library("gslnls"))
suppressPackageStartupMessages(library("cowplot"))



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
  make_option(c("--data"), type="character", default="data/mutationDataSub.tsv.g",
              help="Path to data input file in tidy table format [default %default]", metavar="character"),
  make_option(c("--precomp"), type="character", default="output-variants/globalFittedData.csv",
              help="Path to file with precomputed variant quantification. [default %default]", metavar="character"),
  make_option(c("--plotwidth"), type="double", default=8,
              help="Base size of plot width [default %default]", metavar="character"),
  make_option(c("--plotheight"), type="double", default=4.5,
              help="Base size of plot height [default %default]", metavar="character"),
  make_option(c("--TimecoursePlot"), type="character", default="last",
              help="Select for which weeks outbreak-Info Plot schould be produced. Allowed values are 'all', 'last', 'none' [default %default]", metavar="character"),
  make_option(c("--ninconsens"), type="double", default="0.4",
              help="Minimal fraction of genome covered by reads to be considered (0-1) [default %default]", metavar="character"),
  make_option(c("--zero"), type="double", default=0.01,
              help="Minimal allele frequency to be considered [default %default]", metavar="double"),
  make_option(c("--depth"), type="integer", default=50,
              help="Minimal depth at mutation locus to be considered [default %default]", metavar="character"),
  make_option(c("--removeLongIndels"), type="integer", default=9,
              help="Should indels with a length longer than <INT> nucleotides be ignored [default %default]", metavar="character"),
  make_option(c("--geogrowthlimit"), type="double", default=0.05,
              help="Mutations with smaller growth per week in geographic spread are ignored [default %default]", metavar="character"),
  make_option(c("--wmefgrowthlimit"), type="double", default=0.02,
              help="Mutations with smaller growth per week in weighted mean excess frequency are ignored [default %default]", metavar="character"),
  make_option(c("--alphaprime"), type="double", default=3.16,
              help="Alpha prime value, used together with lineage frequency to calculate Beta to model a Beta distribution with empirically observed variance [default %default]", metavar="character"),
  make_option(c("--pval"), type="double", default=0.05,
              help="Significants level for betareg test [default %default]", metavar="character"),
  make_option(c("--mingeospread"), type="double", default=0.15,
              help="Minimal fraction of positive wwtp per week in latest timepoint to evaluate a mutation [default %default]", metavar="character"),
  make_option(c("--minexcessfreq"), type="double", default=0.05,
              help="Minimal weighted mean excess frequency across all wwtp per week to evaluate a mutation [default %default]", metavar="character"),
  make_option(c("--minpositivesample"), type="integer", default=2,
              help="Only mutations with at least <INT> timepoints with a positive signal are considered [default %default]", metavar="character"),
  make_option(c("--minweeksample"), type="integer", default=6,
              help="Only weeks with at least <INT> samples are considered [default %default]", metavar="character"),
  make_option(c("--indels"), type="logical", default=FALSE,
              help="Should be indels be considered [default %default]", metavar="character"),
  make_option(c("--filstrat"), type="character", default="and",
              help="Filter strategy if mutations growing geographically AND/OR in mean excess frequency. Allowed values are 'or', 'and'. Other values will be ignored and defaulted to 'and' [default %default]", metavar="character"),
  make_option(c("--nw"), type="integer", default="2",
              help="Number of conecutive weeks with excess growth to be plotted [default %default]", metavar="character"),
  make_option(c("--start"), type="character", default="2000-01-01",
              help="Start of time course [default %default]", metavar="character"),
  make_option(c("--end"), type="character", default="2100-12-31",
              help="End of time course [default %default]", metavar="character"),
  make_option(c("--periodlength"), type="integer", default=2*30.5,
              help="Duration of analysis perid in days [default %default]", metavar="character"),
  make_option(c("--graph"), type="character", default="pdf",
              help="Fileformate of produced graphics. Select from pdf, png [default %default]", metavar="character"),
  make_option(c("--ignoreregion"), type="character", default="0-0",
              help="Regions to be ignored from analysis. Useful if certain amplicons are producing unreliable AF. Specify in the a format like <100-120|21955-23288> to specify two regions. [default %default]", metavar="character"),
  make_option(c("--onlyLastRun"), type="logical", default="TRUE",
              help="Set to FALSE if all samples and not just the samples which where included in the last run should be considered. [default %default]", metavar="character"),
  make_option(c("--debug"), type="logical", default="FALSE",
              help="Toggle to use run with provided example data [default %default]", metavar="character"),
  make_option(c("--verbose"), type="logical", default="FALSE",
              help="Toggle to report more output, in terms of log messages but also graphical content, which slows down run time significantly [default %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#####################################################
####  parameter setting for interactive debugging
if(opt$debug){
  opt$dir = "output-mutations"
  opt$metadata = "data/metaData_general.csv"
  opt$data="data/mutationData_DB_NationMonitoringSites.tsv.gz"
  opt$inputformat = "tidy"
  opt$marker="VaQuERo/resources/mutations_list_grouped_pango_codonPhased_2024-09-04_Europe.csv"
  opt$mutstats  = "VaQuERo/resources/mutations_stats_pango_codonPhased_2024-09-04.csv.gz"
  opt$group2var = "VaQuERo/resources/groupMembers_pango_codonPhased_2024-09-04_Europe.csv"
  opt$pmarker="VaQuERo/resources/mutations_problematic_2023-11-23.csv"
  opt$precomp = "output-variants/globalFittedData.csv"
  #opt$ninconsens = 0.2
  opt$zero=0.01
  opt$depth=500
  opt$pval = 0.01
  opt$minpositivesample=2
  opt$removeLongIndels=1
  opt$periodlength = 61
  opt$indels = FALSE
  opt$verbose = FALSE
  opt$graph = "pdf"
  opt$TimecoursePlot="all"
  opt$geogrowthlimit = 0.04
  opt$wmefgrowthlimit = 0.04
  opt$mingeospread=0.04
  opt$minexcessfreq=0.04
  opt$nw   = 2
  opt$alphaprime = 2.2
  opt$minweeksample = 8
  opt$filstrat = "and"
  opt$start = Sys.Date() - 6*30
  opt$end = Sys.Date()
  opt$onlyLastRun = FALSE
  #opt$ignoreregion = "21955-23288"
  writeLines("Warning: command line option overwritten")
}
#####################################################

## define functions
options(warn=-1)

## loading external function library, expected to be in same dir as execution script
sourceFileBase = "VaQueR_functions.r"
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
  sourceFile <- list.files(
      ".",
      pattern = sourceFileBase,
      recursive = TRUE,
      full.names = TRUE
    )
}

if(length(sourceFile) == 1 && file.exists(sourceFile)){
  writeLines(paste("LOG: loading function source file", sourceFile))
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
    writeLines(paste("LOG: loading function source file", sourceFile))
    opt$source = sourceFile
    source(sourceFile)
  } else if(length(sourceFile) > 1 && file.exists(sourceFile[1])){
    writeLines(paste0("WARNING: several function source files found (", paste(sourceFile, collapse = ", ", sep = ","), ")"))
    writeLines(paste("LOG: loading function source file", sourceFile[1]))
    opt$source = sourceFile[1]
    source(sourceFile[1])
  } else{
    writeLines(paste("ERROR: source file", sourceFileBase, "not found. Please double check if it is in the same directory as the analysis script VaQuERo_v2.R"))
    quit(save="no")
  }
}

## writeLines parameter to Log
writeLines("##~LOG~PARAMETERS~####################")
print(opt)
writeLines("##~LOG~PARAMETERS~####################")
writeLines("\n\n\n")

## read in config file to overwrite all para below

## set variable parameters
outdir = opt$dir
summaryDataFile <- paste0(outdir, "/summary.csv")
markermutationFile  <- opt$marker
mutationstatsFile  <- opt$mutstats
group2variantFile  <- opt$group2var
problematicmutationFile <- opt$pmarker
alphaprime <- opt$alphaprime
plotWidth  <- opt$plotwidth
plotHeight <- opt$plotheight
minpossamp_th <- opt$minpositivesample
pval_th <- opt$pval
wmef_th <- opt$minexcessfreq
r_th <- opt$mingeospread
ignoreregionstart <- c()
ignoreregionend <- c()

if (opt$ignoreregion != "0-0"){
  ignoreregions <- unlist(strsplit(opt$ignoreregion, "\\|"))
  for (ir in ignoreregions){
    if(grepl("^\\d+-\\d+$", ir, perl = TRUE)){
        irp <- unlist(strsplit(ir, "-"))
        ignoreregionstart <- c(ignoreregionstart, irp[1])
        ignoreregionend <- c(ignoreregionend, irp[2])
    } else{
      writeLines("Warning: comand line defined regions to be ignored <ir> does not fit format like <start-stop>; region will not be ignored. please double check.")
    }
  }
  writeLines("LOG: mutations in following regions will not be considered for analysis:")
  for (i in 1:length(ignoreregionstart)){
    writeLines(paste(ignoreregionstart[i], "-", ignoreregionend[i]))
  }
   writeLines("")
}



## create directory to write plots
if(opt$verbose){
  writeLines(paste0("PROGRESS: create directory "))
}
if( ! dir.exists(outdir)){
  dir.create(outdir, showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs"))){
  dir.create(paste0(outdir, "/figs"), showWarnings = FALSE)
}
if( opt$verbose & !dir.exists(paste0(outdir, "/figs/singleModels"))){
  dir.create(paste0(outdir, "/figs/singleModels"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/nucleotides"))) {
  dir.create(paste0(outdir, "/figs/nucleotides"), showWarnings = FALSE)
}
file.create(summaryDataFile)

## read alias file
aliases <- fromJSON(file ="https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json")
as.list(names(aliases)) -> dealiases
names(dealiases) = as.character(aliases)

## set definition of failed/detected/passed
### how many none-N in consensus may be seen to be called detected
### how many non-N in consensus may be seen to be called "pass"
N_in_Consensus_filter <- 29903 - 29781 * opt$ninconsens

zeroo <- opt$zero      # marker below this value are considered to be zero
min.depth = opt$depth     # mutations with less read support are ignored
minUniqMarker <- opt$minuniqmark # minmal absolute number of uniq markers that variant is considered detected
minUniqMarkerRatio <- opt$minuniqmarkfrac # minmal fraction of uniq markers that variant is considered detected (0.25)
minMarker <- opt$minqmark # minmal absolute number of sensitive markers that variant is considered detected
minMarkerRatio <- opt$minmarkfrac # minmal fraction of sensitive markers that variant is considered detected (0.25)
colorSets <- c("Blues", "Greens", "Purples", "Oranges", "Reds", "Greys")
removeLongIndels <- opt$removeLongIndels


## define global variables to fill
globalAFdataList  <- list()
globalAFdata      <- data.table(sampleID = character(), LocationID = character(), LocationName  = character(), sample_date = character(), nuc_mutation = character(), aa_mutation = character(), observed = numeric(), expected = numeric(), excess = numeric(), pvalue = numeric() )

## read mutation stats
if(opt$verbose){
  writeLines(paste0("PROGRESS: read mutation per variant stats "))
}
mstat <- fread(file = mutationstatsFile)
group2var <- fread(file = group2variantFile)
left_join(x = mstat, y = group2var, by = c("ID" = "variants_alias")) %>% group_by(pos, ref, alt, AA_change, nucc, groupName) %>% summarize(value = sum(value, na.rm = TRUE), count = sum(count, na.rm = TRUE), .groups = "keep") %>% ungroup() %>% mutate(sensitivity = value/count) -> mstatgroup
if (removeLongIndels > 0){
  mstatgroup %>% filter(nchar(alt) <= removeLongIndels & nchar(ref) <= removeLongIndels) -> mstatgroup
}

## read pre computed variant quantifications.
if(file.exists(opt$precomp)){
    precomputed_var <- fread(opt$precomp)
} else{
    writeLines(paste0("ERROR: expecting a file with variant quantification provided via <--precomp FILE>; but specified (or default) file '", opt$precomp, "' does not exist"))
    quit(save="no")
}


## read problematic mutations from file
if(opt$verbose){
  writeLines(paste0("PROGRESS: read file with problematic mutation sites, to be ommited throughout the analysis "))
}
poi <- fread(file = problematicmutationFile)
unite(poi, NUC, c(4,3,5), ALT, sep = "", remove = FALSE) -> poi
poi %>% dplyr::select("Postion","NUC") -> problematic

## read in meta data
if(opt$verbose){
  writeLines(paste0("PROGRESS: read and process meta data "))
}
metaDT       <- fread(file = opt$metadata)
unique(metaDT$BSF_sample_name) -> sampleoi

# take care of "same date, same location" issue
# introduce artificial decimal date
metaDT %>% group_by(LocationID, sample_date) %>% mutate(n = rank(BSF_sample_name)-1)  %>% ungroup() %>% rowwise() %>% mutate(sample_date_decimal = decimalDate(sample_date, n)) %>% dplyr::select(-n) -> metaDT

# remove sample before start and after end date specified via command line
metaDT %>% filter(sample_date >= opt$start) %>% filter(sample_date <= opt$end) -> metaDT
opt$start <- min(metaDT$sample_date)
opt$end <- max(metaDT$sample_date)


## determine last sequencing batch ID
## based on BSF_start_date specified in metadata
as.character(metaDT %>% filter(!is.na(BSF_run)) %>% filter(!is.na(BSF_start_date)) %>% group_by(BSF_run, BSF_start_date) %>% summarize(.groups = "keep") %>% ungroup() %>% filter(as.Date(BSF_start_date) == max(as.Date(BSF_start_date))) %>% dplyr::select("BSF_run")) -> last_BSF_run_id
writeLines(paste("LOG: last_BSF_run_id", last_BSF_run_id))


## read in mutations data
if(opt$verbose){
  writeLines(paste0("PROGRESS: read AF data "))
}
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

## remove problematic sites
sewage_samps.dt %>% filter(POS %notin% problematic$Postion) -> sewage_samps.dt

## removeLongIndels
if (removeLongIndels > 0){
  befi <- length(unique(sewage_samps.dt$NUC))
  sewage_samps.dt %>% filter(nchar(ALT) <= removeLongIndels & nchar(REF) <= removeLongIndels) -> sewage_samps.dt
  afi <- length(unique(sewage_samps.dt$NUC))
  writeLines(paste0("LOG: indels of length greater than or equal to ",removeLongIndels, " nt are removed; in total ", befi-afi, " mutations were ignored; set --removeLongIndels 0 if undesired."))
}

## set AF := 0 if DP <= min.depth
sewage_samps.dt <- sewage_samps.dt %>% filter(value.depth >= min.depth)

# remove regions specified via --ignoreregion
if(length(ignoreregionend) > 0){
    mnb <- length(unique(sewage_samps.dt$NUC))
    for (ir in (1:length(ignoreregionend))){
      sewage_samps.dt <- sewage_samps.dt %>% filter(POS < ignoreregionstart[ir] | POS > ignoreregionend[ir])
    }
    mna <- length(unique(sewage_samps.dt$NUC))
    writeLines(paste("LOG: a total of", mnb-mna, "(out of", mnb, ") were removed since they overlapped with a region as specified with the --ignoreregion parameter"))
}

## add location to sewage_samps.dt
sewage_samps.dt %>% mutate(RNA_ID_int = gsub("_S\\d+$","", ID)) -> sewage_samps.dt
sewage_samps.dt %>% mutate(RNA_ID_int = gsub("_\\d$","", RNA_ID_int))  -> sewage_samps.dt

if(any(is.na(sewage_samps.dt$LocationID))){
  writeLines(paste("WARNING: samples without assigned location were found and removed."))
  sewage_samps.dt %>% filter(!is.na(LocationID)) -> sewage_samps.dt
}

metaDT %>% filter(! is.na(LocationID) ) %>% dplyr::select("RNA_ID_int", "LocationID", "LocationName") -> sample_location
left_join(x = sewage_samps.dt, y = sample_location, by = "RNA_ID_int") -> sewage_samps.dt

## remove samples which are not include_in_report == TRUE
metaDT %>% filter(is.na(include_in_report) | include_in_report == TRUE ) %>% dplyr::select("RNA_ID_int", "BSF_sample_name") -> sample_includeInReport
sewage_samps.dt %>% filter(ID %in% sample_includeInReport$BSF_sample_name) -> sewage_samps.dt

## remove samples with >N_in_Consensus_filter N
metaDT  %>% filter(N_in_Consensus < N_in_Consensus_filter) %>% dplyr::select("RNA_ID_int", "BSF_sample_name") -> passed_samples
sewage_samps.dt <- sewage_samps.dt[sewage_samps.dt$ID %in% passed_samples$BSF_sample_name,]

## add sampling date
metaDT %>% dplyr::select("RNA_ID_int", "sample_date", "sample_date_decimal") -> sample_dates
left_join(x=sewage_samps.dt, y=sample_dates, by = "RNA_ID_int") -> sewage_samps.dt

## get most recent sampling date from last Run
metaDT %>% filter(BSF_start_date == sort(metaDT$BSF_start_date)[length(sort(metaDT$BSF_start_date))]) %>% dplyr::select("BSF_run", "RNA_ID_int", "BSF_sample_name", "sample_date") -> RNA_ID_int_currentRun

RNA_ID_int_currentRun %>% ungroup() %>% summarize(latest = max(as.Date(sample_date)), .groups = "keep") -> latestSample
RNA_ID_int_currentRun %>% ungroup() %>% summarize(earliest = min(as.Date(sample_date)), .groups = "keep") -> earliestSample

writeLines(paste("LOG: current run ID:", unique(RNA_ID_int_currentRun$BSF_run)))
writeLines(paste("LOG: earliest sample in current run:", earliestSample$earliest))
writeLines(paste("LOG: latest sample in current run:", latestSample$latest))

## generate for each mutation a nice printable label
sewage_samps.dt %>% mutate(label = paste(NUC, paste0("[", paste(ANN.GENE,ANN.AA, sep =":"), "]"))) %>% dplyr::select(NUC, label) %>% distinct() -> nuc2label


## check for each mutation if DP explains increase better then time
if(FALSE){
  DPartefakts <- c()
  for (mif in unique(sewage_samps.dt$NUC)) {
    mif.dt <- sewage_samps.dt %>% filter(NUC == mif)

    if(length(unique(mif.dt$LocationID)) > 1){
      modeldpt <- lm(log10(value.freq) ~ LocationID + log10(value.depth) + sample_date_decimal, data = mif.dt)
      modeldp <- lm(log10(value.freq) ~ LocationID + log10(value.depth), data = mif.dt)

      modelcmp <- anova( modeldp, modeldpt, test="Chisq")
      mif.pval <- modelcmp$`Pr(>Chi)`[2]
      mif.pval <- ifelse(is.na(mif.pval), 1, mif.pval)
      if(mif.pval > pval_th){
        DPartefakts <- c(DPartefakts, mif)
      }
    } else{
      DPartefakts <- c(DPartefakts, mif)
    }
  }
  sewage_samps.dt <- sewage_samps.dt %>% filter(NUC %notin% DPartefakts)
}

## kill if no mutations are found
if (length(unique(sewage_samps.dt$ID)) < 1){
  writeLines(paste("WARNING: no mutation data found in input file. Program ends here!"))
  quit(save="no")
}

### FROM HERE LOOP OVER EACH SAMPLE LOCATION
writeLines(paste("PROGRESS: start to loop over sample location"))
# r <- grep("ATTP_9", unique(sewage_samps.dt$LocationID))
for (r in 1:length(unique(sewage_samps.dt$LocationID))) {
    roi = unique(sewage_samps.dt$LocationID)[r]
    roiname = unique(sewage_samps.dt$LocationName)[r]
    if(opt$verbose){
      writeLines(paste("PROGRESS: start processing", roiname, roi, format(Sys.time(), "%X %Y-%m-%d")))
    }

    ### filter data set for current location
    sewage_samps.dt %>% filter(LocationID == roi) %>% dplyr::select("ID", "sample_date", "sample_date_decimal", "value.freq", "value.depth", "LocationID", "LocationName", "NUC") -> sdt

    ## count how many sample from this plant were in the last BSF run
    metaDT %>% filter(BSF_sample_name %in% sdt$ID) %>% filter(BSF_run %in% last_BSF_run_id) %>% dplyr::select(BSF_run) %>% group_by(BSF_run) %>% summarize(n = n(), .groups = "keep") -> count_last_BSF_run_id

    ## skip roi if not in the most recent run
    if(opt$onlyLastRun == TRUE & identical(count_last_BSF_run_id$n, integer(0))){
        writeLines(paste0("WARNING: sample location <", roiname, "> skipped from analysis since not part of the latest sequencing run <", last_BSF_run_id, ">"))
        next;
    }

    ### find excessive mutation frequencies in each timepoint
    sdt %>% dplyr::select(sample_date_decimal, sample_date) %>% distinct() -> timePointsCombinations
    timePoints <- timePointsCombinations$sample_date_decimal[order(timePointsCombinations$sample_date_decimal)]
    timePoints_classic <- timePointsCombinations$sample_date[order(timePointsCombinations$sample_date_decimal)]

    for (t in 1:length(timePoints)) {
            timepoint <- timePoints[t]
            timepoint_classic <- timePoints_classic[t]
            ref_timepoint_classic <- timepoint_classic
            timepoint_day <- decimalDate(timepoint_classic,0)

            sdt %>% filter(sample_date_decimal == timepoint ) -> ssdt
            sampleID <- paste0(unique(ssdt$ID))
            sewage_samps.dt %>% filter(ID == sampleID) -> sewage_samps.sdt
            if(opt$verbose){
              writeLines(paste("PROGRESS:", sampleID, "@", roiname, "@", timepoint_classic, "(", signif(timepoint, digits = 10), ")"))
            }

            ## get specified lineage from pre computed values externally provided
            specifiedLineages <- precomputed_var %>% filter(LocationName == roiname) %>% filter(sample_date == timepoint_classic) %>% filter(value > 0) %>% pull(variant) %>% unique()
            sampleFittedData <- precomputed_var %>% filter(LocationName == roiname) %>% filter(sample_date == timepoint_classic) %>% mutate(norm.value = value) %>% filter(sample_id != "")
            if(length(specifiedLineages) == 0){
              writeLines(paste0("WARNING: there are no lineages detected for ", roiname, " @ ", timepoint_classic, "; this might be due to bad quality of the sample, wrong use of variant quantification input file (", opt$precomp, "), or since a unknown variant is already dominant. this sample will be skipped, please manually double check!"))
              next;
            }
            writeLines(paste("LOG: (", t, ")", timepoint_classic, roiname, paste("(", length(specifiedLineages), "detected lineages)")))

            ## further analyse if successfully deduced variants frequencies
            if(!exists("sampleFittedData")){
                writeLines(paste("LOG: no variant deduced", sampleID, roiname, timePoints_classic[t]))
                next;
            }
            if(opt$verbose){
              writeLines(paste("PROGRESS: infer excess mutations", sampleID, roiname, timePoints_classic[t]))
            }

            ##################
            ## infer expected AF based on AF per variant (as seen in pango) and deduced variant frequency
            full_join(x = (sampleFittedData %>% filter(value > 0)) , y = (mstatgroup %>% filter(groupName %in% sampleFittedData$variant)), by = c("variant" = "groupName"), suffix = c(".var", ".mut")) %>% dplyr::select(variant, value.var, sensitivity, AA_change, nucc) %>%  ungroup() %>% mutate(varAFcontribution = value.var*sensitivity) %>% group_by(nucc, AA_change) %>% summarize(AF = sum(varAFcontribution, na.rm = TRUE), .groups = "keep") %>% mutate(AF = ifelse(AF > 1, 1, AF)) -> expectedAF

            ## extract observed AF
            sewage_samps.sdt %>% filter(ID == sampleID) -> observedAF
            right_join(x= expectedAF, y = observedAF, by = c("nucc" = "NUC")) %>% dplyr::select(nucc, AA_change, AF, value.freq, value.depth) %>% rename(expected.freq = AF, observed.freq = value.freq) -> expected_vs_observed_AF

            ##################
            ## contrast observed and expected AF
            expected_vs_observed_AF %>% mutate(expected.freq = ifelse(is.na(expected.freq), 0, expected.freq)) %>% mutate(observed.freq = ifelse(is.na(observed.freq), 0, observed.freq)) -> expected_vs_observed_AF

            ## filter sig. excess mutations based on beta distriubtion
            # transform 0 and 1 values to [0,1]
            refDepth <- min(expected_vs_observed_AF$value.depth) # was set to 5000 at some point
            minRefFreq <- qbeta(c(pval_th), alphaprime, betaParamFromMean((1-(refDepth - 1 + 0.5)/refDepth), alphaprime), ncp = 0, lower.tail = FALSE, log.p = FALSE)
            expected_vs_observed_AF  %>%
                filter(expected.freq > minRefFreq | observed.freq > minRefFreq) %>%
                mutate(expected.freq.trans = ifelse(expected.freq < (1-(refDepth - 1 + 0.5)/refDepth), (1-(refDepth - 1 + 0.5)/refDepth), expected.freq)) %>%
                mutate(expected.freq.trans = (expected.freq.trans * (refDepth - 1 + 0.5)/refDepth) ) %>%
                mutate(observed.freq.trans = ifelse(observed.freq < (1-(refDepth - 1 + 0.5)/refDepth), (1-(refDepth - 1 + 0.5)/refDepth), observed.freq)) %>%
                mutate(observed.freq.trans = (observed.freq.trans * (refDepth - 1 + 0.5)/refDepth) ) -> expected_vs_observed_AF_transformed

            # test if observed and expected frequencies are congruent
            # using one.sided beta distribution test
            expected_vs_observed_AF_transformed %>%
                mutate(pval = pbeta(observed.freq.trans, alphaprime, betaParamFromMean(expected.freq.trans, alphaprime), ncp = 0, lower.tail = FALSE, log.p = FALSE)) -> expected_vs_observed_AF_transformed

            # adjust for multiple testing
            expected_vs_observed_AF_transformed$qval <- p.adjust(expected_vs_observed_AF_transformed$pval, method = "fdr")

            # filter mutations based on significance level pval_th
            expected_vs_observed_AF_transformed %>% filter(qval < pval_th) -> expected_vs_observed_AF_transformed_filtered

            ## collect globally data on expected mutations
            if(dim(expected_vs_observed_AF_transformed_filtered)[1]>0){
                localAFdata <- data.table(
                                  sampleID      = sampleID,
                                  LocationID    = roi,
                                  LocationName  = roiname,
                                  sample_date   = as.character(timepoint_classic),
                                  nuc_mutation  = expected_vs_observed_AF_transformed$nucc,
                                  aa_mutation   = expected_vs_observed_AF_transformed$AA_change,
                                  observed      = expected_vs_observed_AF_transformed$observed.freq,
                                  expected      = expected_vs_observed_AF_transformed$expected.freq,
                                  excess        = expected_vs_observed_AF_transformed$observed.freq - expected_vs_observed_AF_transformed$expected.freq,
                                  pvalue        = expected_vs_observed_AF_transformed$qval,
                                  dp            = expected_vs_observed_AF_transformed$value.depth
                              )
                globalAFdataList[[length(globalAFdataList)+1]] <- localAFdata

            }
    }
}
globalAFdata <- rbindlist(globalAFdataList)

writeLines(paste("PROGRESS: loop over WWTP finished"))


if(dim(globalAFdata)[1] == 0){
    writeLines(paste("WARNING: no mutation data found in excess. Program ends here!"))
    quit(save="no")
}

####################
## add meta data and consistend midweek data to globalAFdata

writeLines(paste("PROGRESS: add meta data and consistend midweek data to globalAFdata"))

metaDT %>% dplyr::select(LocationID, LocationName, connected_people, dcpLatitude, dcpLongitude) %>% distinct() -> wwtp_meta
metaDTs <- metaDT %>% dplyr::select(LocationID, LocationName, connected_people, sample_date) %>% mutate(midweek_date = sample_date + 3 - as.numeric(strftime(sample_date, format = "%u"))) %>% distinct()
metaDTs <- metaDTs %>% group_by(midweek_date) %>% mutate(n = n()) %>% filter(n > opt$minweeksample) %>% dplyr::select(-n, -sample_date)
metaDTs$midweek_date <- as.Date(metaDTs$midweek_date)

globalAFdata_m <- left_join(x = globalAFdata, y = wwtp_meta, by = c("LocationID", "LocationName"))
globalAFdata_m$sample_date <- as.Date(globalAFdata_m$sample_date, tryFormats = c("%Y-%m-%d"))
globalAFdata_m <- globalAFdata_m %>% mutate(midweek_date = sample_date + 3 - as.numeric(strftime(sample_date, format = "%u")))


####################
## loop over time and deduce growing excess mutation till there
writeLines(paste("PROGRESS: start to loop over time course"))
timecourseMutationSet_List <- list()
timecourseMutationSet <- data.table(nuc_mutation = character(), label = character(), midweek_date = character(), r = numeric(), excess_mw = numeric(), growth_pred.geo = numeric(), growth_pred.wmef = numeric())

for (periodend in as.character(seq( from = min(globalAFdata_m$midweek_date) + opt$periodlength, to = max(globalAFdata_m$midweek_date), by = 7))){

    periodend <- as.Date(periodend)
    periodstart <- as.Date(periodend-opt$periodlength)
    if(opt$verbose){
      writeLines(paste("PROGRESS: processing", periodstart, "to", periodend))
    }

    currentAFdata <- globalAFdata_m %>% filter(midweek_date > periodstart & midweek_date <= periodend)
    currentmetaDTs <- metaDTs %>% filter(midweek_date > periodstart & midweek_date <= periodend)

    ### next if no excess mutations are found
    if (dim(currentAFdata)[1] < 1){
      writeLines(paste("WARNING: no excess mutations found; skipping timepoint", periodend))
      next
    }

    ####################
    ## complete record so that for each week with measurments all locations are add with observation 0
    ## remove weeks with less then opt$minweeksample samples
    ## set none significant excess values to 0
    currentAFdata <- currentAFdata %>% mutate(excess = ifelse(is.na(excess), 0, excess))
    currentAFdata <- currentAFdata %>% rowwise() %>% mutate(excess = ifelse(pvalue < pval_th, excess, 0))
    currentAFdata_completed <- data.table::dcast(currentAFdata, nuc_mutation+midweek_date~LocationID, fun.aggregate = mean, value.var = "excess", fill = 0) %>% data.table::melt(id.vars = c("nuc_mutation", "midweek_date"), variable.name = "LocationID", value.name = "excess")
    currentDPdata_completed <- data.table::dcast(currentAFdata, nuc_mutation+midweek_date~LocationID, fun.aggregate = mean, value.var = "dp", fill = 0) %>% data.table::melt(id.vars = c("nuc_mutation", "midweek_date"), variable.name = "LocationID", value.name = "dp")
    currentAFdata_completed <- inner_join(x = currentAFdata_completed, y = currentmetaDTs, by = c("LocationID", "midweek_date"))
    currentAFdata_completed <- left_join(x = currentAFdata_completed, y = currentDPdata_completed, by = c("LocationID", "midweek_date", "nuc_mutation"))


    ####################
    ## calculate population weighted mean excess mutation frequency per mutation
    if(opt$verbose){
      writeLines(paste("PROGRESS: calculate population weighted mean excess mutation frequency per mutation"))
    }
    dt_wmef <- currentAFdata_completed %>%
        group_by(nuc_mutation, midweek_date) %>%
        summarize(
          excess_mw = weighted.mean(excess, w = connected_people),
          dp_mw = weighted.mean(dp[dp>0], w = connected_people[dp>0]),
          .groups = "keep"
        )

    ####################
    ## deduce nls growth of timecourse population weighted mean excess mutation frequency per mutation
    ## filter mutations with a max weighed mean excess freq > wmef_th and with at least minpossamp_th weighed mean excess freq > 0 across all timepoints
    if(opt$verbose){
      writeLines(paste("PROGRESS: deduce nls growth of timecourse population weighted mean excess frequency"))
    }
    ind.mois <- dt_wmef %>%
      group_by(nuc_mutation) %>%
      mutate(max = max(excess_mw)) %>%
      mutate(n = sum(excess_mw > 0)) %>%
      rowwise() %>%
      filter(max >= wmef_th & n >= minpossamp_th) %>%
      pull(nuc_mutation) %>%
      unique()

    prediction_collector_wmef_List <- list()
    prediction_collector_wmef <- data.table(mutation = character(), growth_pred = numeric(), inflection_pred = numeric())

    for (ind.moi in ind.mois){
      dtoi <- dt_wmef %>% filter(nuc_mutation == ind.moi)
      dtoi_weights <- dtoi$dp_mw
      colnames(dtoi)[2] <- "t"
      colnames(dtoi)[3] <- "y"
      subset(dtoi, select = c("t", "y")) -> dtoi
      dtoi$d <- as.numeric(as.Date(dtoi$t)-as.Date(dtoi$t[1]))
      a <- 1
      fo <- y ~ a / (1 + exp(-b * (d-c)))
      model <- tryCatch(gsl_nls(fo, start = list(b = 0.05, c = 2*max(dtoi$d)/3), data = dtoi, weights = dtoi_weights),error=function(e) NA, warning=function(w) NA)
      # if gsl_nls regession did not converge try different starting parameters
      if(1){
          if(identical(NA, model)){
            model <- tryCatch(gsl_nls(fo, start = list(b = 0.01, c = max(dtoi$d)/2), data = dtoi, weights = dtoi_weights),error=function(e) NA, warning=function(w) NA)
          }
          if(identical(NA, model)){
            model <- tryCatch(gsl_nls(fo, start = list(b = 0.1, c = 2*max(dtoi$d)/3), data = dtoi, weights = dtoi_weights),error=function(e) NA, warning=function(w) NA)
          }
          if(identical(NA, model)){
            model <- tryCatch(gsl_nls(fo, start = list(b = 0.2, c = max(dtoi$d)/2), data = dtoi, weights = dtoi_weights),error=function(e) NA, warning=function(w) NA)
          }
      }

      if(identical(NA, model)){
         writeLines(paste("WARNING: no convergence for", ind.moi, "assessing weighted AF mean"))
      } else{
         # predict the relative weekly growth in the last 7 days
         tdays <- 7
         dtoi$predict <- predict(model)
         lastbiweek_prediction <- predict(model, newdata = data.frame(t = c(as.Date(periodend)-tdays, periodend), d = as.numeric(c(as.Date(periodend)-tdays, periodend)-as.Date(dtoi$t[1]))))
         lastweek_growth <- (lastbiweek_prediction[2] - lastbiweek_prediction[1])/(tdays/7)
         lastweek_growth -> growth_pred
         coefficients(model)[2] -> inflection_pred

         ## set inflection_pred and growth_pred to NA growth_pred confidence interval overlaps 0
         summary(model) -> model_summary
         model_summary$coefficients[ , 4][1] -> pvalue_growth_pred
         model_summary$coefficients[ , 4][2] -> pvalue_inflection_pred
         CI_growth_pred <- tryCatch(suppressMessages(confintd(model, "b", level = 0.90)),error=function(e) NA, warning=function(w) NA)
         CI_inflection_pred <- tryCatch(suppressMessages(confintd(model, "c", level = 0.90)),error=function(e) NA, warning=function(w) NA)

         if (any(is.na(CI_growth_pred)) | (min(CI_growth_pred[2:3]) < 0 & max(CI_growth_pred[2:3]) > 0) | growth_pred < opt$wmefgrowthlimit){
           growth_pred = NA
           inflection_pred = NA
         } else{
            if (!is.na(growth_pred) & opt$verbose == TRUE){
              melt(dtoi[,-3], id.vars = "t")  %>% mutate(t = as.Date(t)) %>%
                  ggplot(aes(x = t, y = value, color = variable)) +
                  geom_point(size = 3, alpha = 0.66, aes(shape = variable, size = variable)) +
                  geom_smooth(formula = y ~ x, method = "glm", method.args = list(family = "binomial"), se = FALSE) +
                  theme_bw() +
                  ggtitle(paste("Till", periodend, nuc2label$label[nuc2label$NUC == ind.moi], sep = ": ")) +
                  theme(legend.position="bottom", legend.title = element_blank()) +
                  ylab("mutation frequency [1/1]") +
                  xlab("Sample Date") +
                  #geom_vline(xintercept=as.Date(dtoi$t[1])-dtoi$d[1]+as.numeric(inflection_pred), linetype = "dotted")  +
                  scale_color_brewer(palette = "Set2", breaks = c("y", "predict"), labels = c("Observations", "Model Prediction")) +
                  scale_shape(breaks = c("y", "predict"), labels = c("Observations", "Model Prediction")) +
                  annotate(x = median(dtoi$t), y = max(dtoi$predict, dtoi$y), geom = "text", label = paste("Growth per week:\n", signif(growth_pred, digits = 2))) -> p
              ggsave(filename = paste0(outdir, "/figs/singleModels/", "plot_wmef_", periodend, "_", ind.moi, ".", opt$graph), p, width = 5, height = 5)
            }
        }
        prediction_wmef <-  data.table(mutation = ind.moi, growth_pred = growth_pred, inflection_pred = inflection_pred)
        prediction_collector_wmef_List[[length(prediction_collector_wmef_List)+1]] <- prediction_wmef
      }
      rm(dtoi, model)
    }
    prediction_collector_wmef <- rbindlist(prediction_collector_wmef_List)
    wmef_growing_nucs <- prediction_collector_wmef %>% filter(!is.na(growth_pred))

    ####################
    ## calculate weekly percentage of all samples with significant, positive excess freq per mutation
    if(opt$verbose){
     writeLines(paste("PROGRESS: calculate weekly percentage of postive samples per mutation"))
    }

    dt_geo1 <- currentmetaDTs %>% group_by(midweek_date) %>% summarize(N = n(), .groups = "keep")
    dt_geo2 <- currentAFdata_completed %>% filter(excess > 0) %>% group_by(midweek_date, nuc_mutation) %>% summarize(n = n(), .groups = "keep")
    dt_geo <- full_join(x = dt_geo1, y = dt_geo2, by = "midweek_date") %>% rowwise() %>% mutate(r = n / N)
    dt_geo <- data.table::dcast(dt_geo, nuc_mutation~midweek_date, fun.aggregate = mean, value.var = "r", fill = 0)  %>% data.table::melt(id.vars = c("nuc_mutation"), variable.name = "midweek_date", value.name = "r")


    ####################
    ## deduce nls growth of timecourse geographic spread per mutation
    ## filter  mutations with at least minpossamp_th timepoint with a r>0 and an max(r) > r_th
    if(opt$verbose){
      writeLines(paste("PROGRESS: deduce nls growth of timecourse geographic spread per mutation"))
    }

    dt_geo$midweek_date <- as.Date(dt_geo$midweek_date)
    ind.mois <- dt_geo %>% group_by(nuc_mutation) %>% mutate(max = max(r)) %>% mutate(n = sum(r > 0))  %>% rowwise() %>% filter(max >= r_th & n >= minpossamp_th) %>% pull(nuc_mutation) %>% unique()

    prediction_collector_geo_List <- list()
    prediction_collector_geo <- data.table(mutation = character(), growth_pred = numeric(), inflection_pred = numeric(), tdata = character(), odata = character(), pdata = character())

    for (ind.moi in ind.mois){
      dtoi <- dt_geo %>% filter(nuc_mutation == ind.moi)
      colnames(dtoi)[2] <- "t"
      colnames(dtoi)[3] <- "y"
      subset(dtoi, select = c("t", "y")) -> dtoi
      dtoi$d <- as.numeric(as.Date(dtoi$t)-as.Date(dtoi$t[1]))
      a <- 1
      fo <- y ~ a / (1 + exp(-b * (d-c)))
      model <- tryCatch(gsl_nls(fo, start = list(b = 0.05, c = 2*max(dtoi$d)/3), data = dtoi),error=function(e) NA, warning=function(w) NA)
      # if gsl_nls regession did not converge try different starting parameters
      if(1){
          if(identical(NA, model)){
            model <- tryCatch(gsl_nls(fo, start = list(b = 0.01, c = max(dtoi$d)/2), data = dtoi),error=function(e) NA, warning=function(w) NA)
          }
          if(identical(NA, model)){
            model <- tryCatch(gsl_nls(fo, start = list(b = 0.1, c = 2*max(dtoi$d)/3), data = dtoi),error=function(e) NA, warning=function(w) NA)
          }
          if(identical(NA, model)){
            model <- tryCatch(gsl_nls(fo, start = list(b = 0.2, c = max(dtoi$d)/2), data = dtoi),error=function(e) NA, warning=function(w) NA)
          }
      }

      if(identical(NA, model)){
         writeLines(paste("WARNING: no convergence for", ind.moi, "assessing geographic spead"))
      } else{
         # predict the relative weekly growth in the last 14days
         tdays <- 14
         dtoi$predict <- predict(model)
         lastbiweek_prediction <- predict(model, newdata = data.frame(t = c(as.Date(periodend)-tdays, periodend), d = as.numeric(c(as.Date(periodend)-tdays, periodend)-as.Date(dtoi$t[1]))))
         lastweek_growth <- (lastbiweek_prediction[2] - lastbiweek_prediction[1])/(tdays/7)
         lastweek_growth -> growth_pred
         coefficients(model)[2] -> inflection_pred

         ## set inflection_pred and growth_pred to NA growth_pred confidence interval overlaps 0
         summary(model) -> model_summary
         model_summary$coefficients[ , 4][1] -> pvalue_growth_pred
         model_summary$coefficients[ , 4][2] -> pvalue_inflection_pred
         CI_growth_pred <- tryCatch(suppressMessages(confintd(model, "b", level = 0.90)),error=function(e) NA, warning=function(w) NA)
         CI_inflection_pred <- tryCatch(suppressMessages(confintd(model, "c", level = 0.90)),error=function(e) NA, warning=function(w) NA)

         if (any(is.na(CI_growth_pred)) | (min(CI_growth_pred[2:3]) < 0 & max(CI_growth_pred[2:3]) > 0) | growth_pred < opt$geogrowthlimit){
           growth_pred = NA
           inflection_pred = NA
         } else{
            if (!is.na(growth_pred) & opt$verbose == TRUE){
              melt(dtoi[,-3], id.vars = "t") %>% mutate(t = as.Date(t)) %>%
                  ggplot(aes(x = t, y = value, color = variable)) +
                  geom_point(size = 3, alpha = 0.66, aes(shape = variable, size = variable)) +
                  geom_smooth(formula = y ~ x, method = "glm", method.args = list(family = "binomial"), se = FALSE) +
                  theme_bw() +
                  ggtitle(paste("Till", periodend, nuc2label$label[nuc2label$NUC == ind.moi], sep = ": ")) +
                  theme(legend.position="bottom", legend.title = element_blank()) +
                  ylab("mutation frequency [1/1]") +
                  xlab("Sample Date") +
                  #geom_vline(xintercept=as.Date(dtoi$t[1])-dtoi$d[1]+as.numeric(inflection_pred), linetype = "dotted")  +
                  scale_color_brewer(palette = "Set2", breaks = c("y", "predict"), labels = c("Observations", "Model Prediction")) +
                  scale_shape(breaks = c("y", "predict"), labels = c("Observations", "Model Prediction")) +
                  annotate(x = median(dtoi$t), y = max(dtoi$predict, dtoi$y), geom = "text", label = paste("Growth per week:\n", signif(growth_pred, digits = 2))) -> p
              ggsave(filename = paste0(outdir, "/figs/singleModels/", "plot_geospread_", periodend, "_", ind.moi, ".", opt$graph), p, width = 5, height = 5)

            }
        }

        prediction_geo <- data.table(mutation = ind.moi, growth_pred = growth_pred, inflection_pred = inflection_pred, tdata = paste(dtoi$t, collapse="|"), odata = paste(dtoi$y, collapse="|"), pdata = paste(dtoi$predict, collapse="|"))
        prediction_collector_geo_List[[length(prediction_collector_geo_List)+1]] <- prediction_geo
      }
      rm(dtoi, model)
    }
    prediction_collector_geo <- rbindlist(prediction_collector_geo_List)
    geo_growing_nucs <- prediction_collector_geo %>% filter(!is.na(growth_pred))

    ####################
    ## join mutations growing geographically and/or by weighted mean excess freq.
    if(opt$verbose){
      writeLines(paste("PROGRESS: join mutations growing geographically", opt$filstrat, "by weighted mean excess freq"))
    }

    if(is.null(geo_growing_nucs$mutation) | is.null(wmef_growing_nucs$mutation)){
      next;
    }
    growing_nucs_or  <- full_join(x = geo_growing_nucs, y = wmef_growing_nucs, by = "mutation", suffix = c(".geo", ".wmef"))
    growing_nucs_and <- inner_join(x = geo_growing_nucs, y = wmef_growing_nucs, by = "mutation", suffix = c(".geo", ".wmef"))

    dt <- full_join(x = dt_geo, y = dt_wmef, by = c("nuc_mutation", "midweek_date"))
    dt <- left_join(x = dt, y = nuc2label, by = c("nuc_mutation" = "NUC"))
    dt_filtered_or <- dt %>% filter(nuc_mutation %in% growing_nucs_or$mutation)
    dt_filtered_and <- dt %>% filter(nuc_mutation %in% growing_nucs_and$mutation )

    ## set here if final filtered should AND or OR linked geo wmef results
    if(opt$filstrat == "or"){
      dt_filtered <- dt_filtered_or
      growing_nucs <- growing_nucs_or
    } else{
      dt_filtered <- dt_filtered_and
      growing_nucs <- growing_nucs_and
    }

    if(dim(dt_filtered)[1]==0){
      writeLines(paste("LOG: no mutations found to grow in frequency and geographically at", periodend, "; no results to be written in output tables and plots; timepoint skipped!"))
      next;
    }
    currentMutationSet <- dt_filtered %>% mutate(latest = max(midweek_date)) %>% filter(midweek_date == latest)  %>% left_join(y = growing_nucs, by = c("nuc_mutation" = "mutation")) %>% dplyr::select(nuc_mutation, label, midweek_date, r, excess_mw, growth_pred.geo, growth_pred.wmef)
    timepointMutationSet <- currentMutationSet %>% mutate(midweek_date = as.character(midweek_date))
    timecourseMutationSet_List[[length(timecourseMutationSet_List)+1]] <- timepointMutationSet

    #####################
    ## produce kinetic plot of filtered mutations
    ## filter mutations based on geographical spead and/or by weighted mean excess freq. in the latest time point
    if(dim(dt_filtered)[1]>0){
      dt_filtered <- dt_filtered %>% mutate(midweek_date = as.Date(midweek_date)) %>% group_by(nuc_mutation, label) %>% mutate(latest_excess_mw = excess_mw[which(max(midweek_date) == midweek_date)]) %>% mutate(latest_r = r[which(max(midweek_date) == midweek_date)])
      prediction_data <- cbind(geo_growing_nucs %>% filter(mutation %in% dt_filtered$nuc_mutation) %>% mutate(pdata = strsplit(as.character(pdata), "\\|")) %>% unnest(pdata) %>% dplyr::select(mutation, pdata), geo_growing_nucs %>% filter(mutation %in% dt_filtered$nuc_mutation) %>% mutate(tdata = strsplit(as.character(tdata), "\\|")) %>% unnest(tdata) %>% dplyr::select(tdata))
      prediction_data <- left_join(x = prediction_data, y = nuc2label, by = c("mutation" = "NUC"))
    }

    if(dim(dt_filtered)[1]==0){
      writeLines(paste("LOG: no mutations found to grow in frequency and geographically at", periodend, "; no results to be written in output tables and plots; timepoint skipped!"))
      next;
    }

    if( opt$TimecoursePlot == "all" | (opt$TimecoursePlot == "last" & abs(periodend - max(globalAFdata_m$midweek_date)) < 7) ){

        if(dim(dt_filtered)[1]>0){
          if(opt$verbose){
            writeLines(paste("PROGRESS: produce kinetic plot of filtered mutations"))
          }

          kineticPlot <- ggplot(dt_filtered, aes(x = midweek_date, y = r, group = label)) +
              geom_smooth(data = prediction_data, aes(x = as.Date(tdata), y = as.numeric(pdata)), formula = y ~ x, method = "glm", method.args = list(family = "binomial"), se = FALSE) +
              geom_point(alpha = .9, aes(size = excess_mw, color = excess_mw)) +
              theme_bw() +
              scale_colour_viridis_b(begin = 0.1, end = 0.9, option = "D", direction = +1, name = "Weighted mean\nexcess frequency") +
              scale_size(name = "Weighted mean\nexcess frequency", range = c(2, 8)) +
              facet_wrap(~label) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
              ylab("fraction of positive wwtp [1/1]") +
              xlab("") +
              theme(legend.position="bottom")

            filename <- paste0(outdir, "/figs/",  paste('/kineticPlot', as.character(latestSample$latest), sep="_"), ".", opt$graph)
            plot.width <- 5 + 1.2*ceiling(sqrt(length(unique(dt_filtered$nuc_mutation))))
            plot.height <- 5 + ceiling(length(unique(dt_filtered$nuc_mutation))/ceiling(sqrt(length(unique(dt_filtered$nuc_mutation)))))
            plot.width <- ifelse(plot.width > 49, 49, plot.width)
            plot.height <- ifelse(plot.height > 49, 49, plot.height)
            ggsave(filename = filename, plot = kineticPlot, width = plot.width, height = plot.height)
            if( abs(periodend - max(globalAFdata_m$midweek_date)) < 7 ){
              fwrite(as.list(c("excessmutations", "kineticPlot", as.character(periodend), filename)), file = summaryDataFile, append = TRUE, sep = "\t")
            }
        }

        #####################
        ## produce outbreak.info like heatmap of filtered mutations
        ## generate outbreak.info like heatmap for excess mutations which are also sig. growing over time
        if(dim(dt_filtered)[1]>0){
          if(opt$verbose){
            writeLines(paste("PROGRESS: generate outbreak.info like heatmap for growing excess mutations"))
          }

          unique(dt_filtered$nuc_mutation) -> growing_mutations
          globalAFdata %>% group_by(nuc_mutation) %>% mutate(latestTP = max(sample_date)) %>% rowwise() %>% filter(latestTP == sample_date) %>% filter(nuc_mutation %in% growing_mutations) %>% group_by(nuc_mutation) %>% summarize(aa_mutation = aa_mutation, expected = expected, observed = observed, excess = excess, pvalue = pvalue, .groups = "keep") %>% left_join(y = nuc2label, by = c("nuc_mutation" = "NUC")) %>% arrange(desc(excess)) -> outbreak.selection
          outbreak.selection <- outbreak.selection %>% group_by(nuc_mutation) %>% summarize(aa_mutation = unique(aa_mutation), expected = mean(expected), observed = mean(observed), excess = mean(excess), pvalue = mean(pvalue), label = unique(label), .groups = "keep")
          data.table::melt(outbreak.selection, id.vars = c("label"), measure.vars = c("expected", "observed")) -> outbreak.freqs
          mstat %>% filter(nucc %in% outbreak.selection$nuc_mutation) -> outbreak.dt
          if(dim(outbreak.dt)[1]==0){
            next;
          }
          left_join(x = outbreak.dt, y = nuc2label, by = c("nucc" = "NUC")) %>% data.table::dcast(formula = ID ~ label, value.var = "sensitivity", fill = 0) -> outbreak.dt
          data.table::melt(outbreak.dt, variable.name = "label", value.name = "AF.freq") -> outbreak.dt
          outbreak.dt %>% group_by(ID) %>% mutate(max = max(AF.freq, na.rm = TRUE)) %>% filter(max >= min(max(outbreak.dt$AF.freq), .9)) %>% dplyr::select(-"max") -> outbreak.dt
          outbreak.dt %>% mutate(dealiasID = dealias(ID)) -> outbreak.dt
          # remove variants which are represented by a ancestor with same fingerwriteLines
          outbreak.dt %>% group_by(ID) %>%arrange(label) %>% mutate(fingerwriteLines = paste(ifelse(AF.freq > 0.5, 1, 0), collapse = "")) %>% group_by(fingerwriteLines) %>% mutate(groupIds = paste(unique(dealiasID), sep = ";", collapse = ";")) %>% rowwise() %>% mutate(keep = collapse2mother(dealiasID, groupIds)) %>% filter(keep) -> outbreak.dt


          outbreakInfoPlot <- ggplot(data = outbreak.dt, aes(x = label, y = ID, fill = AF.freq))
          outbreakInfoPlot <- outbreakInfoPlot + geom_tile(color = "white", alpha = 0.8)
          if(dim(outbreak.dt)[1] <= 1){
            outbreakInfoPlot <- outbreakInfoPlot + scale_fill_viridis_c(name = "Allele frequency", trans = "sqrt", option = "B", begin = 0.1, end = 0.9, guide = guide_colorbar(direction = "horizontal", title.position = "top", label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, label.theme = element_text(angle = 90)))
          } else if(length(unique(outbreak.dt$AF.freq)) > 1){
            outbreakInfoPlot <- outbreakInfoPlot + scale_fill_viridis_b(name = "Allele frequency", trans = "sqrt", option = "B", begin = 0.1, end = 0.9, guide = guide_colorbar(direction = "horizontal", title.position = "top", label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, label.theme = element_text(angle = 90)))
          } else{
            outbreakInfoPlot <- ggplot(data = outbreak.dt, aes(x = label, y = ID, fill = AF.freq))
            outbreakInfoPlot <- outbreakInfoPlot + geom_tile(color = "white", alpha = 0.8)
            outbreakInfoPlot <- outbreakInfoPlot + scale_fill_viridis_c(name = "Allele frequency", option = "B", begin = 0.1, end = 0.9, guide = guide_colorbar(direction = "horizontal", title.position = "top", label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, label.theme = element_text(angle = 90)))
          }

          outbreakInfoPlot <- outbreakInfoPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
          outbreakInfoPlot <- outbreakInfoPlot + xlab("Excess Mutations") + ylab("Variants")
          outbreakInfoPlot <- outbreakInfoPlot + coord_flip()
          outbreakInfoPlot <- outbreakInfoPlot + theme(legend.position="bottom", legend.direction="horizontal")
          outbreakInfoPlot <- outbreakInfoPlot + ggtitle("Mutation associated variants")
          filename <- paste0(outdir, "/figs/",  paste('/outbreakInfoPlot', as.character(latestSample$latest), sep="_"), ".", opt$graph)
          plot.width <- 5 + 1.2*ceiling(sqrt(length(unique(dt_filtered$nuc_mutation))))
          plot.height <- 5 + ceiling(length(unique(dt_filtered$nuc_mutation))/ceiling(sqrt(length(unique(dt_filtered$nuc_mutation)))))
          plot.width <- ifelse(plot.width > 49, 49, plot.width)
          plot.height <- ifelse(plot.height > 49, 49, plot.height)
          ggsave(filename = filename, plot = outbreakInfoPlot, width = plot.width, height = plot.height)
          if( abs(periodend - max(globalAFdata_m$midweek_date)) < 7 ){
            fwrite(as.list(c("excessmutations", "outbreakInfoPlot", as.character(periodend), filename)), file = summaryDataFile, append = TRUE, sep = "\t")
          }

        }

        #####################
        ## make table of all excess growing mutations
        if(dim(dt_filtered)[1]>0 & exists("outbreak.selection")){

          if(opt$verbose){
            writeLines(paste("PROGRESS: make table of all excess growing mutations"))
          }

          legendTxt <- paste0("Überschuss-Mutationen die in der Woche um ", periodend, " österreichweit sig. Wachstum und sig. Ausbreitung zeigen")
          filename <- paste0(outdir, "/figs/",  paste('/excessmutationTable', periodend, sep="_"), ".tex")
          nuccToUse <- unique(dt_filtered$nuc_mutation)
          left_join(x = outbreak.selection, y = (growing_nucs %>% dplyr::select(-tdata, -pdata)), by = c("nuc_mutation" = "mutation"))  %>% filter(nuc_mutation %in% nuccToUse) %>% dplyr::select(nuc_mutation, aa_mutation, expected, observed, growth_pred.geo, growth_pred.wmef) -> growing_excessmutationsTable.dt
          growing_excessmutationsTable.dt %>% rowwise() %>% mutate('cov.link' = covspectrumLinkSimple(nuc_mutation)) -> growing_excessmutationsTable.dt
          growing_excessmutationsTable.dt %>% mutate(expected = ifelse(expected < zeroo/100, 0, expected)) -> growing_excessmutationsTable.dt
          growing_excessmutationsTable.dt %>% mutate(observed = ifelse(observed < zeroo/100, 0, observed)) -> growing_excessmutationsTable.dt
          growing_excessmutationsTable.dt %>% mutate(expected = signif(expected, digits = 2), observed = signif(observed, digits = 2), growth_pred.geo = signif(growth_pred.geo, digits = 2), growth_pred.wmef = signif(growth_pred.wmef, digits = 2)) -> growing_excessmutationsTable.dt
          rename_lookup <- c("Nuc Mutation" = "nuc_mutation", "AA Mutation" = "aa_mutation", "Erw. AF" = "expected", "Beob. AF" = "observed", "Ausbreitung [1/w]" = "growth_pred.geo", "Wachstum [1/w]" = "growth_pred.wmef", "cov-spectrum" = "cov.link")
          rename(growing_excessmutationsTable.dt, all_of(rename_lookup)) -> growing_excessmutationsTable.dt
          makeTexTab(filename, growing_excessmutationsTable.dt, legendTxt)
          if( abs(periodend - max(globalAFdata_m$midweek_date)) < 7 ){
            fwrite(as.list(c("excessmutations", "table", as.character(periodend), filename)), file = summaryDataFile, append = TRUE, sep = "\t")
          }
        }

        #####################
        ## compile and writeLines result table
        if(dim(dt_filtered)[1]>0 & exists("wmef_growing_nucs") & exists("geo_growing_nucs") & exists("outbreak.dt")){
          if(opt$verbose){
            writeLines(paste("PROGRESS: compile and writeLines result table"))
          }
          results_compiled <- dt_filtered %>% group_by(nuc_mutation) %>% mutate(latest = max(midweek_date)) %>% filter(latest == midweek_date) %>% dplyr::select(nuc_mutation, label, midweek_date, latest_excess_mw, latest_r)
          results_compiled <- left_join(x = results_compiled, y = (wmef_growing_nucs %>% dplyr::select(mutation, growth_pred)), by = c("nuc_mutation" = "mutation")) %>% rename(wmef_growth = growth_pred)
          results_compiled <- left_join(x = results_compiled, y = (geo_growing_nucs %>% dplyr::select(mutation, growth_pred)), by = c("nuc_mutation" = "mutation")) %>% rename(geo_growth = growth_pred)
          positiveFractionFromString <- function(string_){
            numbs <- as.numeric(unlist(strsplit(string_, split="\\|")))
            return(sum(numbs>zeroo)/length(numbs))
          }
          results_compiled <- left_join(x = results_compiled, y = (geo_growing_nucs %>% rowwise() %>% mutate(positiveFractionFromString = positiveFractionFromString(odata)) %>% dplyr::select(mutation, positiveFractionFromString)), by = c("nuc_mutation" = "mutation"))
          results_compiled <- left_join(x = results_compiled, y = (outbreak.dt %>% rowwise() %>% filter(AF.freq > 0.8) %>% group_by(label) %>% summarize(Variants = paste(ID, sep = ",", collapse = ", "), .groups = "keep")), by = c("label"))
          filename <- paste0(outdir,  paste('/globalExcess', sep="_"), ".csv")
          fwrite(results_compiled, file = filename, sep = "\t")
          if( abs(periodend - max(globalAFdata_m$midweek_date)) < 7 ){
            fwrite(as.list(c("results", "table", as.character(periodend), filename)), file = summaryDataFile, append = TRUE, sep = "\t")
          }
        }
    }
}
timecourseMutationSet <- rbindlist(timecourseMutationSet_List )
timecourseMutationSet$midweek_date <- as.Date(timecourseMutationSet$midweek_date)


#####################################################
##### make plots as done in cuadrilla with data in timecourseMutationSet
## filter mutations which happen to grow in at least opt$nw consecutive weeks

writeLines(paste("LOG: filter mutations which happen to grow in at least", opt$nw, "consecutive weeks"))

span <- seq(from = as.Date(opt$start), by = 7, to = as.Date(opt$end))
span <- data.table(midweek_date = span + 3 - as.numeric(strftime(span, format = "%u")))
span$midweek_date <- as.Date(span$midweek_date)


### kill if no mutations are found
if (!exists("timecourseMutationSet")){
  if(length(unique(timecourseMutationSet$nuc_mutation)) < 1){
    writeLines(paste("WARNING: no mutation data found in <timecourseMutationSet>. Program ends here!"))
    quit(save="no")
  }
}

consecutive_mutation_to_use <- timecourseMutationSet %>%
    dplyr::select(nuc_mutation, label, midweek_date) %>%
    distinct() %>%
    group_by(label, nuc_mutation) %>%
    arrange(midweek_date, .by_group = TRUE) %>%
    mutate(prior = as.numeric(midweek_date - lag(midweek_date, n = opt$nw-1))) %>%
    dplyr::select(nuc_mutation, label, midweek_date, prior) %>%
    mutate(maxprio = min(prior, na.rm = TRUE)) %>%
    filter(maxprio >= (opt$nw-1)*7 & maxprio < (opt$nw)*7) %>%
    pull(nuc_mutation) %>%
    unique()
timecourseMutationSet <- timecourseMutationSet %>% filter(nuc_mutation %in% consecutive_mutation_to_use)

### kill if no mutations are found
if (length(unique(consecutive_mutation_to_use)) < 1){
  writeLines(paste("WARNING: no mutation data found after filtering for", opt$nw, "consecutive mutations. Program ends here!"))
  quit(save="no")
}


####  count mutations per week an detect outlier
if(opt$verbose){
  writeLines(paste("PROGRESS: count mutations per week an detect outlier"))
}
mut_per_week <- left_join(x = span, y = timecourseMutationSet, by = c("midweek_date")) %>%
    group_by(midweek_date) %>%
    summarize(n = sum(!is.na(label)), .groups = "keep")
mean_mut_per_week <- mean(mut_per_week$n)
mut_per_week <- mut_per_week %>% rowwise() %>% mutate(p = ppois(n, lambda = mean_mut_per_week, lower.tail = FALSE))
mut_per_week$q <- p.adjust(mut_per_week$p)


## produce outbreak.info like heatmap of filtered mutations
## generate outbreak.info like heatmap for plotted excess mutations
## read mutation stats
if(opt$verbose){
  writeLines(paste("PROGRESS: generate outbreak.info like heatmap for plotted excess mutations"))
}

mstatf <- mstat %>% filter(nucc %in% consecutive_mutation_to_use)
mstatf <- mstatf %>% mutate(label = paste0(nucc, " [", AA_change, "]"))
mstatf <- mstatf %>% mutate(label = toupper(label))

## identify mutations which are seen in timecourseMutationSet but not in mstatf
timecourseMutationSet %>% mutate(label = toupper(label)) -> timecourseMutationSet
considered_mutations <- timecourseMutationSet %>% dplyr::select(nuc_mutation,  label) %>% distinct()
missed_in_mstatf <- considered_mutations %>% filter(!nuc_mutation %in% mstatf$nucc)

## dcast-complete and timecourseMutationSet-complete mstatf
data.table::melt(mstatf, measure.vars = c("sensitivity"), value.name = "AF.freq") -> outbreak.dt
outbreak.dt %>% group_by(ID) %>% mutate(max = max(AF.freq, na.rm = TRUE)) %>% filter(max >= min(max(outbreak.dt$AF.freq), 0)) %>% dplyr::select(-"max") -> outbreak.dt
outbreak.dt %>% data.table::dcast(formula = ID ~ label, value.var = "AF.freq", fill = 0) -> outbreak.dt
if(dim(missed_in_mstatf)[1] > 0){
  for ( i in missed_in_mstatf$label){
    if(!i %in% colnames(outbreak.dt)){
      outbreak.dt[,(dim(outbreak.dt)[2]+1)] <- rep(0, dim(outbreak.dt)[1])
      colnames(outbreak.dt)[dim(outbreak.dt)[2]] <- i
    }
  }
}

## define similarity between mutation in mstatf
if(is.data.frame(outbreak.dt[,1])){
  tree <- hclust(dist(t(outbreak.dt[,-1])))
  mutation_order <- tree$labels[tree$order]
  tree <- hclust(dist(data.frame(outbreak.dt[,-1], row.names = outbreak.dt[,1])))
  lineage_order <- tree$labels[tree$order]
} else{
  mutation_order <- colnames(outbreak.dt)[-1]
  lineage_order <- outbreak.dt$ID
}
data.table::melt(outbreak.dt, variable.name = "label", value.name = "AF.freq", id.vars = c("ID")) -> outbreak.dt
outbreak.dt %>% rowwise() %>% mutate(dealiasID = dealias(ID)) -> outbreak.dt

outbreak.dt %>% group_by(ID) %>%
    arrange(label) %>%
    mutate(fingerprint = paste(ifelse(AF.freq > 0.5, 1, 0), collapse = "")) %>%
    group_by(fingerprint) %>%
    mutate(groupIds = paste(unique(dealiasID), sep = ";", collapse = ";")) %>%
    rowwise() %>% mutate(keep = collapse2mother(dealiasID, groupIds)) %>% filter(keep) -> outbreak.dt

outbreak.dt$label <- factor(outbreak.dt$label, levels = mutation_order)
outbreak.dt$ID <- factor(outbreak.dt$ID, levels = lineage_order)

outbreakInfoPlot <- ggplot(data = outbreak.dt, aes(x = label, y = ID, fill = AF.freq))
outbreakInfoPlot <- outbreakInfoPlot + geom_tile(color = "white", alpha = 0.8)
if(dim(outbreak.dt)[1] <= 1){
  outbreakInfoPlot <- outbreakInfoPlot + scale_fill_viridis_b(name = "Allele frequency", trans = "sqrt", option = "B", begin = 0.1, end = 0.9)
} else{
  outbreakInfoPlot <- outbreakInfoPlot + scale_fill_viridis_c(name = "Allele frequency", trans = "sqrt", option = "B", begin = 0.1, end = 0.9)
}
outbreakInfoPlot <- outbreakInfoPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
outbreakInfoPlot <- outbreakInfoPlot + xlab("") + ylab("Variants")
outbreakInfoPlot <- outbreakInfoPlot + coord_flip()
outbreakInfoPlot <- outbreakInfoPlot + theme_minimal()
outbreakInfoPlot <- outbreakInfoPlot + theme(panel.grid = element_blank())
outbreakInfoPlot <- outbreakInfoPlot + theme(legend.position="right")
outbreakInfoPlot <- outbreakInfoPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
outbreakInfoPlot <- outbreakInfoPlot + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
outbreakInfoPlot <- outbreakInfoPlot + theme(axis.text.y = element_text(hjust = 0.5))
outbreakInfoPlot <- outbreakInfoPlot + theme(axis.title=element_blank())

plot2 <- plot_spacer()/outbreakInfoPlot + plot_layout(heights = c(1,4),  guides = 'collect') & theme(legend.position = "right")

### make plot mutation vs week
if(opt$verbose){
  writeLines(paste("PROGRESS: generate plot mutation vs week"))
}

timecourseMutationSet$label <- factor(timecourseMutationSet$label, levels = mutation_order)
p <- ggplot(data = timecourseMutationSet, aes(x = midweek_date, y = label)) +
    geom_point(aes(color = growth_pred.geo), size = 6, alpha  =.6) +
    scale_x_date(date_breaks = "4 weeks", limits = c(as.Date(opt$start), as.Date(opt$end))) +
    theme_bw() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    ylab("Mutations") +
    xlab("") +
    scale_colour_viridis_b(begin = .1, end = .9, direction = -1, option = "D", name = "weekly growth\ngeographic\nspead")

if(length(levels(timecourseMutationSet$label)) == length(levels(timecourseMutationSet$label))){
  if(all(levels(timecourseMutationSet$label) == levels(timecourseMutationSet$label))){
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
}

q <- ggplot(data = mut_per_week, aes(x = midweek_date, y = n)) +
    geom_col(aes(fill = ifelse(q<pval_th, "sig.", "n.s."))) +
    theme_minimal() +
    scale_x_date(date_breaks = "2 weeks", limits = c(as.Date(opt$start), as.Date(opt$end))) +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    ylab("Number of\nmutations") +
    xlab("") +
    scale_fill_brewer(palette = "Pastel2", name = "mutation\nenrichment")

plot1 <- q/p + plot_layout(heights = c(1,3),  guides = 'collect') & theme(legend.position = "left")


## save plot
design <- "
13
24
24
24
"
plot <- q + p  + plot_spacer() + outbreakInfoPlot + plot_layout(design=design,  guides = 'collect') & theme(legend.position = "left")

filename <- paste0(opt$dir, "/",  paste('figs/cuadrilla', sep="_"), ".", opt$graph)
plot.width <- 5 + 0.2*((length(unique(mut_per_week$midweek_date))/4) + length(unique(outbreak.dt$ID)))
plot.height <- 3 + 0.33*length(consecutive_mutation_to_use)
plot.width <- plot.width * 0.87
plot.height <- plot.height * 0.77
plot.width <- ifelse(plot.width > 49, 49, plot.width)
plot.height <- ifelse(plot.height > 49, 49, plot.height)
ggsave(filename = filename, plot = plot, width = plot.width, height = plot.height)
fwrite(as.list(c("cuadrilla", "plot", as.character(max(timecourseMutationSet$midweek_date)), filename)), file = summaryDataFile, append = TRUE, sep = "\t")


##########################
## print all enriched mutation

if(opt$verbose){
  writeLines(paste("PROGRESS: print all enriched mutation to file"))
}

rdt <- left_join(x = timecourseMutationSet, y = mut_per_week, by = "midweek_date") %>%
    filter(p <= pval_th) %>%
    dplyr::select(nucMutation = nuc_mutation, label, midweekDate = midweek_date, excessFreq=excess_mw, spread=r, growth_excessFreq=growth_pred.wmef, growth_spread=growth_pred.geo, mutation_enrichment_pval=q)
filename <- paste0(opt$dir, "/",  paste('cuadrilla', sep="_"), ".csv")
fwrite(x = rdt, file = filename)

##### make plots as done in NucList_surveillance.r for mutations in last months data

if(opt$verbose){
  writeLines(paste("PROGRESS: generate time course plots for excessive growing mutations"))
}

moi.list <- rdt %>% dplyr::select(nucMutation, label) %>% distinct()

moi.label.list <- as.vector(moi.list$label)
moi.nuc.list <- as.vector(moi.list$nuc)

mutation.dt <- sewage_samps.dt %>% filter(NUC %in% moi.nuc.list)

left_join(x = mutation.dt, y = metaDT, by = c("ID" = "BSF_sample_name", "RNA_ID_int", "LocationID", "LocationName", "sample_date", "sample_date_decimal")) -> mutation.dtm
mutation.dtm %>% filter(!is.na(sample_date)) -> mutation.dtm

mutation.dtm %>% group_by(LocationID, LocationName, state, sample_date) %>% summarize(n = sum(!is.na(unique(NUC))), NUCS = paste(unique(NUC), sep = " & ", collapse = " & "),  AFS = paste(unique(value.freq), sep = " & ", collapse = " & "), .groups = "keep") -> mutation.dtm


if(opt$verbose){
  writeLines(paste("PROGRESS: generate co-occurence plot for excessive growing mutations"))
}

if(dim(mutation.dtm)[1] > 1){
  cooccurencePlot <- ggplot(data = mutation.dtm, aes(x = sample_date, y = LocationID, fill = n)) +
          geom_point(shape = 21, aes(size = n), color = "grey75") +
          geom_text(aes(label = ifelse(n>0, n, ""))) +
          theme_bw() + facet_grid(state~., scale = "free_y", space = "free_y") +
          theme(axis.text.y = element_blank(), strip.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank()) +
          scale_x_date(date_breaks="1 months", date_labels = "%b %Y") +
          ylab("") + xlab("Sample Date") +
          ggtitle(paste("Marker co-occurrence, out of", length(moi.nuc.list), "mutations considered")) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_line(color = "grey75"))

  if(length(unique(mutation.dtm$NUCS)) > 1){
    cooccurencePlot <- cooccurencePlot + scale_fill_fermenter(palette = "Oranges", direction = 1)
  } else{
    cooccurencePlot <- cooccurencePlot + scale_fill_distiller(palette = "Oranges", direction = 1)
  }

  filename <- paste0(outdir, "/figs/nucleotides",  paste('/Nuc', "cooccurence", sep="_"), ".", opt$graph)
  ggsave(filename = filename, plot = cooccurencePlot, height = 9, width = 9)
}

if(length(moi.nuc.list) > 0){
  for (moi.idx in 1:length(moi.nuc.list)){
        moi.nuc <- moi.nuc.list[moi.idx]
        moi.label <- moi.label.list[moi.idx]
        if(opt$verbose){
          writeLines(paste("PROGRESS: generate time course plots for excessive growing mutations:", moi.nuc))
        }


        ######################################################
        # plot AF per wwtp sample of mutation of interest
        mutation.dt %>% filter(NUC %in% moi.nuc) -> mutation.dtm


        full_join(x = mutation.dtm, y = metaDT, by = c("ID" = "BSF_sample_name", "RNA_ID_int", "LocationID", "LocationName", "sample_date", "sample_date_decimal")) %>% rename(AF = value.freq) -> mutation.dtm
        mutation.dtm %>% filter(!is.na(sample_date)) -> mutation.dtm

        p <- ggplot(data = mutation.dtm, aes(x = sample_date, y = LocationID, fill = AF)) +
            geom_point(shape = 21, aes(size = AF), color = "grey75") +
            theme_bw() +
            facet_grid(state~., scale = "free_y", space = "free_y") +
            theme(axis.text.y = element_blank(), strip.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank()) +
            scale_x_date(date_breaks="1 months", date_labels = "%b %Y") +
            ylab("") +
            xlab("Sample Date") +
            ggtitle(paste(moi.label)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_line(color = "grey75"))


        ######################################################
        # plot number of moi positive samples per two weeks (apply sliding window: week+previous_week)
        mutation.dt %>% filter(NUC %in% moi.nuc) %>% pull(ID) -> positive_samples.dt
        metaDT %>% mutate(mutpos = ifelse(BSF_sample_name %in% positive_samples.dt, 1, 0)) -> dt1
        dt1 %>% mutate(midweek_date = as.Date(sample_date, tryFormats = c("%Y-%m-%d")) + 3 - as.numeric(strftime(as.Date(sample_date, tryFormats = c("%Y-%m-%d")), format = "%u"))) -> dt1
        dt1 %>% arrange() %>% group_by(midweek_date) %>% summarize(n = n(), p = sum(mutpos==1)) %>% mutate(nw = n+lag(n), pw = p+lag(p)) %>% mutate(r = p/n, rw = pw/nw) %>% filter(!is.na(r)) -> dt1
        dt1$state <- "Combined"

        q1 <- ggplot(data = dt1, aes(x = midweek_date, y = pw)) +
            geom_point() + geom_step(direction = "mid") +
            scale_x_date(date_breaks="2 months", date_labels = "%b %Y") +
            ylab("") +
            xlab("Sample Date") +
            theme_bw() +
            ggtitle(paste(moi.nuc, "positive samples [1]")) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5), panel.grid.major.x = element_line(color = "grey75")) +
            facet_grid(state~., scale = "free_y", space = "free_y")
        q2 <- ggplot(data = dt1, aes(x = midweek_date, y = rw)) +
            geom_step(direction = "mid") +
            geom_point(shape = 21, size = 3.6, fill = "white") +
            geom_text(aes(label = pw), color = "steelblue", nudge_y = 0, size = 2) +
            scale_x_date(date_breaks="2 months", date_labels = "%b %Y") +
            ylab("") +
            xlab("Sample Date") +
            theme_bw() +
            ggtitle(paste(moi.nuc, "positive samples fraction [1/1]")) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5), panel.grid.major.x = element_line(color = "grey75")) +
            facet_grid(state~., scale = "free_y", space = "free_y")

        #####################################################
        # plot populated weighted mean mutation AF per  week
        left_join(x = metaDT, y = mutation.dt %>% filter(NUC %in% moi.nuc), by = c("BSF_sample_name" = "ID", "LocationID", "LocationName", "RNA_ID_int", "sample_date", "sample_date_decimal")) %>%
            mutate(midweek_date = as.Date(sample_date, tryFormats = c("%Y-%m-%d")) + 3 - as.numeric(strftime(as.Date(sample_date, tryFormats = c("%Y-%m-%d")), format = "%u"))) %>%
            mutate(value.freq = ifelse(is.na(value.freq), 0, value.freq)) %>%
            group_by(midweek_date) %>%
            summarize(mw = weighted.mean(value.freq, w = connected_people)) -> dt2
        dt2$state <- "Combined"
        dt2 %>% rename(AF = mw, sample_date = midweek_date) -> dt2
        r <- ggplot(data = dt2, aes(x = sample_date, y = AF, fill = AF)) +
            geom_step(color = "grey75") +
            geom_point(shape = 21, aes(size = AF), color = "grey75") +
            theme_bw() + facet_grid(state~., scale = "free_y", space = "free_y") +
            theme(strip.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),  panel.grid.minor.y = element_blank()) +
            scale_x_date(date_breaks="2 months", date_labels = "%b %Y") +
            scale_y_continuous(n.breaks=3) +
            ylab("") + xlab("Sample Date") +
            ggtitle(paste(moi.nuc)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_line(color = "grey75"))

        #####################################################
        # combine all plot
        plot <- p / plot_spacer() / (r + ggtitle("population weighted average per week")) / plot_spacer() / (q2 + ggtitle("ratio of positve wwtp sample per sliding bi-week")) + plot_layout(guides = "collect", heights = c(10, -.5, .75, -.5, 1.5)) & scale_fill_fermenter(limits = range(c(mutation.dtm$value.freq, dt2$AF), na.rm = TRUE), na.value = "white", palette = "YlOrRd", direction = 1) & scale_size_binned_area(limits = range(c(mutation.dtm$value.freq, dt2$AF), na.rm = TRUE), na.value = .5)

        filename <- paste0(outdir, "/figs/nucleotides",  paste('/Nuc', moi.nuc, sep="_"), ".", opt$graph)
        ggsave(filename = filename, plot = plot, height = 9, width = 9)

  }
}

writeLines(paste("PROGRESS: DONE"))
timestamp()
