
## load libraries
suppressMessages(library("tidyverse"))
suppressMessages(library("data.table"))
suppressMessages(library("patchwork"))
suppressMessages(library("ggpubr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("lubridate"))
suppressMessages(library("optparse"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("dendextend"))
suppressMessages(library("circlize"))
suppressMessages(library("NbClust"))
suppressMessages(library("gslnls"))
suppressMessages(library("cowplot"))
suppressMessages(library("stylo"))
suppressMessages(library("odbc"))
suppressMessages(library("DBI"))

timestamp()

# get location of wwtp

# connect to data base
con <- dbConnect(odbc(),
  Driver = "ODBC Driver 17 for SQL Server",
  Server = "SARSCOV2DB\\SARSCOV2;",
  Database = "SARS_CoV2_SeqDB",
  UID = "CeMM_Reporter",
  PWD = "cemm_reporter"
)
#dbListTables(con)

# get tables
sewagePlants_db = setDT(dbReadTable(con, "sewagePlants"))

# define sites of special intersts
poi = c(22598:22600, 22628:22630, 22892:22894, 22895:22897, 22910:22912, 22898:22900, 22916:22918, 22940:22942, 23018:23020, 23030:23032, 23042:23044, 23039:23041)
#poi = 22553:23155

poi.alt = c(23603:23605)

# get Options
option_list = list(
  make_option(c("--dir"), type="character", default="ExampleOutput",
              help="Directory to write results [default %default]", metavar="character"),
  make_option(c("--metadata"), type="character", default="data/metaDataSub.tsv",
              help="Path to meta data input file [default %default]", metavar="character"),
  make_option(c("--data2"), type="character", default="data/mutationDataSub_sparse.tsv",
              help="Path to data input file in deprecated sparse matrix format [default %default]", metavar="character"),
  make_option(c("--data"), type="character", default="data/mutationDataSub.tsv.gz",
              help="Path to data input file in tidy table format [default %default]", metavar="character"),
  make_option(c("--inputformat"), type="character", default="sparse",
              help="Input data format. Should be 'tidy' or 'sparse' [default %default]", metavar="character"),
  make_option(c("--growthlimit"), type="double", default=0.005,
              help="Mutations with smaller growth rate are ignored [default %default]", metavar="character"),
  make_option(c("--periodend"), type="character", default=as.character(Sys.Date()),
              help="End of analysis period as date in %Y-%M-%D format [default %default]", metavar="character"),
  make_option(c("--graph"), type="character", default="pdf",
              help="Fileformate of produced graphics. Select from pdf, png [default %default]", metavar="character"),
  make_option(c("--periodlength"), type="integer", default=2*30.5,
              help="Duration of analysis perid in days [default %default]", metavar="character"),
  make_option(c("--indels"), type="logical", default=FALSE,
              help="Should be indels be considered [default %default]", metavar="character"),
  make_option(c("--debug"), type="logical", default="FALSE",
              help="Toggle to use run with provided example data [default %default]", metavar="character"),
  make_option(c("--verbose"), type="logical", default="FALSE",
              help="Toggle to report more output [default %default]", metavar="character")

);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (0){

  opt$metadata = "data/metaData_general.csv";
  opt$data = "data/mutationData_DB_NationMonitoringSites.tsv.gz";
  opt$growthlimit = 0.001
  opt$periodend = "2022-12-12"
  opt$dir = paste0("Analysis_", opt$periodend)
  opt$inputformat = "tidy";
  opt$periodlength = 61;
  opt$indels = FALSE
  opt$verbose = TRUE

  ## Rscript --vanilla scripts/vaquera.r
  ## --dir Analysis_${DATE}
  ## --metadata data/metaData_general.csv
  ## --data data/mutationData_DB_NationMonitoringSites.tsv.gz
  ## --inputformat tidy
  ## --growthlimit $GROWTHLIMIT #
  ## --periodlength $PERIODLENGTH
  ## --indels FALSE --verbose FALSE
}


#####################################################
####  parameter setting for interactive debugging
if(opt$debug){
  print("Warning: command line option overwritten")
}
#####################################################

## define functions

options(warn=-1)

`%notin%` <- Negate(`%in%`)

# z-score normalisation function
zscore_norm <- function(x) {
    (x-mean(x))/sd(x)
}

n_topest <- function(x, n) {
  #Sort the wages in descending order
  x1 <- sort(x, decreasing = TRUE)
  return(x1[n])
}

decimalDate <- function(x, d){
  strftime(as.POSIXct(as.Date(x)+(d/96)), format="%Y-%m-%d %H:%M", tz="UCT") -> dx
  daysinyear <- leapYear(year(dx))
  year(dx) + yday(dx)/daysinyear + hour(dx)/(daysinyear*24) + minute(dx)/(daysinyear*24*60) -> ddx
  return(ddx)
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

## print parameter to Log

print("##~LOG~PARAMETERS~####################")
print(opt)
print("##~LOG~PARAMETERS~####################")
writeLines("\n\n\n")

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
if( ! dir.exists(paste0(outdir, "/figs/cluster"))){
  dir.create(paste0(outdir, "/figs/cluster"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/singleModels"))){
  dir.create(paste0(outdir, "/figs/singleModels"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/heatmap"))){
  dir.create(paste0(outdir, "/figs/heatmap"), showWarnings = FALSE)
}
if( ! dir.exists(paste0(outdir, "/figs/mutations"))){
  dir.create(paste0(outdir, "/figs/mutations"), showWarnings = FALSE)
}


## read in meta data
print(paste0("PROGRESS: read and process meta data "))
metaDT       <- fread(file = opt$metadata)
unique(metaDT$BSF_sample_name) -> sampleoi

## read in mutations data
if(opt$inputformat == "sparse"){
  print(paste0("PROGRESS: read AF data (deprecated file format!) "))
  timestamp()
  sewage_samps <- fread(opt$data2 , header=TRUE, sep="\t" ,na.strings = ".", check.names=TRUE)

  ## get sample names
  sample_names = grep("\\.AF$", names(sewage_samps), value = TRUE)
  sample_names = gsub("\\.AF$","",sample_names)

  ## remove all positions which are not mutations of interest
  unite(sewage_samps, NUC, c(3,2,4), ALT, sep = "", remove = FALSE) -> sewage_samps

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
  sewage_samps.dt$value.depth[is.na(sewage_samps.dt$value.depth)] <- 0
  sewage_samps.dt$value.qual[is.na(sewage_samps.dt$value.qual)] <- 0

  ## remove all samples which are not specified in metadata
  sewage_samps.dt %>% filter(ID %in% sampleoi) -> sewage_samps.dt

  rm(sewage_samps)
}

## generate mapping between nucc mutation and mutation label
sewage_samps.dt %>% mutate(label = paste(NUC, paste0("[", paste(ANN.GENE,ANN.AA, sep =":"), "]"))) %>% select(NUC, label) %>% distinct() -> nuc2label

## remove indes if specified by opt$indels
if (!opt$indels){
  print(paste0("LOG: indels removed due to default --indels FALSE option"))
  sewage_samps.dt %>% filter(nchar(REF) == nchar(ALT)) -> sewage_samps.dt
}

## add shorten meta data to mutation data
left_join(x = sewage_samps.dt, y = metaDT %>% filter(grepl("pass", status)) %>% select(BSF_sample_name, LocationID, sample_date), by = c("ID" = "BSF_sample_name")) -> sewage_samps.dt

## check if all samples in mutation data is also included in meta data
sewage_samps.dt %>% filter(is.na(sample_date)) %>% group_by(ID) %>% summarize(n = n()) -> metadatalessSamples
if (dim(metadatalessSamples)[1] > 0){
  print(paste0("WARNING: ", dim(metadatalessSamples)[1], " samples in data file not specified in metadata file; will be ignored."))
  sewage_samps.dt %>% filter(ID %notin% metadatalessSamples$ID) -> sewage_samps.dt
}

## check if periodend is <= latest date in data; if not reset
if (max(as.Date(sewage_samps.dt$sample_date)) < as.Date(opt$periodend)){
  print(paste0("WARNING: ", opt$periodend, " later than latest sample_time in data; reset to ", max(as.Date(sewage_samps.dt$sample_date))))
  opt$periodend <- max(as.Date(sewage_samps.dt$sample_date))
}

## remove all samples which are not in the (periodend-periodlength):periodend interval
metaDT %>% filter(sample_date <= as.Date(opt$periodend) & sample_date >= as.Date(opt$periodend) - ceiling(opt$periodlength)) -> metaDT

## remove all mutation data which are not in the (periodend-periodlength):periodend interval
sewage_samps.dt %>% filter(sample_date <= as.Date(opt$periodend) & sample_date >= as.Date(opt$periodend) - ceiling(opt$periodlength)) -> sewage_samps.dt

## plot kinetics of all mutations in regions of interst
## filter for mutations which are
##    seen in at >2 wwtp at >2 timepoints
##    AF value between 5% and 95% in at least 1 wwtp
##
## use https://codon2nucleotide.theo.io/ to map between nuc and codons
sewage_samps.dt %>% filter(POS %in% poi) %>% group_by(NUC, LocationID) %>% mutate(n = n(), ma = max(value.freq), mi = min(value.freq)) %>% filter(n > 2 & ma > 0.05 & mi < 0.95) %>% ungroup() %>% group_by(NUC) %>% mutate(n = length(unique(LocationID))) %>% filter(n > 2) %>% ggplot(aes(x = as.Date(sample_date), y = value.freq, color = paste(ANN.GENE,ANN.AA,sep=":"))) + geom_point(size = 2, alpha = 0.66) + geom_smooth(formula = y ~ x, method = "glm", method.args = list(family = binomial(link='logit')), se = F, color = "blue") + facet_grid(paste(paste(ANN.GENE, ANN.AA, sep=":"),NUC, sep="\n")~LocationID, scale = "free_y")  + theme_bw() + theme(legend.position="none", strip.text.x = element_text(size = 5), axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) + scale_x_date(date_breaks="2 weeks", date_labels="%d %b") + ylab("Mutation frequency") + xlab("Sampling Date") + ggtitle("Mutations in Codons of Interest per WWTP") + scale_color_viridis_d(name = "AA Mutation", option = "H") -> p
f<-1.8
ggsave(filename = paste0(opt$dir, "/figs/Mutations_of_Interest.", opt$graph), plot = p, width = f*16, height = f*9)


## plot kinetics of all alternative mutatations in regions of interst
## use https://codon2nucleotide.theo.io/ to map between nuc and codons
sewage_samps.dt %>% filter(POS %in% poi.alt) %>% group_by(NUC, LocationID) %>% mutate(n = n()) %>% filter(n >= 1) %>% ungroup() %>% group_by(NUC) %>% mutate(n = length(unique(LocationID))) %>% filter(n >= 1) %>% ggplot(aes(x = as.Date(sample_date), y = value.freq, color = paste(ANN.GENE,ANN.AA,sep=":"))) + geom_point(size = 2, alpha = 0.66) + geom_smooth(formula = y ~ x, method = "glm", method.args = list(family = binomial(link='logit')), se = F, color = "blue") + facet_grid(paste(paste(ANN.GENE, ANN.AA, sep=":"),NUC, sep="\n")~LocationID, scale = "free_y")  + theme_bw() + theme(legend.position="none", strip.text.x = element_text(size = 5), axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) + scale_x_date(date_breaks="2 weeks", date_labels="%d %b") + ylab("Mutation frequency") + xlab("Sampling Date") + ggtitle("Mutations in Codons of Interest per WWTP") + scale_color_viridis_d(name = "AA Mutation", option = "H") -> p
f<-0.75
ggsave(filename = paste0(opt$dir, "/figs/Mutations_of_Alternative_Interest.", opt$graph), plot = p, width = f*16*2.3, height = f*9*0.75)

## data collector global
data.table(mutation = character(), growth_pred = numeric(), inflection_pred = numeric(), label = character(), cluster = numeric(), LocationID = character()) -> cluster_collector
data.table(mutation = character(), growth_pred = numeric(), inflection_pred = numeric(), LocationID = character()) -> prediction_collector_global

## loop over each WWTP
for (llocationID in sort(unique(sewage_samps.dt$LocationID))){
  print(paste0("PROGRESS: analyze ", llocationID))
  timestamp()

  ## data collector per WWTP
  data.table(mutation = character(), growth_pred = numeric(), inflection_pred = numeric()) -> prediction_collector


  sewage_samps.dt %>% filter(LocationID == llocationID) -> dt
  ## take care of mean function here (same data same location issue!!!!)
  data.table::dcast(data = dt, formula = sample_date~NUC, value.var = "value.freq", fill = 0, fun.aggregate = mean)  -> dt

  ## keep only mutations with AF = 0 in the first OR last l=3 timepoints
  dim(dt)
  l = 1
  as.vector(c(TRUE, colSums(dt[1:l,2:dim(dt)[2]]) != 0)) -> selected_columns_first
  as.vector(c(TRUE, colSums(dt[(dim(dt)[1]-l+1):(dim(dt)[1]),2:dim(dt)[2]]) != 0)) -> selected_columns_last
  selected_columns_first | selected_columns_last -> selected_columns
  as.data.table(dt) -> dt
  dt[,..selected_columns] -> dt
  dim(dt)

  print("LOG: ## test each mutation individual for growth")
  ## to play with S:G346T: which(colnames(dt) == "G22599C") -> moiidx
  for (moiidx in seq(from = 2, to = length(colnames(dt)))){
    dtoi <- dt
    mname <- colnames(dtoi)[moiidx]
    colnames(dtoi)[moiidx] <- "y"
    colnames(dtoi)[1] <- "t"
    subset(dtoi, select = c("t", "y")) -> dtoi

    #### remove zeros assuming dropouts
    if(FALSE){
      print("WARNING: zeros are removed to anticipate drop-outs!")
      dtoi[dtoi$y > 0] -> dtoi
    }

    #### ~nls: fixed asymptote to 1
    dtoi$d <- as.numeric(as.Date(dtoi$t)-as.Date(dtoi$t[1]))
    a <- 1
    fo <- y ~ a / (1 + exp(-b * (d-c)))
    #model <- tryCatch(nls(fo, start = list(b = 0.05, c = 2*max(dtoi$d)/3), data = dtoi),error=function(e) NA, warning=function(w) NA)
    model <- tryCatch(gsl_nls(fo, start = list(b = 0.05, c = 2*max(dtoi$d)/3), data = dtoi),error=function(e) NA, warning=function(w) NA)
    if(identical(NA, model)){
      model <- tryCatch(nls(fo, start = list(b = 0.2, c = max(dtoi$d)/2), data = dtoi),error=function(e) NA, warning=function(w) NA)
    }
    if(identical(NA, model)){
      model <- tryCatch(gsl_nls(fo, start = list(b = -0.05, c = 2*max(dtoi$d)/3), data = dtoi),error=function(e) NA, warning=function(w) NA)
    }
    if(identical(NA, model)){
      model <- tryCatch(gsl_nls(fo, start = list(b = -0.2, c = max(dtoi$d)/2), data = dtoi),error=function(e) NA, warning=function(w) NA)
    }

    if(identical(NA, model)){
       #print(paste("Log: error in", mname, "@", llocationID))
       rbind(prediction_collector, data.table(mutation = mname, growth_pred = NA, inflection_pred = NA)) -> prediction_collector
    } else{
       #print(paste("Log: convergence in", mname, "@", llocationID))

       # predict the relative growth in the last 14 days
       dtoi$predict <- predict(model)
       lastweek_prediction <- predict(model, newdata = data.frame(t = c(as.Date(opt$periodend)-14, opt$periodend), d = as.numeric(c(as.Date(opt$periodend)-14, opt$periodend)-as.Date(dtoi$t[1]))))
       lastweek_growth <- (lastweek_prediction[2] - lastweek_prediction[1])/2

       #coefficients(model)[1] -> growth_pred
       lastweek_growth -> growth_pred
       coefficients(model)[2] -> inflection_pred

       ## set inflection_pred and growth_pred to NA growth_pred confidence interval overlaps 0
       summary(model) -> model_summary
       model_summary$coefficients[ , 4][1] -> pvalue_growth_pred
       model_summary$coefficients[ , 4][2] -> pvalue_inflection_pred
       #CI_growth_pred <- tryCatch(suppressMessages(confint(model, "b", level = 0.90)),error=function(e) NA, warning=function(w) NA)
       CI_growth_pred <- tryCatch(suppressMessages(confintd(model, "b", level = 0.90)),error=function(e) NA, warning=function(w) NA)
       CI_inflection_pred <- tryCatch(suppressMessages(confintd(model, "c", level = 0.90)),error=function(e) NA, warning=function(w) NA)

       #if (any(is.na(CI_growth_pred)) | (min(CI_growth_pred) < 0 & max(CI_growth_pred) > 0)){
       if (any(is.na(CI_growth_pred)) | (min(CI_growth_pred[2:3]) < 0 & max(CI_growth_pred[2:3]) > 0)){
         if (opt$verbose){
            print(paste("LOG: growth_pred for", mname, "in", llocationID, "set to NA since CI spanning 0 [", paste(paste(colnames(CI_growth_pred)[2], signif(CI_growth_pred[2], digits = 2), sep = ": "), paste(colnames(CI_growth_pred)[1], signif(CI_growth_pred[1], digits = 2), sep = ": "), paste(colnames(CI_growth_pred)[3], signif(CI_growth_pred[3], digits = 2), sep = ": "), sep = "; "), "]"))
         }
         growth_pred = NA
         inflection_pred = NA
       }


       rbind(prediction_collector, data.table(mutation = mname, growth_pred = growth_pred, inflection_pred = inflection_pred)) -> prediction_collector

       if (opt$verbose & !is.na(growth_pred)){
          melt(dtoi[,-3], id.vars = "t") %>% ggplot(aes(x = t, y = value, color = variable)) + geom_point(size = 3, alpha = 0.66, aes(shape = variable, size = variable)) + geom_smooth(formula = y ~ x, method = "glm", method.args = list(family = "binomial"), se = FALSE) + theme_bw() + ggtitle(paste(llocationID, nuc2label$label[nuc2label$NUC == mname], sep = ": ")) + theme(legend.position="bottom", legend.title = element_blank()) + ylab("mutation frequency [1/1]") + xlab("Sample Date") + geom_vline(xintercept=as.Date(dtoi$t[1])-dtoi$d[1]+as.numeric(inflection_pred), linetype = "dotted")  + scale_color_brewer(palette = "Set1", breaks = c("y", "predict"), labels = c("Observations", "Model Prediction")) + scale_shape(breaks = c("y", "predict"), labels = c("Observations", "Model Prediction")) + annotate(x = min(dtoi$t)+(opt$periodlength/10), y = max(dtoi$predict), geom = "text", label = paste("Growth per week:\n", signif(growth_pred, digits = 2)))-> p
          ggsave(filename = paste0(opt$dir, "/figs/singleModels/", "plot_", llocationID, "_", mname, ".", opt$graph), p, width = 5, height = 5)
       }
    }
    rm(dtoi, model)
  }

  print("LOG: ## cluster according 2d scatter plot")
  prediction_collector %>% select(mutation, growth_pred, inflection_pred) %>% filter(!is.na(growth_pred)) %>% mutate(growth_pred = growth_pred, z_growth_pred = zscore_norm(growth_pred), z_inflection_pred = zscore_norm(inflection_pred) ) -> cDat

  #   cluster if >= 2 mutations and wwtp
  cDat %>% filter(abs(growth_pred) >= opt$growthlimit) -> cDat

  # add mutation label to data
  left_join(x = cDat, y = nuc2label, by = c("mutation" = "NUC")) -> cDat

  if ( all(dim(cDat) >= 2)) {
    d <- dist(cDat[,4:5], method = "binary")
    hls <- tryCatch(hclust(d, method = "ward.D2"),error=function(e) NA, warning=function(w) NA)
    if(identical(NA, hls)){
      cDat$cluster <- 1:dim(cDat)[1]
      cDat$cluster <- as.factor(cDat$cluster)
    } else{
      cluster <- cutree(hls, h = 1)
      cDat <- cbind(cDat, cluster)
      cDat$cluster <- as.factor(cDat$cluster)
    }
  } else{
    cDat$cluster <- 1:dim(cDat)[1]
    cDat$cluster <- as.factor(cDat$cluster)
  }
  if(all(dim(cDat) > 0)){

    cDat$LocationID <- llocationID


    cDat %>% group_by(cluster) %>% mutate(cluster_title = paste(paste("cid", cur_group_id(), sep = ":"), paste("size", length(unique(mutation)), sep = ":"), paste("r", signif(mean(growth_pred), digits = 3), sep=":"))) %>% ungroup() %>% ggplot(aes(x = cluster_title, y=growth_pred, fill=cluster_title)) + geom_boxplot(aes(color = cluster_title), alpha = 0.5) + geom_point(color = "black") + geom_text_repel(aes(label = label, color = cluster_title), size = 1.5, max.overlaps = 40) + theme_bw() + ylab("Growth Rate [1/w]") + xlab("") + theme(legend.position="none") + coord_flip() + scale_fill_viridis_d(option = "H") + scale_color_viridis_d(option = "H") -> p

    left_join(x = melt(dt, id.vars = "sample_date"), y = cDat, by = c("variable" = "mutation")) %>% filter(!is.na(cluster)) %>% group_by(cluster) %>% mutate(cluster_title = paste(paste("cid", cur_group_id(), sep = ":"), paste("size", length(unique(variable)), sep = ":"), paste("r", signif(mean(growth_pred), digits = 3), sep=":"))) %>% ungroup() %>% ggplot(aes(x = sample_date, y = value, group = variable)) + geom_smooth(method = "glm", method.args = list(family=binomial(link='logit')), se = FALSE) + geom_point(aes(color = cluster_title)) + theme_bw() + facet_wrap(.~cluster_title) + ggtitle(llocationID) + theme(legend.position="none") + scale_x_date(date_breaks = paste(ceiling(opt$periodlength/2), "days"), date_labels = "%d-%b")  + ylab("Mutation Frequency [1/1]") + xlab("Sample Date") + scale_fill_viridis_d(option = "H") + scale_color_viridis_d(option = "H") -> q

    plot_grid(p, q) -> pq

    left_join(x = melt(dt, id.vars = "sample_date"), y = cDat, by = c("variable" = "mutation")) %>% filter(!is.na(cluster)) %>% group_by(cluster) %>% mutate(cluster_title = paste(paste("cid", cur_group_id(), sep = ":"), paste("size", length(unique(variable)), sep = ":"), paste("r", signif(mean(growth_pred), digits = 3), sep=":"))) %>% group_by(cluster_title) %>% summarize(n = n()) -> N
    f = 1.5

    #save_plot(filename=paste0(opt$dir, "/figs/cluster/", "plot_", llocationID, "_clustering.", opt$graph), plot = pq)
    ggsave(filename=paste0(opt$dir, "/figs/cluster/", "plot_", llocationID, "_clustering.", opt$graph), p + q, width = 3+(length(N$n)*f), height = 3+(length(N$n)*f))

    cDat$LocationID <- llocationID
    rbind(cDat[,c(1,2,3,6,7,8)], cluster_collector) -> cluster_collector

    prediction_collector$LocationID <- llocationID
    rbind(prediction_collector_global, prediction_collector) -> prediction_collector_global

    cDat %>% group_by(cluster) %>% mutate(cluster_title = paste(paste("cid", cur_group_id(), sep = ":"), paste("size", length(unique(mutation)), sep = ":"), paste("r", signif(mean(growth_pred), digits = 3), sep=":"))) %>% ungroup() %>% group_by(cluster_title) %>% summarize(covstring = paste0("[", n(), "-of:", paste(unique(mutation), collapse = ", "), "]")) -> covspectrum_searchString
    fwrite(file=paste0(opt$dir, "/figs/cluster/", "plot_", llocationID, "_clustering_covspecString.csv"), x = covspectrum_searchString, sep = "\t")
  }
}

## dump data into file
left_join(x = prediction_collector_global, y = nuc2label, by = c("mutation" = "NUC")) -> prediction_collector_global_labeled
fwrite(file = paste0(opt$dir, "/growthrate_prediction.csv"), sep = "\t", prediction_collector_global_labeled)
fwrite(file = paste0(opt$dir, "/clusterAssignment.csv"), sep = "\t", cluster_collector)

## dump data for poi into file
prediction_collector_global_labeled %>% mutate(pos = gsub("\\D", "", mutation, perl = TRUE)) %>% filter(pos %in% poi) %>% filter(!is.na(growth_pred)) %>% group_by(mutation, label) %>% summarize(mean = mean(growth_pred), median = median(growth_pred), max = max(growth_pred), n2_val = n_topest(growth_pred, 2), n3_val = n_topest(growth_pred, 3), nr_growing = sum(growth_pred > 0)) %>% arrange(desc(nr_growing)) -> summarised_mutation_growth_poi
fwrite(file = paste0(opt$dir, "/growthrate_per_mutation_pois.csv"), sep = "\t", summarised_mutation_growth_poi)


prediction_collector_global_labeled %>% filter(!is.na(growth_pred)) %>% group_by(mutation, label) %>% summarize(mean = mean(growth_pred), median = median(growth_pred), max = max(growth_pred), n2_val = n_topest(growth_pred, 2), n3_val = n_topest(growth_pred, 3), nr_growing = sum(growth_pred > 0)) %>% arrange(desc(nr_growing)) -> summarised_mutation_growth
fwrite(file = paste0(opt$dir, "/growthrate_per_mutation.csv"), sep = "\t", summarised_mutation_growth)

# select mutations if:
##    median growth rate > 0.05
##    top3 growthrate > 0.05
##    in more than 2 wwtp positive growth
summarised_mutation_growth %>% filter(nr_growing > 2) %>% filter(median > 0.05 & n3_val > 0.05) -> selected_mutations

# select poi mutations if:
##    median growth rate > 0.01
##    top3 growthrate > 0.01
##    in more than 0 wwtp positive growth
summarised_mutation_growth_poi %>% filter(nr_growing > 0) %>% filter(median > 0.01 & n3_val > 0.01) -> selected_poi_mutations

## add mutation score to selected poi mutations
selected_poi_mutations %>% rowwise() %>% mutate(aa = substr(label, regexpr("S:\\w+", label)[1]+3, regexpr("S:\\w+", label)[1]+attr(regexpr("S:\\w+", label), "match.length")-1)) -> selected_poi_mutations
fread(file="auxData/mutations_scores.csv") -> mutations_Scores
left_join(x = selected_poi_mutations, y = mutations_Scores, by = "aa") -> selected_poi_mutations
selected_poi_mutations %>% ggplot(aes(x = mean_mut_escape, y = delta_bind)) + geom_point(aes(size = nr_growing, fill = median), alpha = 0.5, shape = 21) + geom_text_repel(aes(label = label), max.overlaps = 40, min.segment.length = 0) + theme_bw() + xlab("Mean Mutation Escape to BA.5 convalescents sera (Cao et al.)") + ylab("Delta ACE2 binding in BA.2 background (Bloom et al.)") + scale_size(range = c(4,8), name = "Number of \nWWTP w/\nsig. growth") + scale_fill_viridis_b(option = "E", name = "Median\ngrowth\nper week") + ggtitle("Characterization of growing RBD Mutations") -> p
ggsave(filename = paste0(opt$dir, "/figs/RBD_Mutation_Characterisation.", opt$graph), plot = p, width = 9, height = 9)


# plot AF for each selected mutation across all WWTP
if(opt$verbose == TRUE){
  for (mmutation in unique(selected_mutations$mutation)) {
    sewage_samps.dt %>% filter(NUC == mmutation) -> pDat
    pDat %>% group_by(LocationID) %>% mutate(N = length(unique(sample_date))) -> pDat
    if(max(pDat$N) > 2 ){
      left_join(x = pDat, y = prediction_collector_global, by = c("NUC" = "mutation", "LocationID" = "LocationID")) -> pDat
      ggplot(data = pDat, aes(x = sample_date, y = value.freq, color = LocationID, group = LocationID,)) + geom_point(size = 3, alpha = 0.75) + geom_smooth(formula = y ~ x, method = "glm", method.args = list(family = binomial(link='logit')), se = F, color = "blue") + theme_bw() + theme(legend.position="none") + scale_color_viridis_d(option = "H") + scale_fill_viridis_d(option = "H") + ylab("Allel Frequency [1/1]") + xlab("Sample Date") + facet_wrap(~paste(LocationID, paste0(" | R: ", signif(pDat$growth_pred, digits = 2)))) + ggtitle(paste(paste(unique(paste(pDat$ANN.GENE, pDat$ANN.AA, sep = ":")), collapse = ", "), "|", mmutation)) -> p
      f = 1.2
      ggsave(filename = paste0(opt$dir, "/figs/mutations/plot_", mmutation, "_", paste(unique(paste(pDat$ANN.GENE, pDat$ANN.AA, sep = "_")), collapse = "__"), ".", opt$graph), plot = p, width = f*16, height = f*9)
    }
  }
}



# filter mutation which are seen only once sig. modulated in only one location
prediction_collector_global_labeled %>% filter(!is.na(growth_pred)) %>% group_by(mutation) %>% mutate(n=n()) %>% filter(n > 1) -> prediction_collector_global_labeled

#### produce heatmap of all mutations
if(!any(dim(prediction_collector_global) == 0)){
  data.table::dcast(data = prediction_collector_global_labeled, formula = label~LocationID, value.var = "growth_pred", fun.aggregate = mean)  -> growth_prediction_collector_global.dt

  # generate cluster matrix from data
  growth_prediction_collector_global.dt$label -> mat_rownames

  grep("label", colnames(growth_prediction_collector_global.dt), invert = T) -> selected_columns
  growth_prediction_collector_global.dt[, selected_columns] -> growth_prediction_collector_global.dt
  colnames(growth_prediction_collector_global.dt) -> mat_colnames

  as.matrix(growth_prediction_collector_global.dt) -> mat
  rownames(mat) <- mat_rownames
  colnames(mat) <- mat_colnames

  # cluster all mutations
  if (all(dim(mat) >= 2)){

    # set NA to zero
    mat2cluster <- mat
    mat2cluster[is.na(mat2cluster)] <- 0
    mat2cluster4all <- mat2cluster

    clust_mut <- tryCatch(hclust(dist(mat2cluster, method = "binary"), method="ward.D2"),error=function(e) e, warning=function(w) w)
    if(any(grepl("warning|error", class(clust_mut)))){
      print("Warning: hclust mutation failed with method binary")
      clust_mut <- tryCatch(hclust(dist(mat2cluster), method="ward.D2"),error=function(e) e, warning=function(w) w)
      if(any(grepl("warning|error", class(clust_mut)))){
        print("Warning: hclust mutation failed with dist")
      }
    }
    data.frame(cluster = cutree(clust_mut, h = 2)) ->  clust_mut_split


    # cluster locations based on distance
    sewagePlants_db %>% select(LocationID_coronA, coordinates_east, coordinates_north) %>% filter(!is.na(coordinates_east) & !is.na(coordinates_north))-> wwtp_coordinates
    wwtp_coordinates %>% filter(grepl("Meining", LocationID_coronA)) %>% mutate(LocationID_coronA = "Liechtenstein", coordinates_north = 385597) -> liefoo
    rbind(wwtp_coordinates, liefoo) -> wwtp_coordinates

    wwtp_coordinates %>% filter(LocationID_coronA %in% colnames(mat)) -> wwtp_coordinates
    wwtp_coordinates$LocationID_coronA -> wwtp_rownames
    as.matrix(wwtp_coordinates[,2:3]) -> wwtp_mat
    rownames(wwtp_mat) <- wwtp_rownames

    clust_loc <- tryCatch(hclust(dist(wwtp_mat, method = "euclidean")),error=function(e) e, warning=function(w) w)

    data.frame(cluster = cutree(clust_loc, k = 9)) ->  clust_loc_split

    # make heatmap
    col_fun = colorRamp2(c(min(mat, na.rm = T), 0, max(mat, na.rm = T)), c("blue", "white", "red"))
    if(opt$graph == "pdf"){
      pdf(file = paste0(opt$dir, "/figs/heatmap/", "growthrates.pdf"))
    }
    if(opt$graph == "png"){
      png(file = paste0(opt$dir, "/figs/heatmap/", "growthrates.png"))
    }



    Heatmap(mat2cluster, name = "Growth Rate", col = col_fun, use_raster = TRUE, cluster_rows = clust_mut, cluster_columns = clust_loc, show_row_names=T, show_column_names=T,row_title = "Mutations", column_title = "WWTP", row_split = max(length(unique(clust_mut_split$cluster)), 2), column_split = max(length(unique(clust_loc_split$cluster)), 2), row_names_gp = gpar(fontsize = 2), column_names_gp = gpar(fontsize = 5), border = TRUE)-> hm
    hm <- draw(hm)
    row_order(hm) -> mutation_sets
    lapply(mutation_sets, function(x){ as.data.table(rownames(mat)[x]) %>% separate(V1, c("nuc", "aa")) %>% select(nuc) -> y; paste0("[", ceiling(length(x)/2), "-of:", paste(unique(y$nuc), collapse = ", "), "]")   }) -> mutation_sets_for_covspectrum
    paste0("Set", seq_along(mutation_sets_for_covspectrum)) -> names(mutation_sets_for_covspectrum)
    fwrite(file = paste0(opt$dir, "/figs/heatmap/cluster2covspectrum.csv"), sep = "\t", mutation_sets_for_covspectrum)
    dev.off()
  }
}


#### produce heatmap of all selected mutations
if(!any(dim(prediction_collector_global) == 0)){
  prediction_collector_global_labeled %>% filter(mutation %in% unique(selected_mutations$mutation)) -> prediction_collector_global_labeled_selected
  data.table::dcast(data = prediction_collector_global_labeled_selected, formula = label~LocationID, value.var = "growth_pred", fun.aggregate = mean)  -> growth_prediction_collector_global_selected.dt

  # generate cluster matrix from data
  growth_prediction_collector_global_selected.dt$label -> mat_rownames

  grep("label", colnames(growth_prediction_collector_global_selected.dt), invert = T) -> selected_columns
  growth_prediction_collector_global_selected.dt[, selected_columns] -> growth_prediction_collector_global_selected.dt
  colnames(growth_prediction_collector_global_selected.dt) -> mat_colnames

  as.matrix(growth_prediction_collector_global_selected.dt) -> mat
  rownames(mat) <- mat_rownames
  colnames(mat) <- mat_colnames

  # cluster all mutations
  if (all(dim(mat) >= 2)){

    # set NA to zero
    mat2cluster <- mat
    mat2cluster[is.na(mat2cluster)] <- 0

    clust_mut <- tryCatch(hclust(dist(mat2cluster, method = "binary"), method="ward.D2"),error=function(e) e, warning=function(w) w)
    if(any(grepl("warning|error", class(clust_mut)))){
      print("Warning: hclust mutation failed with method binary")
      clust_mut <- tryCatch(hclust(dist(mat2cluster), method="ward.D2"),error=function(e) e, warning=function(w) w)
      if(any(grepl("warning|error", class(clust_mut)))){
        print("Warning: hclust mutation failed with dist")
      }
    }
    data.frame(cluster = cutree(clust_mut, h = 1.5)) ->  clust_mut_split


    # cluster locations based on distance
    sewagePlants_db %>% select(LocationID_coronA, coordinates_east, coordinates_north) %>% filter(!is.na(coordinates_east) & !is.na(coordinates_north))-> wwtp_coordinates
    wwtp_coordinates %>% filter(grepl("Meining", LocationID_coronA)) %>% mutate(LocationID_coronA = "Liechtenstein", coordinates_north = 385597) -> liefoo
    rbind(wwtp_coordinates, liefoo) -> wwtp_coordinates


    wwtp_coordinates %>% filter(LocationID_coronA %in% colnames(mat)) -> wwtp_coordinates
    wwtp_coordinates$LocationID_coronA -> wwtp_rownames
    as.matrix(wwtp_coordinates[,2:3]) -> wwtp_mat
    rownames(wwtp_mat) <- wwtp_rownames

    clust_loc <- tryCatch(hclust(dist(wwtp_mat, method = "euclidean")),error=function(e) e, warning=function(w) w)

    data.frame(cluster = cutree(clust_loc, k = 9)) ->  clust_loc_split

    # make heatmap
    col_fun = colorRamp2(c(min(mat, na.rm = T), 0, max(mat, na.rm = T)), c("blue", "white", "red"))
    if(opt$graph == "pdf"){
      pdf(file = paste0(opt$dir, "/figs/heatmap/", "growthrates_selected.pdf"))
    }
    if(opt$graph == "png"){
      png(file = paste0(opt$dir, "/figs/heatmap/", "growthrates_selected.png"))
    }

    Heatmap(mat2cluster, name = "Growth Rate", col = col_fun, use_raster = TRUE, cluster_rows = clust_mut, cluster_columns = clust_loc, show_row_names=T, show_column_names=T,row_title = "Mutations", column_title = "WWTP", row_split = max(length(unique(clust_mut_split$cluster)), 2), column_split = max(length(unique(clust_loc_split$cluster)), 2), row_names_gp = gpar(fontsize = 2), column_names_gp = gpar(fontsize = 5), border = TRUE)-> hm
    hm <- draw(hm)
    row_order(hm) -> mutation_sets
    lapply(mutation_sets, function(x){ as.data.table(rownames(mat)[x]) %>% separate(V1, c("nuc", "aa")) %>% select(nuc) -> y; paste0("[", ceiling(length(x)/2), "-of:", paste(unique(y$nuc), collapse = ", "), "]")   }) -> mutation_sets_for_covspectrum
    paste0("Set", seq_along(mutation_sets_for_covspectrum)) -> names(mutation_sets_for_covspectrum)
    fwrite(file = paste0(opt$dir, "/figs/heatmap/clusterselected2covspectrum.csv"), sep = "\t", mutation_sets_for_covspectrum)
    dev.off()
  }
}


#### produce heatmap of all POI mutations
if(!any(dim(prediction_collector_global) == 0)){
  prediction_collector_global_labeled %>% filter(mutation %in% unique(selected_poi_mutations$mutation)) -> prediction_collector_global_labeled_selected
  data.table::dcast(data = prediction_collector_global_labeled_selected, formula = label~LocationID, value.var = "growth_pred", fun.aggregate = mean)  -> growth_prediction_collector_global_selected.dt

  # generate cluster matrix from data
  growth_prediction_collector_global_selected.dt$label -> mat_rownames

  grep("label", colnames(growth_prediction_collector_global_selected.dt), invert = T) -> selected_columns
  growth_prediction_collector_global_selected.dt[, selected_columns] -> growth_prediction_collector_global_selected.dt
  colnames(growth_prediction_collector_global_selected.dt) -> mat_colnames

  as.matrix(growth_prediction_collector_global_selected.dt) -> mat
  rownames(mat) <- mat_rownames
  colnames(mat) <- mat_colnames

  # cluster all mutations
  if (all(dim(mat) >= 2)){

    # set NA to zero
    mat2cluster <- mat
    mat2cluster[is.na(mat2cluster)] <- 0

    clust_mut <- tryCatch(hclust(dist(mat2cluster, method = "binary"), method="ward.D2"),error=function(e) e, warning=function(w) w)
    if(any(grepl("warning|error", class(clust_mut)))){
      print("Warning: hclust mutation failed with method binary")
      clust_mut <- tryCatch(hclust(dist(mat2cluster), method="ward.D2"),error=function(e) e, warning=function(w) w)
      if(any(grepl("warning|error", class(clust_mut)))){
        print("Warning: hclust mutation failed with dist")
      }
    }
    data.frame(cluster = cutree(clust_mut, h = 1)) ->  clust_mut_split


    # cluster locations based on distance
    sewagePlants_db %>% select(LocationID_coronA, coordinates_east, coordinates_north) %>% filter(!is.na(coordinates_east) & !is.na(coordinates_north))-> wwtp_coordinates
    wwtp_coordinates %>% filter(grepl("Meining", LocationID_coronA)) %>% mutate(LocationID_coronA = "Liechtenstein", coordinates_north = 385597) -> liefoo
    rbind(wwtp_coordinates, liefoo) -> wwtp_coordinates

    wwtp_coordinates %>% filter(LocationID_coronA %in% colnames(mat)) -> wwtp_coordinates
    wwtp_coordinates$LocationID_coronA -> wwtp_rownames
    as.matrix(wwtp_coordinates[,2:3]) -> wwtp_mat
    rownames(wwtp_mat) <- wwtp_rownames

    clust_loc <- tryCatch(hclust(dist(wwtp_mat, method = "euclidean")),error=function(e) e, warning=function(w) w)

    data.frame(cluster = cutree(clust_loc, k = 9)) ->  clust_loc_split

    # make heatmap
    col_fun = colorRamp2(c(min(mat, na.rm = T), 0, max(mat, na.rm = T)), c("blue", "white", "red"))
    if(opt$graph == "pdf"){
      pdf(file = paste0(opt$dir, "/figs/heatmap/", "growthrates_selected_poi.pdf"))
    }
    if(opt$graph == "png"){
      png(file = paste0(opt$dir, "/figs/heatmap/", "growthrates_selected_poi.png"))
    }


    Heatmap(mat2cluster, name = "Growth Rate", col = col_fun, use_raster = TRUE, cluster_rows = clust_mut, cluster_columns = clust_loc, show_row_names=T, show_column_names=T,row_title = "Mutations", column_title = "WWTP", row_split = max(length(unique(clust_mut_split$cluster)), 2), column_split = max(length(unique(clust_loc_split$cluster)), 2), row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6), border = TRUE)-> hm
    hm <- draw(hm)
    row_order(hm) -> mutation_sets
    lapply(mutation_sets, function(x){ as.data.table(rownames(mat)[x]) %>% separate(V1, c("nuc", "aa")) %>% select(nuc) -> y; paste0("[", ceiling(length(x)/2), "-of:", paste(unique(y$nuc), collapse = ", "), "]")   }) -> mutation_sets_for_covspectrum
    paste0("Set", seq_along(mutation_sets_for_covspectrum)) -> names(mutation_sets_for_covspectrum)
    fwrite(file = paste0(opt$dir, "/figs/heatmap/clusterselectedpoi2covspectrum.csv"), sep = "\t", mutation_sets_for_covspectrum)
    dev.off()
  }
}

##### EXXPLORE DATA


## manually chosen positions of interest
ppoi <- poi

# examine difference between all mutations and manually selected POI mutations
as.data.table(mat2cluster4all) %>% mutate(label = rownames(mat2cluster4all)) %>% data.table::melt() %>% mutate(cat = ifelse(grepl(paste(ppoi, collapse="|"), label), "ppoi", "bg")) -> edt


if( ! dir.exists(paste0(outdir, "/EXPLORE"))){
  dir.create(paste0(outdir, "/EXPLORE"), showWarnings = FALSE)
}

## make a pca plot
as.data.table(mat2cluster4all) %>% mutate(label = rownames(mat2cluster4all)) %>% mutate(cat = ifelse(grepl(paste(poi, collapse="|"), label), "poi", "bg")) -> pdt

library(ggfortify)
pdt[,1:(length(pdt)-2)] -> ppdt
pdt[,(length(pdt)-1):length(pdt)] -> dpdt

prcomp(ppdt) -> pca
cbind(pca$x[,which(colnames(pca$x) == "PC1" | colnames(pca$x) == "PC2")], dpdt) %>% ggplot(aes(x = PC1, y = PC2, color = cat)) + geom_point(aes(alpha = cat)) + theme_bw() + scale_alpha_discrete(range = c(0.4, .8)) -> p
pfilenamebase = "pca"
ggsave(filename = paste0(outdir, "/EXPLORE/POIplot_", pfilenamebase, ".", opt$graph), plot = p)



####

## plot kinetics of all mutatations above TH filter on 3_quantile
TH <- 0.01
edt %>% group_by(label) %>% summarize(value2 = n_topest(value, 3), cat = cat) %>% filter(value2 >= TH) %>% select(label) %>% distinct() %>% separate(col = "label", into = c("nnuc", "ggene", "aaa")) -> loi

sewage_samps.dt %>% filter(NUC %in% loi$nnuc) %>% group_by(NUC, LocationID) %>% mutate(n = n()) %>% filter(n > 1) %>% ungroup() %>% group_by(NUC) %>% mutate(n = length(unique(LocationID))) %>% filter(n > 1) %>% ungroup %>% group_by(ANN.GENE)  %>% summarize(n = n()) %>% mutate(ANN.GENE = ifelse(is.na(ANN.GENE), "intergenic", ANN.GENE)) -> genes
for ( gene in genes$ANN.GENE){
  sewage_samps.dt %>% filter(NUC %in% loi$nnuc) %>% group_by(NUC, LocationID) %>% mutate(n = n()) %>% filter(n > 1) %>% ungroup() %>% group_by(NUC) %>% mutate(n = length(unique(LocationID))) %>% filter(n > 1)  %>% filter(ANN.GENE == gene | (is.na(ANN.GENE) & "intergenic" == gene)) %>% ggplot(aes(x = as.Date(sample_date), y = value.freq, color = LocationID, group = LocationID)) + geom_point(size = 2, alpha = 0.66) + geom_smooth(formula = y ~ x, method = "glm", method.args = list(family = binomial(link='logit')), se = F) + facet_wrap(~paste(paste(ANN.GENE, ANN.AA, sep=":"),NUC, sep="\n"), scale = "free_y")  + theme_bw() + theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) + scale_x_date(date_breaks="2 weeks", date_labels="%d %b") + ylab("Mutation frequency") + xlab("Sampling Date") + ggtitle(paste("Selected Mutations of Interest:", gene)) + scale_color_viridis_d(name = "AA Mutation", option = "H") + guides(color=guide_legend(ncol=1)) -> p

  sewage_samps.dt %>% filter(NUC %in% loi$nnuc) %>% group_by(NUC, LocationID) %>% mutate(n = n()) %>% filter(n > 1) %>% ungroup() %>% group_by(NUC) %>% mutate(n = length(unique(LocationID))) %>% filter(n > 1)  %>% filter(ANN.GENE == gene | (is.na(ANN.GENE) & "intergenic" == gene)) %>% ungroup %>% select(LocationID) %>% distinct() %>% group_by(1) %>% summarize(n = ceiling(sqrt(n()))) -> plotsize
  ggsave(filename = paste0(opt$dir, "/EXPLORE/Mutations_of_Interest_2quantile_", gene, ".", opt$graph), plot = p, width = 3+plotsize$n*1.6, height = 3+plotsize$n*0.9)
}
