  suppressPackageStartupMessages(library(odbc))
  suppressPackageStartupMessages(library(DBI))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(optparse))


  option_list = list(
  make_option(c("--cat"), type="character", default="general",
              help="Report categories to be considered [default %default]", metavar="character"),
  make_option(c("--samples"), type="character", default=NA,
              help="List of samples to be considered [default %default]", metavar="character"),
  make_option(c("--server"), type="character", default="undisclosed",
              help="Server name [default %default]", metavar="character"),
  make_option(c("--db"), type="character", default="undisclosed",
              help="Database name [default %default]", metavar="character"),
  make_option(c("--uid"), type="character", default="undisclosed",
              help="UID for DB [default %default]", metavar="character"),
  make_option(c("--pwd"), type="character", default="undisclosed",
              help="Password for DB [default %default]", metavar="character"),
  make_option(c("--out"), type="character", default="metaData_general.csv",
              help="path to outputfile [default %default]", metavar="character"),
  make_option(c("--sampleAfter"), type="character", default="2020-01-01",
              help="Report samples with sample date after specified value [default %default]", metavar="character"),
  make_option(c("--wwtpAfter"), type="character", default="2020-01-01",
              help="Report all samples from WWTP with at least one sample from after specified value [default %default]", metavar="character")
  );
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);

  # samples considered if associated with any of the following tags
  reportcategories <- unique(unlist(str_split(opt$cat, pattern=c(";",","))))

  # samples considered by sample ID
  reportsamples <- unique(unlist(str_split(opt$samples, pattern=c(";", ","))))
  reportsamples <- reportsamples[!grepl(paste(c(";", ","), collapse="|"), reportsamples)]

  # set time limit from when samples shall be considered
  sampleStart = opt$sampleAfter
  # set time limit from when wwtp shall be considered
  wwtpAfter = opt$wwtpAfter

  # connect to data base
  con <- dbConnect(odbc(),
    Driver = "ODBC Driver 17 for SQL Server",
    Server = opt$server,
    Database = opt$db,
    UID = opt$uid,
    PWD = opt$pwd
  )


  # get tables
  samples_db = setDT(dbReadTable(con, "samples"))
  seqdata_db = setDT(dbReadTable(con, "seqData"))
  sewagePlants_db = setDT(dbReadTable(con, "sewagePlants"))

  ## make colnames more generic
  colnames(seqdata_db) <- gsub("_coronA", "", colnames(seqdata_db))
  colnames(sewagePlants_db) <- gsub("_coronA", "", colnames(sewagePlants_db))
  colnames(samples_db) <- gsub("_coronA", "", colnames(samples_db))

  left_join(x = samples_db, y = seqdata_db, by = "RNA_ID_int", suffix = c(".samples",".seqData"), multiple = "all") -> samples_with_seqdata_db

  left_join(x = samples_with_seqdata_db, y = sewagePlants_db, by = "LocationID", suffix = c(".db",".location")) -> samples_with_seqdata_with_wwplants_db

  if (all(is.na(reportsamples))){
    samples_with_seqdata_with_wwplants_db %>% filter(host == "wastewater")  %>% filter(include_in_report == TRUE | is.na(include_in_report)) %>% mutate(include_in_report = TRUE) %>% filter(!is.na(sample_date)) -> samples_with_seqdata_with_wwplants_db

    samples_with_seqdata_with_wwplants_db  %>% filter(grepl(paste(reportcategories,collapse="|"), report_category)) -> samples_with_seqdata_with_wwplants_db
  } else{
    samples_with_seqdata_with_wwplants_db %>% mutate(include_in_report = TRUE) %>% filter(!is.na(sample_date)) -> samples_with_seqdata_with_wwplants_db

    samples_with_seqdata_with_wwplants_db  %>% filter(RNA_ID_int %in% reportsamples | BSF_sample_name %in% reportsamples ) -> samples_with_seqdata_with_wwplants_db

  }

  ## replace unset but essential fields
  samples_with_seqdata_with_wwplants_db %>% mutate(adress_town = ifelse(is.na(adress_town) & !is.na(connected_towns), connected_towns, adress_town)) -> samples_with_seqdata_with_wwplants_db
  samples_with_seqdata_with_wwplants_db %>% mutate(adress_town = ifelse(is.na(adress_town) & !is.na(LocationName), LocationName, adress_town)) -> samples_with_seqdata_with_wwplants_db
  samples_with_seqdata_with_wwplants_db %>% mutate(adress_town = ifelse(is.na(adress_town), "unknown", adress_town)) -> samples_with_seqdata_with_wwplants_db

  samples_with_seqdata_with_wwplants_db %>% mutate(connected_people = ifelse(is.na(connected_people), 0, connected_people)) -> samples_with_seqdata_with_wwplants_db

  samples_with_seqdata_with_wwplants_db %>% mutate(dcpLatitude = ifelse(is.na(dcpLatitude), 0, dcpLatitude)) -> samples_with_seqdata_with_wwplants_db
  samples_with_seqdata_with_wwplants_db %>% mutate(dcpLongitude = ifelse(is.na(dcpLongitude), 0, dcpLongitude)) -> samples_with_seqdata_with_wwplants_db

  samples_with_seqdata_with_wwplants_db %>% mutate(state = ifelse(is.na(state), "Other", state)) -> samples_with_seqdata_with_wwplants_db


  ## filter for data set fulfilling all criteria
  samples_with_seqdata_with_wwplants_db %>% filter(! is.na(BSF_run) & ! is.na(BSF_sample_name) & ! is.na(BSF_start_date) & ! is.na(LocationID) & ! is.na(LocationName) & ! is.na(N_in_Consensus) & ! is.na(RNA_ID_int) & ! is.na(adress_town) & ! is.na(connected_people) & ! is.na(dcpLatitude) & ! is.na(dcpLongitude) & ! is.na(include_in_report) & ! is.na(report_category) & ! is.na(sample_date) & ! is.na(status)) -> samples_with_seqdata_with_wwplants_db

  ## set status to pass if regex fit with pass (e.g. passed_qc; in former pango version)
  ## remove all none-passed samples
  samples_with_seqdata_with_wwplants_db %>% mutate(status = ifelse(grepl("^pass", status), "pass", status)) -> samples_with_seqdata_with_wwplants_db
  #samples_with_seqdata_with_wwplants_db %>% filter(status == "pass") -> samples_with_seqdata_with_wwplants_db

  ## select columns of interest
  coi <- c("BSF_run", "BSF_sample_name", "BSF_start_date", "LocationID", "LocationName", "N_in_Consensus", "RNA_ID_int", "additional_information", "adress_town", "state", "connected_people", "dcpLatitude", "dcpLongitude", "include_in_report", "report_category", "sample_date", "status")
  samples_with_seqdata_with_wwplants_db %>% select(all_of(coi)) -> metaData

  ## filter samples after sampleStart
  metaData %>% filter(as.Date(sample_date) >= as.Date(sampleStart)) -> metaData

  ## filter samples from wwtp with at least one sample after wwtpAfter
  metaData %>% group_by(LocationID) %>% mutate(latest = max(sample_date, na.rm = TRUE)) %>% rowwise() %>% filter(latest >= as.Date(wwtpAfter)) %>% dplyr::select(-latest) -> metaData

  ## write data to file
  fwrite(metaData, file = opt$out, sep = "\t", col.names = T, row.names = F)
