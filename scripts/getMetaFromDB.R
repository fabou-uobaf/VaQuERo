  suppressPackageStartupMessages(library(odbc))
  suppressPackageStartupMessages(library(DBI))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(optparse))


  option_list = list(
  make_option(c("--cat"), type="character", default="general", 
              help="Report categories to be considered [default %default]", metavar="character"),
  make_option(c("--server"), type="character", default="undisclosed", 
              help="Server name [default %default]", metavar="character"),
  make_option(c("--db"), type="character", default="undisclosed", 
              help="Database name [default %default]", metavar="character"),
  make_option(c("--uid"), type="character", default="undisclosed", 
              help="UID for DB [default %default]", metavar="character"),
  make_option(c("--pwd"), type="character", default="undisclosed", 
              help="Password for DB [default %default]", metavar="character"),
  make_option(c("--out"), type="character", default="metaData_general.csv", 
              help="path to outputfile [default %default]", metavar="character")
  );
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  
  # samples considered if asociated with any of the following tags
  reportcategories <- unlist(str_split(opt$cat, pattern=c(";",",")))

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

  left_join(x = samples_db, y = seqdata_db, by = "RNA_ID_int", suffix = c(".samples",".seqData")) -> samples_with_seqdata_db
  
  left_join(x = samples_with_seqdata_db, y = sewagePlants_db, by = "LocationID", suffix = c(".db",".location")) -> samples_with_seqdata_with_wwplants_db
  
  samples_with_seqdata_with_wwplants_db %>% filter(host == "wastewater")  %>% filter(include_in_report == TRUE | is.na(include_in_report)) %>% mutate(include_in_report = TRUE) %>% filter(!is.na(sample_date)) -> samples_with_seqdata_with_wwplants_db
  
  samples_with_seqdata_with_wwplants_db  %>% filter(grepl(paste(reportcategories,collapse="|"), report_category)) -> samples_with_seqdata_with_wwplants_db
  
  ## replace unset but essantial fields
  samples_with_seqdata_with_wwplants_db %>% mutate(adress_town = ifelse(is.na(adress_town) & !is.na(connected_towns), connected_towns, adress_town)) -> samples_with_seqdata_with_wwplants_db
  samples_with_seqdata_with_wwplants_db %>% mutate(adress_town = ifelse(is.na(adress_town) & !is.na(LocationName), LocationName, adress_town)) -> samples_with_seqdata_with_wwplants_db
  samples_with_seqdata_with_wwplants_db %>% mutate(adress_town = ifelse(is.na(adress_town), "unknown", adress_town)) -> samples_with_seqdata_with_wwplants_db
  
  samples_with_seqdata_with_wwplants_db %>% mutate(connected_people = ifelse(is.na(connected_people), 1234, connected_people)) -> samples_with_seqdata_with_wwplants_db
  
  samples_with_seqdata_with_wwplants_db %>% mutate(dcpLatitude = ifelse(is.na(dcpLatitude), 0, dcpLatitude)) -> samples_with_seqdata_with_wwplants_db
  
  samples_with_seqdata_with_wwplants_db %>% mutate(dcpLongitude = ifelse(is.na(dcpLongitude), 0, dcpLongitude)) -> samples_with_seqdata_with_wwplants_db
  
  ## filter for data set fullfiling all criteria
  samples_with_seqdata_with_wwplants_db %>% filter(! is.na(BSF_run) & ! is.na(BSF_sample_name) & ! is.na(BSF_start_date) & ! is.na(LocationID) & ! is.na(LocationName) & ! is.na(N_in_Consensus) & ! is.na(RNA_ID_int) & ! is.na(adress_town) & ! is.na(connected_people) & ! is.na(dcpLatitude) & ! is.na(dcpLongitude) & ! is.na(include_in_report) & ! is.na(report_category) & ! is.na(sample_date) & ! is.na(status)) -> samples_with_seqdata_with_wwplants_db
  
  coi <- c("BSF_run", "BSF_sample_name", "BSF_start_date", "LocationID", "LocationName", "N_in_Consensus", "RNA_ID_int", "additional_information", "adress_town", "connected_people", "dcpLatitude", "dcpLongitude", "include_in_report", "report_category", "sample_date", "status")

  samples_with_seqdata_with_wwplants_db %>% select(all_of(coi)) -> metaData  
  fwrite(metaData, file = opt$out, sep = "\t", col.names = T, row.names = F)

