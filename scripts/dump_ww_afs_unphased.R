library(readr)
library(stringr)
library(tidyr)
library(openxlsx)
library(odbc)
library(DBI)
library("optparse")
library(data.table)
library(dplyr)


option_list = list(
  make_option(c("--sample_date"), type="character", default="2020-06-01",
              help="earliest sample date timepoint [default= %default]"),
  make_option(c("--include_in_report"), type="character", default="NA,TRUE",
              help="samples with following include_in_report flag are considered [default= %default]"),
  make_option(c("--bsf_runs"), type="character", default="",
              help="comma sep. list of bsf runs, empty string for all [default= %default]"),
  make_option(c("--bsf_samples"), type="character", default="",
              help="comma sep. list of bsf sample ids, empty string for all [default= %default]"),
  make_option(c("--n_in_cons"), type="integer", default=20000,
              help="maximal Ns in consensus [default= %default]"),
  make_option(c("--minfreq"), type="double", default=0.001,
              help="minimal allele frequency [default= %default]"),
  make_option(c("-o", "--output_file_base"), type="character", default="sewagevariants_afs_dps",
              help="output name base [default= %default]"),
  make_option(c("--nmOnly"), type="logical", default=TRUE,
              help="Should only be samples from WWTP which were sampled after Sept 2022 at least once be considered [default= %default]"),
  make_option(c("--server"), type="character", default="undisclosed",
              help="Server name [default %default]", metavar="character"),
  make_option(c("--db"), type="character", default="undisclosed",
              help="Database name [default %default]", metavar="character"),
  make_option(c("--uid"), type="character", default="undisclosed",
              help="UID for DB [default %default]", metavar="character"),
  make_option(c("--pwd"), type="character", default="undisclosed",
              help="Password for DB [default %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# connect to data base
con <- dbConnect(odbc(),
  Driver = "ODBC Driver 17 for SQL Server",
  Server = opt$server,
  Database = opt$db,
  UID = opt$uid,
  PWD = opt$pwd
)



# get all WWTP which were sequenced after Sep 1 2022
if (opt$nmOnly){
  print("WARNING: only samples from WWTP were at least one sample was drawn after Sept 2022 are considered.")

  samples_db = setDT(dbReadTable(con, "samples"))
  seqdata_db = setDT(dbReadTable(con, "seqData"))

  colnames(seqdata_db) <- gsub("_coronA", "", colnames(seqdata_db))
  colnames(samples_db) <- gsub("_coronA", "", colnames(samples_db))

  left_join(x = samples_db, y = seqdata_db, by = "RNA_ID_int", suffix = c(".samples",".seqData")) -> samples_with_seqdata_db

  samples_db %>% filter(sample_date > as.Date("2022-09-01")) %>% select(LocationID) %>% distinct() -> wwtpOI

  samples_with_seqdata_db %>% filter(LocationID %in% wwtpOI$LocationID) %>% select(BSF_sample_name) %>% filter(!is.na(BSF_sample_name)) -> sampleOI
}

# create include in report query:
include_in_rep = strsplit(opt$include_in_report,",",fixed = TRUE)[[1]]
include_query = {}
if("NA" %in% include_in_rep){ include_query = c(include_query,"INCLUDE_IN_REPORT is null") }
if("TRUE" %in% include_in_rep){ include_query = c(include_query,"INCLUDE_IN_REPORT = 1 ") }
if("FALSE" %in% include_in_rep){ include_query = c(include_query,"INCLUDE_IN_REPORT = 2 ") }
if (length(include_query) > 0){
  include_query = paste("(", paste(include_query,collapse = " OR "),") AND ")
}

bsf_runs = strsplit(opt$bsf_runs,",",fixed = TRUE)[[1]]
bsf_runs_query  = ""
if (length(bsf_runs) > 0){
  bsf_runs_query = paste(" AND BSF_RUN in (", paste(sapply(bsf_runs, function(x) paste0("'",x,"'")), collapse = " , "), " ) ")
}


bsf_samples = strsplit(opt$bsf_samples,",",fixed = TRUE)[[1]]
bsf_samples_query  = ""
if (length(bsf_samples) > 0){
  bsf_samples_query = paste(" AND BSF_SAMPLE_NAME in (", paste(sapply(bsf_samples, function(x) paste0("'",x,"'")), collapse = " , "), " ) ")
}


sql = paste0("select BSF_sample_name, chrom, position, ref, alt, ann_gene, ann_aa, allele_freq, call_depth, phred_quality, hrun
from allele_freq as afs
inner join ( select ID, BSF_sample_name
             from seqData
             where ", include_query ,
             "N_in_Consensus is not null AND N_in_Consensus < ", opt$n_in_cons , bsf_runs_query, bsf_samples_query,
             "AND RNA_ID_int IN
                                        ( select RNA_ID_int
                                        from Samples
                                        where host = 'wastewater' AND sample_date > '", opt$sample_date,"' ) ) as ids
on afs.Sample_ID = ids.ID
left join ( select nuc_sub_ID, chrom, ref, alt
             from nuc_substitution
) as ns
on afs.nuc_substitution_ID = ns.nuc_sub_ID
left join ( select ann_ID, ann_gene, ann_featureid, ann_effect, ann_aa
             from annotations
) as anno
on afs.annotation_1_ID = anno.ann_ID
where allele_freq >= ", opt$minfreq)

ww_afs=setDT(dbGetQuery(con, sql))
colnames(ww_afs) = c("SAMPLEID", "CHROM", "POS", "REF", "ALT", "GENE", "AA", "AF", "DP", "PQ", "HRUN")

if(opt$nmOnly){
  ww_afs %>% filter(SAMPLEID %in% sampleOI$BSF_sample_name) -> ww_afs
  as.data.table(ww_afs) -> ww_afs
}

ww_afs[, AA := gsub("^p\\.","",AA)]
gzfh <- gzfile(paste0(opt$output_file_base,".tsv.gz"), "w")
write.table(ww_afs, file=gzfh, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
close(gzfh)

quit(save = "no")
