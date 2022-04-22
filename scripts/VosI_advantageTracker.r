  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(optparse)

  shorty <- function(x){
    n = 17
    if(nchar(x)>n){
      x <- paste0(substr(x, start = 0, stop = n), "...")
    }
    return(x)
  }

  option_list = list(
    make_option(c("--sannoFH"), type="character", default="resources/mutations_special.csv", 
              help="Path to special mutation definition file [default %default]", metavar="character"),
    make_option(c("--filebase"), type="character", default="specialmutations", 
              help="Identifier included in outputfile name [default %default]", metavar="character"),
    make_option(c("--vvariant"), type="character", default="BA.2", 
              help="Background Lineage to be examined [default %default]", metavar="character"),
    make_option(c("--sdtFH"), type="character", default="output-general/globalSpecialmutData.csv", 
              help="Path to special mutation quantification file [default %default]", metavar="character"),
    make_option(c("--vdtFH"), type="character", default="output-general/globalFittedData.csv", 
              help="Path to variant quantification file [default %default]", metavar="character")
  );
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);


  
  sanno <- fread(opt$sannoFH)
  sanno %>% mutate(mmarker = paste(NUC, paste(Gene,AA, sep=":"),paste0("[",Variants,"]"), sep = "|")) -> sanno
  filebase <- opt$filebase 
  vvariant  <- opt$vvariant
  sanno %>% filter(grepl(vvariant, mmarker)) -> sanno_2



  llocation <- c()
  #llocation <- c("ATTP_9-HKA-Simmering")
    
  sdt <- fread(opt$sdtFH)
  vdt <- fread(opt$vdtFH)
  
  if(length(llocation) > 0){
    sdt %>% filter( LocationID %in% llocation) %>% filter(marker %in% sanno_2$mmarker) -> sdt
    vdt %>% filter( LocationID  %in% llocation) %>% filter(variant == vvariant) -> vdt
  } else{
    sdt %>%  filter(marker %in% sanno_2$mmarker) -> sdt
    vdt %>%  filter(variant == vvariant) -> vdt
  }
  
  left_join(x = vdt, y = sdt, by = c("LocationID", "LocationName", "sample_date"), suffix = c(".variant", ".marker")) -> dt

  ## write table
  fwrite(x = dt, file = paste0("Markerproportion_", filebase, "_", gsub("/","_", vvariant), ".csv"), sep = "\t", col.names=TRUE, row.names = FALSE)

  
  dt %>% filter(value.variant > 0) %>% mutate(marker.proportion = value.marker/value.variant) -> dt

  ## remove WWTp with less then four data points
  dt %>% group_by(LocationName)  %>% mutate(N = length(unique(sample_date))) %>% filter(N > 3) -> dt
  
  ## remove marker which 98% quantile is below 10%
  dt %>% group_by(marker) %>% mutate(max.marker = quantile(value.marker, probs = 0.98)) %>% filter(max.marker > 0.1) -> dt

  if(dim(dt)[1] > 1){  
    ## shorten name
    dt %>% rowwise() %>% mutate(LocationName = shorty(LocationName)) -> dt
  
    ## make plot      
    facetSquares <- ceiling(sqrt(length(unique(c(dt$LocationName)))))
    facetSquares_c <- min(facetSquares, 10)
    facetSquares_r <- ceiling(length(unique(c(dt$LocationName)))/facetSquares_c)
      


    ggplot(data = dt, aes(x = as.Date(sample_date))) + geom_line(aes(y = value.variant, linetype = variant.variant), alpha = 0.5, ) + geom_smooth(aes(y = marker.proportion, color = marker), span = 2 , se=FALSE) + geom_point(aes(y = marker.proportion, shape = marker), size = 1, alpha = 0.5) + theme_bw() + xlab("Sample Date") + ylab("Relative abundance [1/1]") + scale_x_date(date_breaks = "2 weeks", date_labels =  "%Y-%m-%d") + facet_wrap(~LocationName, nrow = facetSquares_r, ncol =facetSquares_c) +  theme(legend.position="bottom") + scale_shape(name = "Sub-Lineage Marker") + scale_color_discrete(name = "Sub-Lineage Marker") + scale_linetype(name = "Background Lineage") + scale_y_continuous(limits = c(0,1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))-> p
  
    ggsave(filename = paste0("Markerproportion_", filebase, "_", gsub("/","_", vvariant), ".pdf"), plot = p, width = 4/3*facetSquares_c, heigh = 0.5+1.2*facetSquares_r)
  }
