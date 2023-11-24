## load libraries
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("lubridate"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggridges"))


# get Options
option_list = list(
  make_option(c("--data"), type="character", default="output-variants/globalFullData.csv.gz",
              help="Path to file containing lineage and mutation freqs, output of VaQuERo/scripts/VaQuERo_v2.r [default %default]", metavar="character"),
  make_option(c("--label"), type="character", default="LABEL",
              help="String to used in filenames [default %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###
  #opt$data = "globalFullData.csv.gz"
##


# define functions
## calculate for beta distribution beta (shape2) from alpha (shape1) and mean
sa2b <- function(s, a) {
   #b <- (s*a-a)/(-1*s)
   b <- (a*(1-s))/s
   return(b)
}

## calculate for beta distribution theoretical variance from alpha (shape1) and beta (shape2)
ab2var <- function(a, b){
  var <- (a * b)/( ((a + b)^2)*(a + b + 1))
  return(var)
}

## objective function to find alpha (shape1) which produces with given mean (hence beta or shape2) observed variance closest to theoretical variance
meanVarResidual <- function(data, par) {
  means    <- data[,1]
  vars     <- data[,2]
  alpha    <- par[1]
  betasEst <- (means*alpha-alpha)/(-1*means)
  varsEst  <- (alpha * betasEst)/(((alpha + betasEst)^2)*(alpha + betasEst + 1))
  residual <- sum( ((varsEst - vars)^2) )
  return(residual)
}


## bin data (in an overlapping fashion) and sample (with replacement from each bin equal number of data points)
## take a table with columns "lineage" and "mutation", specifying the inferred lineage frequency and the observed mutation freq
## produces alpha*, which allows to calculate for each mean a beta value which produces a similar variance as the observed data
## alpha* is chosen in a way that it minimizes the sum of square roots between theoretical and observed variance in 20 bins (if populated) across the (0,1) intervall
## additionally, a few diagnostics are shared too
mutationfreq2alphaprime <- function(dt_, binsize = 1/20){
    # dt_= adt; binsize = 1/bincount
    binnings <- seq(from=0,to=1,by=binsize)

    dt2_ <- data.table(lineage = numeric(), mutation = numeric(), bin=character())
    #dt_$bin <- cut(dt_$lineage, breaks = binnings)
    binoverlap <- round((1/binsize)/20)
    print("")
    print(paste("Averaging over", binoverlap, "out of", 1/binsize, "bins"))
    print("")
    for (bin in 1:length(binnings)) {
      low_bin <- ifelse( (bin-floor(binoverlap/2))>=1, bin-floor(binoverlap/2), 1 )
      up_bin  <- ifelse( (bin+ceiling(binoverlap/2))<=length(binnings), bin+ceiling(binoverlap/2), length(binnings) )
      low_limit <- binnings[low_bin]
      up_limit <- binnings[up_bin]
      dt2_ <- dt_ %>% filter(lineage >= low_limit & lineage < up_limit) %>%
                mutate(bin=paste0("(", low_limit, ",", up_limit, "]")) %>%
                rbind(dt2_)
    }
    dt_ <- dt2_
    rm(dt2_)

    ## remove outlier per group flag
    ## remove if outside iqr_margin*IQR +/- 0.05
    iqr_margin <- 1.5
    print(paste("LOG: number of mutations before outlier removal", dim(dt_)[1]))
    dt_ <- dt_ %>%
            group_by(bin) %>%
            mutate(iqr = IQR(mutation)) %>%
            mutate(upperbond = quantile(mutation, 0.75) + iqr_margin * iqr, lowerbond = quantile(mutation, 0.25) - iqr_margin * iqr) %>%
            filter(mutation <= (upperbond+0.05) & mutation >= (lowerbond-0.05)) %>%
            ungroup() %>% dplyr::select(-"iqr", -"upperbond", -"lowerbond")
    print(paste("LOG: number of mutations after outlier removal", dim(dt_)[1]))

    dt_ <- dt_ %>%
            group_by(bin) %>%
            mutate(n = n()) %>%
            filter(n>3) %>%
            dplyr::select(-"n") %>%
            mutate(varianceperbin = var(mutation))
    odt_ <- dt_ %>%
            group_by(bin) %>%
            sample_n(size = 1000, replace = TRUE) %>%
            ungroup() %>%
            dplyr::select("lineage", "varianceperbin")

    optim_output <- optim(par = c(2),
                          fn = meanVarResidual,
                          data = odt_,
                          method = "Brent",
                          lower = 0.001,
                          upper = 1000)

    vdt_ <- dt_ %>% ungroup() %>% dplyr::select(lineage, varianceperbin)
    vdt_ <- vdt_ %>% distinct()
    vdt_$mode <- "observed"
    vdt_$beta = NA
    for (s in binnings){
      beta = sa2b(s, optim_output$par)
      vals <- rbeta(n = 1000, shape1 = optim_output$par, shape2 = beta)
      vdt_ <- rbind(vdt_, data.table(lineage = s, varianceperbin = var(vals), mode = "theoretical", beta = beta))
    }
    vdt_ <- vdt_ %>% group_by(varianceperbin, mode, beta) %>% summarise(lineage = mean(lineage))

    ldt_ <- vdt_ %>% ungroup() %>% filter(!is.na(beta)) %>% dplyr::select(lineage, beta) %>% distinct()
    normalizer <- max(vdt_$varianceperbin) / max(dt_$mutation)

    bin2color <- dt_ %>% select(bin) %>% distinct() %>% arrange(bin) %>% mutate(idx = cur_group_id()) %>% mutate(fill = case_when(idx %% 3 == 0 ~ "grey25", idx %% 3 == 1 ~ "grey50", idx %% 3 == 2 ~ "grey75"))
    left_join(x = dt_, y = bin2color, by = join_by(bin)) -> dt_


    p <- ggplot() +
         geom_point(data = dt_, aes(x = lineage, y = mutation, fill = fill), color = "white", shape = 21, size = 1.3, alpha = 0.3) +
         geom_smooth(data = vdt_, aes(x = lineage, y = varianceperbin/normalizer, color = mode), alpha = .6, linewidth = 2, se = FALSE) +
         #geom_point(data = vdt_, aes(x = lineage, y = varianceperbin/normalizer, color = mode, shape = mode), alpha = .6, size = 3) +
         theme_bw() +
         scale_y_continuous(name = "Observed Mutation Freq [points]", sec.axis = sec_axis(trans= ~.*normalizer, name = 'Variance [lines]'), limits = c(-0.05,1)) +
         scale_color_brewer(name = "Variance\ncalculation\nmethod", palette = "Set1") +
         scale_fill_identity(name = "Lineage\nFreq. Bin") +
         ggtitle(paste("Deduced Alpha*", round(optim_output$par, digits = 4))) +
         geom_point(data = ldt_, aes(x = lineage, y = -0.05, size = beta), alpha = 0.1) +
         scale_size(guide = "legend", name = "Beta per bin") +
         coord_fixed()

    q <- ggplot(data = dt_, aes(x = mutation, y = bin, fill = after_stat(quantile))) +
         stat_density_ridges(lwd = 0, linewidth = 0, size = 0, rel_min_height = 0.01, scale = 5, alpha = 0.3, quantile_lines = TRUE, calc_ecdf = TRUE, geom = "density_ridges_gradient", quantiles = c(0.05, 0.4, 0.6, 0.95)) +
         theme_bw() +
         ylab("binned lineage frequency") +
         xlab("mutation frequency") +
         scale_fill_manual(name = "Prob.", values = c("#8BD3E6", "#FF6D6A", "#E9EC6B", "#EFBE7D", "#B1A2CA"), labels = c("(0, 5%]", "(5, 40%]", "(40, 60%]","(60%, 95%]", "(95%, 1]")) +
         theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

    return(list(alphaprime = optim_output$par, plot1 = p, plot2 = q))
}

## read in file
dt <- fread(opt$data)

## calculate median per group, if group > 4, and use this as lineage freqs
if(FALSE){
  dt <- dt %>%
          group_by(all_variants, LocationID, sample_date) %>%
          mutate(n = n()) %>%
          filter(n>4) %>%
          mutate(value_model = value, value = median(singlevalue))
  ggplot(data = dt, aes(x = value_model, y = value)) + geom_point() + geom_smooth(method = "lm") + geom_abline(color = "red")
}

# collapse entries per variant
# remove none uniq markers
# remove fixed mutations
# bring in required shape
adt <- dt %>% dplyr::select(-"variant") %>%
          distinct() %>%
          filter(!grepl(";", all_variants)) %>%
          dplyr::select(value, singlevalue) %>%
          rename("lineage" = "value", "mutation" = "singlevalue")

bincount <- 500 # should be larger 11
adt <- adt %>% filter(mutation < 1-1/(10*bincount))
alphaprime <- mutationfreq2alphaprime(dt_ = adt, binsize = 1/bincount)
print(paste("RESULT: ", alphaprime$alphaprime))
#alphaprime$plot1


ggsave(filename=paste0("diagnostics_mutationFreq_distro_", opt$label, ".pdf"), plot = alphaprime$plot2, width = 6, height = 7)

ggsave(filename=paste0("diagnostics_variance_compare_", opt$label, ".pdf"), plot = alphaprime$plot1, width = 7, height = 6)

ggsave(filename=paste0("diagnostics_mutationFreq_distro_", opt$label, ".png"), plot = alphaprime$plot2, width = 7, height = 7)

ggsave(filename=paste0("diagnostics_variance_compare_", opt$label, ".png"), plot = alphaprime$plot1, width = 7, height = 6)
