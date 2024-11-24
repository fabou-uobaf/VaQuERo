## installs the dependancies w/ respect to R packages for VaQuERo and VaQera

inst_packages <-  installed.packages()

need_packages <- c("tidyr", "ggplot2", "reshape2", "dplyr", "data.table", "gamlss", "ggmap", "tmaptools", "ggrepel", "scales", "betareg", "ggspatial", "sf", "rnaturalearth", "rnaturalearthdata", "optparse", "stringr", "lubridate", "rjson", "RColorBrewer", "ggsankey", "cowplot", "viridis", "devtools", "patchwork", "ggpubr", "dendextend", "circlize", "NbClust", "gslnls")
need_dev_packages <- c("davidsjoberg/ggsankey")
need_bioc_packages <- c("ComplexHeatmap")

counter = 0

need_to_install_packages <- need_packages[!need_packages %in% inst_packages]
for (need_to_install_package in  need_to_install_packages){
  print(paste("installing", need_to_install_package, "from http://cran.us.r-project.org"))
  install.packages(need_to_install_package, repos="http://cran.us.r-project.org")
  counter <- counter+1
}


for (maybe_need_to_install_dev_package in  need_dev_packages){
  split_package_name <- unlist(strsplit(maybe_need_to_install_dev_package, split ="/"))
  if(!split_package_name[2] %in% inst_packages){
    print(paste("installing", maybe_need_to_install_dev_package, "from github"))
    devtools::install_github(maybe_need_to_install_dev_package)
    counter <- counter+1
  }
}

need_to_install_bioc_packages <- need_bioc_packages[!need_bioc_packages %in% inst_packages]
for (need_to_install_bioc_package in  need_to_install_bioc_packages){
  print(paste("installing", need_to_install_bioc_package, "from BiocManager"))
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install(need_to_install_bioc_package)
  counter <- counter+1
}


print(paste("LOG: all set;", counter, "packages installed"))
