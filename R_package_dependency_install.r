## istalls the dependancies w/ respect to R packages for VaQuERo

inst_packages <-  installed.packages()

need_packages <- c("tidyr", "ggplot2", "reshape2", "dplyr", "data.table", "gamlss", "ggmap", "tmaptools", "ggrepel", "scales", "betareg", "ggspatial", "sf", "rnaturalearth", "rnaturalearthdata", "optparse", "stringr", "lubridate")

need_to_install_packages <- need_packages[!need_packages %in% inst_packages]

for (need_to_install_package in  need_to_install_packages){
  print(need_to_install_package)
  install.packages(need_to_install_package, repos="http://cran.us.r-project.org")
}
