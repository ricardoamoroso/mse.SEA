## Multispecies MSE

### Run online
https://ricardoamoroso.shinyapps.io/mseSEA/

### Run locally
Install R, then run:

install.packages(c("shiny", "shinyjs", "ggplot2", "dplyr",
                   "tidyr", "data.table", "ggpubr", "zip",
                   "scales", "DT"))

shiny::runGitHub("mse.SEA", "ricardoamoroso", ref = "master")
