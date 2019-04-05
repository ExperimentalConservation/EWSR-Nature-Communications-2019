##run three different instances of the full model, with different seeding points

##set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('./From the start sims/1.1) From the start 1.R', chdir = TRUE)
source('./From the start sims/1.2) From the start 2.R', chdir = TRUE)
source('./From the start sims/1.3) From the start 3.R', chdir = TRUE)