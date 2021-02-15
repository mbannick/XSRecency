
library(data.table)
library(magrittr)

# Which version to read in?
dir <- "~/Documents/FileZilla/xs-recent/14-02-21-Feb/"
f <- list.files(dir, full.names=T)
df <- lapply(f, fread) %>% rbindlist(fill=T)

