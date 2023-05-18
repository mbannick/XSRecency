## code to prepare `cephia` dataset goes here

cephia <- fread("data-raw/cephia_public_use_dataset_20210604.csv")
usethis::use_data(cephia, overwrite=TRUE, internal=FALSE)
