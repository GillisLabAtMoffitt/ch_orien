# Prep overall ORIEN data
# Import Library
library(tidyverse)

################################################################################# I ### Load data
################################################################################# 2023----
# Load each file in its own dataframe
# filenames <- list.files(path = paste0(path, "/23PRJ127MCC_NormalizedFiles"),
#                         pattern = "23PRJ127MCC_20230620_.*.csv")
# names <- str_match(filenames, "23PRJ127MCC_20230620_(.*)_V4.csv")[,2]
#
# for(i in names){
#     filepath <- file.path(paste0(path, "/23PRJ127MCC_NormalizedFiles/",paste("23PRJ127MCC_20230620_", i, "_V4.csv",sep="")))
#     assign(i, read.csv(filepath, na.strings = ""))
# }
#
# rm(filenames, names, filepath, i)
# save.image(file = "initial_files.RData)

################################################################################# 2024----
# # Load each file in its own dataframe
# path_raw <- fs::path("", "Volumes", "Gillis_Research",
#                      "Lab_Data", "CHinORIEN")
# filenames <- list.files(path = paste0(path_raw,
#                                       "/RawData",
#                                       "/23PRJ127MCC_NormalizedFiles_20240715"),
#                         pattern = "23PRJ127MCC_20240715_.*.csv")
# names <- str_match(filenames, "23PRJ127MCC_20240715_(.*)_V4.csv")[,2]
# 
# for(i in names){
#     filepath <- file.path(paste0(path_raw,
#                                  "/RawData",
#                                  "/23PRJ127MCC_NormalizedFiles_20240715/",
#                                  paste("23PRJ127MCC_20240715_", i, "_V4.csv",sep="")))
#     assign(i, read.csv(filepath, na.strings = ""))
# }
# 
# rm(filenames, names, filepath, i)
# # Saved environment image in a .Rdata
# save.image(file = paste0(path_raw,
#                          "/RawData",
#                          "/23PRJ127MCC_NormalizedFiles_20240715/",
#                          "CHinORIEN_AllTumorTypes_Files.RData"))


################################################################################# 2026----
# Load each file in its own dataframe
path_raw <- fs::path("", "Volumes", "Gillis_Research",
                     "Lab_Data", "CHinORIEN")
filenames <- list.files(path = paste0(path_raw,
                                      "/RawData",
                                      "/23PRJ127MCC_NormalizedFiles_20260730"),
                        pattern = "23PRJ127MCC_20260303_.*.csv")
names <- str_match(filenames, "23PRJ127MCC_20260303_(.*)_V4.csv")[,2]

for(i in names){
    filepath <- file.path(paste0(path_raw,
                                 "/RawData",
                                 "/23PRJ127MCC_NormalizedFiles_20260730/",
                                 paste("23PRJ127MCC_20260303_", i, "_V4.csv",sep="")))
    assign(i, read.csv(filepath, na.strings = ""))
}

rm(filenames, names, filepath, i)
# Saved environment image in a .Rdata
save.image(file = paste0(path_raw,
                         "/RawData",
                         "/23PRJ127MCC_NormalizedFiles_20260730/",
                         "CHinORIEN_AllTumorTypes_Files_20260730.RData"))




