## Test data generation and subsetting

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(rjson))

# Find root file of the project and load functions
root <- rprojroot::find_rstudio_root_file()
source(file.path(root, "scripts/prepare_mae.R"))
source(file.path(root, "unit_tests/helpers.R"))

# Find datasets used in current pipeline run
current_datasets <- get_current_datasets()

test_that("cells in subsets are assigned the right class", {
  for (ds in current_datasets) {
    message(ds)
    config <- fromJSON(file = file.path(root, paste0("config/", ds, ".json")))
    subsets <- readRDS(file.path(root, config$subfile))
    data <- readRDS(file.path(root, config$mae))
    data <- clean_mae(mae = data, groupid = config$groupid)
    pdt <- pData(data)
    for (i in config$sizes) {
      for (j in nrow(subsets$out_condition[[as.character(i)]])) {
        datasub <- subset_mae(data, subsets$keep_samples, sz = i, i = j, 
                              imposed_condition = subsets$out_condition, filt = "")
        
        ## Check that the condition in "out_condition" is the same as that given
        ## in the phenodata slot of the data set for the same sample
        expect_equal(as.character(pdt[match(subsets$keep_samples[[as.character(i)]][j, ], 
                                            rownames(pdt)), paste(config$groupid, collapse = ".")]), 
                     as.character(gsub("\\.[1-2]$", "", subsets$out_condition[[as.character(i)]][j, ])))

        ## Check that all retained samples (in "keep_samples") belong to one of
        ## the groups that are intended to be kept
        expect_equal(as.character(sort(unique(pdt[match(subsets$keep_samples[[as.character(i)]][j, ], 
                                                        rownames(pdt)), 
                                                  paste(config$groupid, collapse = ".")]))), 
                     as.character(sort(config$keepgroups)))
        
        ## Check that the condition vector of the subsetted data set matches the
        ## "out_condition"
        expect_equivalent(datasub$condt[match(subsets$keep_samples[[as.character(i)]][j, ],
                                              names(datasub$condt))],
                          subsets$out_condition[[as.character(i)]][j, ])
        
        ## Check that the cells in the subsetted objects are all in the same
        ## order
        expect_equal(colnames(datasub$count), colnames(datasub$tpm))
        expect_equal(colnames(datasub$count), names(datasub$condt))
      }
    }
  }
})
