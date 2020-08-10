## Test COBRAData generation

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(iCOBRA))

# Find root file of the project and load functions
root <- rprojroot::find_rstudio_root_file()
source(file.path(root, "scripts/prepare_mae.R"))
source(file.path(root, "unit_tests/helpers.R"))

# Find datasets and filters used in current pipeline run
current_datasets <- get_current_datasets()
current_filters <- get_current_filters()

get_method <- function(x) sapply(strsplit(x, "\\."), .subset, 1)
get_nsamples <- function(x) sapply(strsplit(x, "\\."), .subset, 2)
get_repl <- function(x) sapply(strsplit(x, "\\."), .subset, 3)

test_that("COBRAData object is correctly assembled", {
  for (ds in current_datasets) {
    config <- fromJSON(file = file.path(root, paste0("config/", ds, ".json")))
    subsets <- readRDS(file.path(root, config$subfile))
    data <- readRDS(file.path(root, config$mae))
    data <- clean_mae(mae = data, groupid = config$groupid)
    
    for (f in current_filters) {
      if (f != "") f <- paste0("_", f)
      cobra <- readRDS(
        file.path(root, "output/cobra_data", paste0(ds, f, "_cobra.rds"))
      )
      ngenes <- readRDS(
        file.path(root, "output/cobra_data", paste0(ds, f, "_nbr_called.rds"))
      )
      
      all_methods <- unique(get_method(colnames(padj(cobra))))
      for (mth in all_methods) {
        ## Test that adjusted p-values in COBRAData object are the same as in
        ## the individual result files
        res <- readRDS(file.path(root, "results", paste0(ds, "_", mth, ".rds")))
        for (nm in names(res)) {
          if ("padj" %in% colnames(res[[nm]]$df)) {
            expect_equal(res[[nm]]$df$padj[match(rownames(padj(cobra)), rownames(res[[nm]]$df))], 
                         padj(cobra)[, nm])
          }
          datasub <- subset_mae(data, subsets$keep_samples, sz = get_nsamples(nm), 
                                i = as.numeric(as.character(get_repl(nm))), 
                                imposed_condition = subsets$out_condition, 
                                filt = gsub("^_", "", f))
          ## Test that the calculated number of tested, called and significant
          ## genes are correct
          if (nm %in% colnames(padj(cobra))) {
            ntested <- nrow(datasub$count)
            ncalled <- sum(!is.na(padj(cobra)[, nm]))
            nsign <- sum(padj(cobra)[, nm] <= 0.05, na.rm = TRUE)
            tmp <- subset(ngenes, method == mth & ncells == get_nsamples(nm) & 
                            repl == get_repl(nm) & dataset == ds & 
                            filt == gsub("^_", "", f))
            expect_equal(ntested, tmp$nbr_tested)
            expect_equal(ncalled, tmp$nbr_called)
            expect_equal(nsign, tmp$nbr_sign_adjp0.05)
          }
        }
        
        ## Test that truth is equivalent to adjusted p-values from largest
        ## sample set
        maxn <- max(as.numeric(as.character(get_nsamples(names(res)))))
        if (paste0(mth, ".", maxn, ".1") %in% colnames(padj(cobra))) {
          tmp <- padj(cobra)[, paste0(mth, ".", maxn, ".1")]
          tmp[is.na(tmp)] <- 1
          expect_equal(as.numeric(tmp <= 0.05),
                       truth(cobra)[match(rownames(padj(cobra)), rownames(truth(cobra))), 
                                    paste0(mth, ".truth")])
        }
      }
    }        
  }
})
