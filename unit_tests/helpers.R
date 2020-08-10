# Find results directory
find_res_dir <- function() {
  root <- rprojroot::find_rstudio_root_file()
  res_dir <- file.path(root, "results")
  if (!dir.exists(res_dir)) {
    stop("No directory called `results` found.",
         "Make sure to first run the pipeline using `make`",
         call. = FALSE)
  }
  res_dir
}

# Find datasets used in current pipeline run
get_current_datasets <- function() {
  res_dir <- find_res_dir()
  res_files <- list.files(res_dir)
  
  # Extract and return current datasets
  unique(sub(pattern = "_.+\\.rds$", "", res_files))
}

# Find methods used in current pipeline run
get_current_methods <- function() {
  res_dir <- find_res_dir()
  res_files <- list.files(res_dir)
  
  # Extract and return current methods
  unique(sub("^[^_]+_([^_,.]+).*\\.rds", "\\1", res_files))
}

# Function to read in pipeline results
get_results <- function(res_dir = find_res_dir()) {
  res_files <- list.files(res_dir)
  names(res_files) <- gsub("\\.rds$", "", res_files)
  
  # Read in results
  lapply(file.path(res_dir, res_files), readRDS)
}
