## Bias, robustness and scalability in differential expression analysis of single-cell RNA-seq data

This repository contains all the necessary code to perform the evaluation of differential expression analysis methods in single-cell RNA-seq data, available in 

* C Soneson & MD Robinson: [Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612). Nature Methods 15:255-261 (2018).

In this paper, we compare the performance of more than 30 approaches to differential gene expression analysis in the context of single-cell RNA-seq data. The main results can be further browsed in a [shiny app](http://imlspenticton.uzh.ch:3838/scrnaseq_de_evaluation). 

**Note:** The purpose of the `conquer_comparison` repository is to provide a public record of the exact code that was used for our publication ([Soneson & Robinson, Nature Methods 2018](https://www.nature.com/articles/nmeth.4612)). In particular, it is not intended to be a software package or a general pipeline for differential expression analysis of single-cell data. As a consequence, running the code requires the same software and package versions that were used for our analyses (all versions are indicated in the paper). As the analysis involved running a large number of methods on many data sets and over an extended period of time, we cannot guarantee that it will run successfully with new releases of the software, or that exactly the same results will be obtained with newer versions of the packages. While the repository will not be updated to ensure that it runs with every new version of the used packages, the issues can be used to post questions and/or solutions as they arise. 

**Update (13/07/2020, [Milan Malfait](https://github.com/milanmlft))**: it is now possible to run the pipeline in a Docker container (see [below](#docker-container) for details), using the (as close as possible) same package versions as originally used for the publication. This should make it easier to reproduce the results and to add more methods. Note that, currently, the Docker container is built with __R 3.3.2__ (the R version used for most of the methods from the original publication). So any additional methods added should be compatible with that version of R.


The repository contains the following information:

* `config/` contains configuration files for all the data sets that we considered. The configuration files detail the cell populations that were compared, as well as the number of cells per group used in each comparison.
* `data/` contains some of the raw data that was used for the comparison. All data sets that were used can be downloaded as a bundle from [http://imlspenticton.uzh.ch/robinson_lab/conquer\_de\_comparison/](http://imlspenticton.uzh.ch/robinson_lab/conquer_de_comparison/)
* `export_results/` contains results for the final figures, in tabular format
* `scripts/` contains all R scripts used for the evaluation
* `shiny/` contains the code for a shiny app built to browse the results ([http://imlspenticton.uzh.ch:3838/scrnaseq\_de\_evaluation](http://imlspenticton.uzh.ch:3838/scrnaseq_de_evaluation))
* `unit_tests/` contains unit tests that were used to check the calculations
* `Makefile` is the master script, which outlines the entire evaluation and calls all scripts in the appropriate order
* `include_filterings.mk`, `include_datasets.mk`, `include_methods.mk` and `plot_methods.mk` are additional makefiles listing the filter settings, data set and differential expression methods used in the comparison 
* `Dockerfile`: build instructions for the Docker container
* `install.R`: used by the Docker container to install the necessary R packages to run the pipeline
 

## Running the comparison
Assuming that all prerequisites are available, the comparison can be run by simply typing 

```$ make```

from the top directory (note, however, that this will take a **significant** amount of time!). The Makefile reads the three files *include_filterings.mk*, *include_datasets.mk* and *include_methods.mk* and performs the evaluation using the data sets, methods and filterings defined in these. The file *plot_methods.mk* detail the methods included in the final summary plots. For the code to execute properly, an *.rds* file containing a *MultiAssayExperiment* object for each data set must be provided in the `data/` directory. Such files can be downloaded, e.g., from the [`conquer`](http://imlspenticton.uzh.ch:3838/conquer/) database. The files used for the evaluation are bundled together in an archive that can be downloaded from [here](http://imlspenticton.uzh.ch/robinson_lab/conquer_de_comparison/)

## Adding a differential expression method
To add a differential expression method to the evaluation, construct a script in the form of the provided `apply_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the method. Then add the name of the method to `include_methods.mk`. To make it show up in the summary plots, add it to `plot_methods.mk` and assign it a color in `scripts/plot_setup.R`.

## Adding a data set
To add a data set, put the *.rds* file containing the *MultiArrayExperiment* object in the `data/` folder and construct a script in the form of the provided `generate_config_*.R` scripts (in the `scripts/` directory), where `*` should be the name of the data set. Then add the name of the dataset to the appropriate variables in `include_datasets.mk`. Also, add the data set to the `data/dataset_type.txt` file, indicating the type of values in each data set.

## A note on the data sets
Most data sets in the published evaluation are obtained from the [`conquer`](http://imlspenticton.uzh.ch:3838/conquer/) repository. The RPM values for the Usoskin dataset was downloaded from [http://linnarssonlab.org/drg/](http://linnarssonlab.org/drg/) on December 18, 2016. The 10X data set was downloaded from [https://support.10xgenomics.com/single-cell-gene-expression/datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets) on September 17, 2017.

## Cell cycle genes
The list of mouse cell cycle genes was obtained from [http://www.sabiosciences.com/rt_pcr_product/HTML/PAMM-020A.html](http://www.sabiosciences.com/rt_pcr_product/HTML/PAMM-020A.html) on March 9, 2017.

## Unit tests
To run all the unit tests, start `R`, load the `testthat` package and run 
``source("scripts/run_unit_tests.R")``. Alternatively, to run just the unit tests in a given file, do e.g. ``test_file("unit_tests/test_trueperformance.R", reporter = "summary")``.


## Docker container

A [Docker](https://www.docker.com/) container is provided to run the pipeline with the original versions of the methods used in the paper (see [Supplementary Table 2](https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4612/MediaObjects/41592_2018_BFnmeth4612_MOESM1_ESM.pdf) for details).

Currently the Docker container is built on R 3.3.2 and so only supports methods compatible with that version of R. This excludes *scDD*, *DEsingle* and the *zinbwave*-related methods. Furthermore, the *powsim* package used for simulations has *scDD* as a dependency and so is also not available in the Docker container. However, simulated data sets are included in the archive available [here](http://imlspenticton.uzh.ch/robinson_lab/conquer_de_comparison/).

__Note__: some of the methods in the pipeline require quite some memory (depending on the data sets being used), so make sure to assign __at least 10GB__ of memory to your Docker instance (default is 2GB). Instructions for this can be found [here](https://stackoverflow.com/a/44533437/11801854).


To run the pipeline inside the Docker container, follow these steps:

1. Clone the repository

    ```
    git clone https://github.com/csoneson/conquer_comparison.git
    cd conquer_comparison
    ```
   
2. Add data sets according to the instructions given [above](#adding-a-data-set) to the `data/` directory

3. Optional: exclude any methods and/or add new ones according to the instructions [above](#adding-a-differential-expression-method) Build image by running (you can tag it with any name you want with the `-t` option). Don't forget the `.` at the end.

    ```
    docker build -t conquer .
    ```
    
    This can take some time but needs to be performed only once. Afterwards you can just re-access the container using step 5.

5. Run the container with 

    ```
    docker run -it --rm \
    	--name conquer \
    	-v ${PWD}:/home/conquer \
    	conquer:latest
    ```
    
    The `-v` argument mounts your local *conquer* directory inside the Docker container. This is essential to make the scripts necessary to run the pipeline available inside the container and to get the results back. `-it` runs the container interactively and `--rm` automatically removes the container when you stop it.

Launching the Docker should open a bash shell inside the container, from where the pipeline can be run using `make`. Note that any changes you make in your local *conquer* directory, outside Docker, will propagate to the container because of the mounted volume. When you're done, just type `exit` in the container shell. The container will be stopped and removed but the installed image will still be present, so to re-access you only have to repeat step 5.
