## piTracer - Automatic reconstruction of molecular cascades for the identification of synergistic drug targets

------------------------------------------------------------------------

The piTracer framework reconstructs molecular cascades between two genes, two metabolites, or a gene and a metabolite based on an underlying multi-omics pathway network. Its drug ranking functionality can identify synthetic lethal candidates based on a differential metabolomics experiment after single agent treatment.

## Citation

Gomari D, Achkar I, Benedetti E, Tabling J, Halama A, Krumsiek J. "piTracer - Automatic reconstruction of molecular cascades for the identification of synergistic drug targets". *biorXiv* [link]

## Hosted app

A hosted version of the Shiny app version of piTracer is can be found [here](https://cbsunemo.biohpc.cornell.edu/).

## Code repository

Clone this repository to your local machine

`git clone https://github.com/krumsieklab/piTracer.git`

### Package setup

piTracer was tested on R 4.2.1 with the corresponding latest packages as of early 2023. We recommend using the package versions stored using [renv](https://rstudio.github.io/renv/articles/renv.html). To set this up, open the piTracer project file `piTracer.proj` in RStudio, and run:

`renv::restore()`

### Precalculated data

The app requires the precalculated binary files from the `precalculated_data` folder to run. To clone those files as part of the repository, you need to have [git-lfs](https://git-lfs.com/) set up. If you don't have access to git-lfs, you will need to download the files manually.


### Running the Shiny app locally

Run `run_piTracer.R` or excute to launch the Shiny app in your local browser.

### Code version of gene ranking algorithm

The script `/src/gene_ranking/execute_generanking.R` runs the ranking algorithms programmatically based on an input spreadsheet.

### Source code overview

| Folder              | Description                                                                                           |
|:--------------------|:------------------------------------------------------------------------------------------------------|
| `src/global.R`      | used by various scripts to set up the tracing engine, Shiny app modules, and precalculated data files |
| `src/gene_ranking/` | contains the drug combination prediction functionality                                                |
| `src/shiny_app/`    | contains UI and server modules of the Shiny app                                                       |
| `src/tracing/`      | contains scripts that run the core tracing steps                                                      |
