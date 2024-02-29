# Statistical Models in Computational Biology

This repo contains the solution to the projects and some practice code (not necessarily in `R`) for the course Statistical Models in Computational Biology at ETH Zurich.

## File Structure

The directory names are self-explanatory. The `projects` directory contains the solutions to the course projects, while the `practice` directory contains some practice code (authored by the owner of the repo).

> For the curious ones, `create_proj_dir.py` is a utility script to create the project directories. It takes the project number as an argument and creates the project directory with the necessary files and directories (the most useful file is a proper skeleton solution.Rmd).

## How to run the course project code yourself

As requested by the lecturer, the projects were implemented in `R` and consequently it's best to use [Rstudio](https://posit.co/download/rstudio-desktop/) as the IDE to view and run the code.

Double-clicking on the `.Rproj` file will open the project in Rstudio. Or you can open Rstudio and use **File > Open Projects** to open the `.Rproj` file, if you have already opened Rstudio.

The dependencies are managed by `renv` and the projects are set up to use it. To replicate the environment, simply type:

```R
renv::restore()
```

in the R console. This will install all the necessary packages and their dependencies.

> **Note**: As of now (Feb 2024), Microsoft has terminated MRAN. Should you encounter any `curl` errors while trying to install the packages, please try `options(renv.config.mran.enabled = FALSE)` to disable MRAN.

The projects are not implemented using the same R version. It is recommended to use the following R versions for each project. Using a different version may result in errors due to package incompatibility:

- Project 01-07: R 4.2.3
- Project 08-11: R 4.3.0

## Useful Links

- [Recorded Lectures (2018: **Open to ALL**)](https://video.ethz.ch/lectures/d-bsse/2018/spring/636-0702-00L.html)
- [Quick Intro to Parallel Computing in R](https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html)
