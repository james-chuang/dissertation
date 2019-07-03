# My PhD dissertation and defense slides

This is the code I used to generate my PhD dissertation (`thesis.pdf`) and the slides I used at my defense (`presentation.pdf`).
The raw data is not included since some of it is not yet published, so this repo won't run by itself and is provided mostly for reference.
In the future, I may upload a complete archive including raw data to [Zenodo](https://zenodo.org/).

How it works:

    0. Figures in PDF format are generated by R scripts called through a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow:
        - Paths to input data files are specified in `config.yaml`
        - Snakemake rules for generating individual figures are specified in `rules/`
        - Conda virtual environments that Snakemake uses to run rules are specified in `envs/`
    1. After all of the figures have been generated, Snakemake runs `compile.sh`, which uses TeX Live to typeset the thesis and presentation slides
        - The TeX files are `thesis.tex`, `presentation.tex`, and the TeX files `tex/`
        - Citation is done using the natbib package, with `.bib` files in `references/`
        - Layout of the dissertation is handled by options set in `thesis.tex` as well as `bu_ece_thesis_minimal.sty`, a stripped down version of `bu_ece_thesis.sty` from the [BU Engineering LaTeX template](http://www.bu.edu/eng/departments/ece/resourcesforcurrentstudent/ece-ms-and-phd-thesis-prep/ms-thesis-phd-dissertation/).

