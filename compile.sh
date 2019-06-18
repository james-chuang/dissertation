#!/usr/bin/bash
module load texlive
latexmk -pdfxe -f thesis.tex
latexmk -c
latexmk -pdfxe -f presentation.tex
latexmk -c
# latexmk -pdfxe -f acknowledgements.tex
# latexmk -c
