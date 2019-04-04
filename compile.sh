#!/usr/bin/bash
module load texlive
latexmk -pdfxe -f thesis.tex
latexmk -c
