#! /usr/bin/bash

gnuplot 'plot.gnuplot'
printf "\n--- done gnuplot ---\n\n"

for i in {1..12}; do
  printf "\\\\documentclass[a4paper]{article}\n\\\\usepackage{graphicx}\n\\\\usepackage{xcolor}\n\\\\begin{document}\n\\\\pagenumbering{gobble}\n\\\\centering\\\\include{output_%d.tex}\n\\\\end{document}" "$i" > plot.tex;
  latex plot.tex;
  printf "\n--- done latex ---\n\n";
  pdflatex plot.tex;
  printf "\n--- done ---\n\n";
  pdfcrop --margins 10 plot.pdf plot.pdf;
  mv plot.pdf plot_$i.pdf;
done

latex correlation.tex;
pdflatex correlation.tex;
pdfcrop --margins 10 correlation.pdf correlation.pdf
