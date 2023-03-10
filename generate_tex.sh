#!/bin/bash

> ./results.tex

echo "\documentclass{article}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{float}
\usepackage[export]{adjustbox}

\title{Results}
\date{}

\begin{document}
\maketitle

" >> ./results.tex

for entry in ./output/*
do

    meshname=$(echo "$entry" | grep -o -P '(?<=\.\/output/).*')
    meshname=$(echo "$meshname" | sed 's/_/\\_/g')
    echo "$meshname"
    #TODO better read
    line1=$(sed -n '1,1p;1q' < ${entry}/results-uv-at.dat)
    nv=$(echo "$line1" | sed 's/\(.*\),.*/\1/')
    nf=$(echo "$line1" | sed 's/.*, \(.*\)/\1/')
    line2=$(sed -n '2,2p;2q' < ${entry}/results-uv-at.dat)
    e_dist=$(echo "$line2" | sed 's/\(.*\),.*,.*/\1/')
    e_cut=$(echo "$line2" | sed 's/.*, \(.*\), .*/\1/')
    duration=$(echo "$line2" | sed 's/.*, .*, \(.*\)/\1/')
    echo "
\begin{figure}[H]
\centering
        \includegraphics[max height=0.1\textheight,max width=0.24\linewidth,keepaspectratio]{${entry}/param.jpg}
        \includegraphics[max height=0.1\textheight,max width=0.24\linewidth,keepaspectratio]{${entry}/mesh.jpg}
        \caption*{\textbf{$meshname} ($nv vertices, $nf faces) \\\\
        \$ E_d = $e_dist \$ \$ E_s = $e_cut \$ \$ d = $duration ms \$
    }
\end{figure}
    " >> ./results.tex
done

echo "\end{document}" >> ./results.tex
