#!/bin/bash

if [[ $1 -eq "d" ]]; then
  grep Procesor tests/testx_methane.out | grep iar | sed s/\,/\ \/g > tfil.txt
elif [[ $2 -eq "d" ]]; then
  grep Procesor\ \ $1 tests/testx_methane.out | grep iar | sed s/\,/\ \/g > tfil.txt
else
  grep Procesor\ \ $1 tests/testx_methane.out | sed s/\,/\ \/g | awk "{if (\$8==$2){print \$0}}" > tfil.txt
fi

awk '{print $9, $14}' tfil.txt > pfil.txt

if [[ $3 ]]; then
    echo "set term postscript enhanced color; set output 'pfil.eps'; plot 'pfil.txt' u 1:2, $3/x;" | gnuplot
else
    echo "set term postscript enhanced color; set output 'pfil.eps'; plot 'pfil.txt' u 1:2;" | gnuplot
fi
