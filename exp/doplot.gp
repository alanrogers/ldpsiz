set terminal pdf
set out "doplot.pdf"
set autoscale
#set logscale y
set xlabel "Recombinational Distance (cM) Between SNPs"
set ylabel "Sigma_d^2"
#set yrange [0:0.12]
#set style line 1 lt 1 lw 2 lc rgb "red"
#set style line 2 lt 1 lw 1 lc rgb "red"
#set style line 3 lt 1 lw 2 lc rgb "green"
set title "Decay of LD with Distance along Chromosome"
#      "macs.preld" using 1:3 title "Equil Epoch 0" with lines,
#      "macs.preld" using 1:4 title "Equil Epoch 1" with lines,
plot  "macs.eld" using 1:2 title "macs" with linespoints, \
      "macs.sald" using 1:2 title "fit" with linespoints, \
      "preld.out" using 1:2 title "Strobeck & Morgan" with linespoints, \
      "preld.out" using 1:3 title "equil 0" with linespoints, \
      "preld.out" using 1:4 title "equil 1" with linespoints
exit
