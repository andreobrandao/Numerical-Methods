set title "Erro relativo - Secante vs Müller"
set xlabel "Iteração"
set ylabel "Erro relativo (%)"
set grid
set key left top
set format y "%.1e"
plot "saida.dat" using 1:2 with linespoints title "Secante", \
     "saida.dat" using 1:3 with linespoints title "Müller"
pause -1
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'grafico.png'
