# mapa.gnu - Gnuplot script para mapa de convergência com 10 cores

set terminal pngcairo size 1000,800 enhanced font 'Verdana,10'
set output 'mapa.png'

set title "Mapa de Convergência - Iterações até convergência"
set xlabel "r"
set ylabel "s"
set cblabel "Número de Iterações"

# Define a paleta com 10 cores - índices de 0 a 9
set palette defined ( \
  0  '#0000FF', \
  3  '#0033FF', \
  6  '#0066FF', \
  9  '#0099FF', \
 12  '#00CCFF', \
 15  '#00FFFF', \
 18  '#33FFCC', \
 21  '#66FF99', \
 24  '#99FF66', \
 27  '#FF0000')

set cbrange [0:30]


# Ajuste os limites do gráfico conforme seu domínio
set xrange [-100:100]
set yrange [-100:100]

# Defina o intervalo da barra de cores conforme o número máximo de iterações esperado (0 a 9 aqui)
set cbrange [0:30]

set pm3d map

# Use o arquivo gerado pelo seu Fortran
splot 'mapa.dat' using 1:2:3 with pm3d notitle
