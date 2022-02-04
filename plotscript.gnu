set xrange[0:149]
set yrange[0:149]
set cbrange[-2:2]
set palette defined (-2 'blue', 0 'white', 2 'red')
plot 'Edata0.txt' matrix with image
do for [n=1:10000] { plot "Edata".n.".txt" matrix with image; pause 0.5 }
