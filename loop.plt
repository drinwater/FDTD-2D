set xrange [0:149]
set yrange [0:149]
set cbrange [-0.3:0.3]
set palette defined (-0.3 "blue", 0 "white", 0.3 "red")

plot sprintf(".\\data\\Edata%d.txt",i) matrix with image

i=i+1
print i
set label 1 sprintf("%d",i) at 130,200
if(i<=n) reread

#end of loop.plt