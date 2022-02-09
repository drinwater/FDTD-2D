#start of out.plt file
#This is used to load the loop file
#set up basis config

i=1
n=1000
set title "Near-field movie"
set size ratio 1
set terminal gif animate delay 15 loop 10000
set output "Edata.gif"
load "loop.plt"
set output

#end of out.plt file