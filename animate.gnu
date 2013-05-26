set parametric
set hidden3d
set nokey
set term gif animate delay 10
set output "animate.gif"
set xlabel "East/x-axis (m)"
set ylabel "South/y-axis (m)"
set zlabel "Altitude/z-axis (m)"
xrot=60
zrot=0
splot "coordinate.dat" with lines
load "rot.gnu"
set output
