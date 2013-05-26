i=1
xrot=(xrot+5)%180
zrot=(zrot+10)%360
set view xrot,zrot
i=i+1
replot
if (i<12) reread

