D=../OrionMuse/LineMaps
linelist=$D/mean-*[0-9][0-9][0-9][0-9]-patfixx.fits
for line in $linelist; do
    echo "Processing $line"
    time python multibin-map.py $line > ${line}-multibin.log
done
