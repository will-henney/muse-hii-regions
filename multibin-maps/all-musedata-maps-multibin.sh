D=~/tmp/musedata
for line in $D/{sigma,mean,linesum}-*-[0-9][0-9][0-9][0-9].fits; do
    echo "Processing $line"
    time python multibin-map.py $line > ${line}-multibin.log 2>&1
done
