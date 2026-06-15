#!/bin/bash
for i in {1..100}
do
        cat ./$i/isoc/isoc.bestlhoods | awk '{print $7}' | tail -n 1 >>t
        echo $i >>t2
done
paste t t2 >sum.txt
sort -k1,1 -nr sum.txt
