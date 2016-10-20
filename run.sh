#!/bin/bash
for n in 10 20 30 40 50 60 70
do
	for alfa in 0.2 0.4 0.6 0.8
	do
		./bin/main ./dataset/autralian_post/ap${n}_2.txt ${alfa} 0.75 100
	done
done
#./bin/main ./dataset/autralian_post/ap10_2.txt 0.2 0.75 100
