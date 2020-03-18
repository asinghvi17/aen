#!/bin/zsh

for x in $(ls | grep CSV)
do
		cat $x | cut -d',' -f4-5 > ${str: -7 : -4}.csv
done
