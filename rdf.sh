#!/bin/sh

mol='CO'
filename='ether-g'

nm1=`grep "${mol}Z" ${filename}.gro | tail -1 | awk '{printf ("%d\n",$1)}'`
nm2=`grep "${mol}L" ${filename}.gro | tail -1 | awk -v n=$nm1 '{printf ("%d",$1-n)}'`
nm0=`grep "${mol}L" ${filename}.gro | tail -1 | awk '{printf ("%d",$1)}'`
na1=`grep -c " 1${mol}Z " ${filename}.gro`
na2=`grep -c " ${nm0}${mol}L " ${filename}.gro`

echo "${filename}.xtc	 ${filename}_rdf.csv 	2" > tmp.inp
echo ".1 	0 	500 " >> tmp.inp
echo "$nm1	$na1" >> tmp.inp
echo "$nm2	$na2" >> tmp.inp
echo "$nm0" >> tmp.inp
grep " 1${mol}Z " ${filename}.gro \
| sed -e "s/ H/ 1.008 /" \
| sed -e "s/ C/ 12.011 /" \
| sed -e "s/ N/ 14.007 /" \
| sed -e "s/ O/ 15.999 /" \
| sed -e "s/ F/ 18.998 /" \
| awk '{printf ("%7.3f\n",$2)}'  >> tmp.inp
grep " ${nm0}${mol}L " ${filename}.gro \
| sed -e "s/ H/ 1.008 /" \
| sed -e "s/ C/ 12.011 /" \
| sed -e "s/ N/ 14.007 /" \
| sed -e "s/ O/ 15.999 /" \
| sed -e "s/ F/ 18.998 /" \
| awk '{printf ("%7.3f\n",$2)}'  >> tmp.inp

/Users/leej355/Desktop/ClusterSize/rdf.x < tmp.inp

rm -f tmp.inp
mv rdf_${mol}L.csv ~/Desktop/rdf_analysis/csv/
