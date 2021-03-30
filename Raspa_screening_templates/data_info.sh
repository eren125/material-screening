#!/bin/bash
grep -A 1 "perpendicular cell widths:" Output/System_0/* > boxlengths_temp.txt
sed -n "2~3p" boxlengths_temp.txt > boxlengths.txt
rm boxlengths.txt
grep "volume of the cell" Output/System_0/* > volume_temp.txt
sed -n "1~2p" volume_temp.txt > volume.txt
rm volume.txt
grep -A 1 "Framework Mass:" Output/System_0/* > framework_temp.txt
sed -n "1~6p;2~6p" framework_temp.txt > framework.txt
rm framework_temp.txt
grep -A 3 "Box\[0\]" Output/System_0/* > matrix_temp.txt
sed -n "2~5p;3~5p;4~5p" matrix_temp.txt > matrix.txt
rm matrix_temp.txt

python3 merge_info.py 
