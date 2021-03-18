grep -A 1 "perpendicular cell widths:" Output/System_0/* > boxlengths_temp.txt
sed -n "2~3p" boxlengths_temp.txt > boxlengths.txt
rm boxlengths.txt
grep "volume of the cell" Output/System_0/* > volume_temp.txt
sed -n "1~2p" volume_temp.txt > volume.txt
rm volume.txt
