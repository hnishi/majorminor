#!/bin/bash -eu

fn_out1=out_do1.txt

rm -f $fn_out1

#for i in {0..5000}
for i in {0..2}
do

echo INPUTPDB1 ../../120mM_longeq/02angle/out_do3/conf${i}.pdb > tmp.in 
echo OUTPUTFILE1 out1.txt >> tmp.in 
echo CHAIN_A_a 1  >> tmp.in 
echo CHAIN_A_b 38 >> tmp.in 
echo CHAIN_B_a 39 >> tmp.in 
echo CHAIN_B_b 76 >> tmp.in 
echo CHAIN_C_a 77 >> tmp.in 
echo CHAIN_C_b 114 >> tmp.in 
echo CHAIN_D_a 115 >> tmp.in 
echo CHAIN_D_b 152 >> tmp.in 
echo CHAIN_D_b 152 >> tmp.in 
echo RANGE_COM 0 >> tmp.in 
echo CUTOFF_DIST 25 >> tmp.in 

#./a.out tmp.in > /dev/null 2>&1 
./a.out tmp.in #> /dev/null 2>&1 
cat out1.txt >> $fn_out1 

done

rm tmp.in

cat $fn_out1
