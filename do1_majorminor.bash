#!/bin/bash -eu

make 

fn_out1=out_do1.txt

rm -f $fn_out1

#for i in {16..16}
for i in {0..0}
do

#echo INPUTPDB1 ../../120mM_longeq/02angle/out_do3/conf${i}.pdb > tmp.in 
echo INPUTPDB1 ./conf394.pdb > tmp.in 
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
echo DIST_BASE 8 >> tmp.in #cutoff distance to define major-backbone mode -> distance between DG/O6 and P atoms


#./a.out tmp.in > /dev/null 2>&1 
./a.out tmp.in #> /dev/null 2>&1 
cat out1.txt >> $fn_out1 

done

rm tmp.in

echo mode backbone-backbone_residue_pair1 distance_of_pair1 backbone-backbone_residue_pair2 distance_of_pair2 base-backbone_pair3 distance_of_pair3 input_file
cat $fn_out1
