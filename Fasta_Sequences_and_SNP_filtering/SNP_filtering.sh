##For start and end point comparisons for each family. The code below is an example for one of these comparisons (cross TX1)
#exact allele frequency files were created using the popoolation2 script snp-frequency-diff.pl
grep "pop" AR2-SE_exact_rc > AR2-SE_exact_rc_pop #we want only allele frequencies identified between the start and end points and not the reference chromosome


cut -f1,2,3,4,6,7,8,9,10,11,12,13 AR2-SE_exact_rc_pop > AR2-SE_exact_rc_pop2 #getting rid of columns irrelevant for analysis
sed 's/\//\t/g' AR2-SE_exact_rc_pop2 > AR2-SE_exact_rc_pop3 #splitting numerator and denominator of frequencies into their own columns 

awk -v OFS='\t' '{$17 = $13 / $14}1' AR2-SE_exact_rc_pop3 > AR2-SE_exact_rc_pop4 #producing decimals for allele frequencies by dividing numerator by denominator (start frequency)
awk -v OFS='\t' '{$18 = $15 / $16}1' AR2-SE_exact_rc_pop4 > AR2-SE_exact_rc_pop5 #(end frequency)

awk '{ total += $17 } END { print total/NR }' AR2-SE_exact_rc_pop5 #check average allele frequency for starting samples, should be around 0.25
awk -v OFS='\t' '$17<0.35' AR2-SE_exact_rc_pop5 > AR2-SE_exact_rc_pop6 #filtering starting frequency to be between 0.15 and 0.35
awk -v OFS='\t' '$17>0.15' AR2-SE_exact_rc_pop6 > AR2-SE_exact_rc_pop7

awk -v OFS='\t' '{$19 = $17 - $18}1' AR2-SE_exact_rc_pop7 > AR2-SE_exact_rc_pop8 #calculating difference in frequency between start and end and putting difference into a new column

awk -v OFS='\t' '$19<0' AR2-SE_exact_rc_pop8 > AR2-SE_exact_rc_pop9 #selecting SNPs that only increase in frequency after selection

##fet files were creasted using the popoolation2 Fisher's exact test script fisher-test.pl
cut -c-9 AR2-SE.fet > AR2-SE_edit.fet #the following lines contain code for editing first column to separate gene name and start and end position into their own columns for organization
cut -c10- AR2-SE.fet > AR2-SE_edit2.fet 
paste -d '\t' AR2-SE_edit.fet AR2-SE_edit2.fet > AR2-SE_edit3.fet

cut -c-30 AR2-SE_edit3.fet > AR2-SE_edit4.fet
cut -c31- AR2-SE_edit3.fet > AR2-SE_edit5.fet
paste -d '\t' AR2-SE_edit4.fet AR2-SE_edit5.fet > AR2-SE_edit6.fet

cut –f1,2,3 AR2-SE_edit6.fet > AR2-SE_longname.fet
cut -f4,5,6,7,8 AR2-SE_edit6.fet > AR2-SE_longname2.fet
paste -d '-' AR2-SE_longname.fet AR2-SE_longname2.fet > AR2-SE_longname3.fet


cut -c-28 AR2-SE_exact_rc_pop9 > AR2-SE_exact_rc_pop9_edit #doing the same as above, editing first column to separate gene name and start and end position into their own columns for all the files that contain the exact allele frequency differences
cut -c29- AR2-SE_exact_rc_pop9 > AR2-SE_exact_rc_pop9_edit2
paste -d '\t' AR2-SE_exact_rc_pop9_edit AR2-SE_exact_rc_pop9_edit2 > AR2-SE_exact_rc_pop9_edit3

cut -f1,2 AR2-SE_exact_rc_pop9_edit3 > AR2-SE_exact_rc_pop9_edit4
cut -f3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 AR2-SE_exact_rc_pop9_edit3 > AR2-SE_exact_rc_pop9_edit5
paste -d '-' AR2-SE_exact_rc_pop9_edit4 AR2-SE_exact_rc_pop9_edit5 > AR2-SE_exact_rc_pop9_edit6

#the fet file contains the pvalues for allele frequency differences and the exact_rc file contains info on which of these SNPs increased in frequency. We filtered the fet file for SNPs that only increased in frequency after selection.
sig_list=($(cat ./AR2-SE_exact_rc_pop9_edit6 | cut -f2 ))
for i in "${sig_list[@]}"
do
        grep –w $i ../AR2-SE_longname3.fet >> AR2-SE_sig.fet
done

cut -c-45 AR2-SE_sig.fet > AR2-SE_sig_edit.fet #final formatting
cut -c47- AR2-SE_sig.fet > AR2-SE_sig_edit2.fet
paste -d '\t' AR2-SE_sig_edit AR2-SE_sig_edit2 > AR2-SE_sig_edit3


#following works on command line, have not tested over qsub
sed 's/-/\t/g' AR2-SE_sig_edit3.fet > AR2-SE_sig_edit4.fet
sed 's/-/\t/g' AR3-SE_sig_edit3.fet > AR3-SE_sig_edit4.fet
sed 's/-/\t/g' AR4-SE_sig_edit3.fet > AR4-SE_sig_edit4.fet
sed 's/-/\t/g' AR5-SE_sig_edit3.fet > AR5-SE_sig_edit4.fet
sed 's/-/\t/g' SL1-SE_sig_edit3.fet > SL1-SE_sig_edit4.fet
sed 's/-/\t/g' SL2-SE_sig_edit3.fet > SL2-SE_sig_edit4.fet
sed 's/-/\t/g' SL3-SE_sig_edit3.fet > SL3-SE_sig_edit4.fet
sed 's/-/\t/g' SLNS3-SE_sig_edit3.fet > SLNS3-SE_sig_edit4.fet
sed 's/-/\t/g' VB2-SE_sig_edit3.fet > VB2-SE_sig_edit4.fet
sed 's/-/\t/g' VB3-SE_sig_edit3.fet > VB3-SE_sig_edit4.fet
sed 's/-/\t/g' VB3_R2-SE_sig_edit3.fet > VB3_R2-SE_sig_edit4.fet

sed 's/'1:2='//g' VB3_R2-SE_sig_edit4.fet > VB3_R2-SE_sig_edit5.fet
sed 's/'1:2='//g' VB3-SE_sig_edit4.fet > VB3-SE_sig_edit5.fet
sed 's/'1:2='//g' VB2-SE_sig_edit4.fet > VB2-SE_sig_edit5.fet
sed 's/'1:2='//g' SL1-SE_sig_edit4.fet > SL1-SE_sig_edit5.fet
sed 's/'1:2='//g' SL2-SE_sig_edit4.fet > SL2-SE_sig_edit5.fet
sed 's/'1:2='//g' SL3-SE_sig_edit4.fet > SL3-SE_sig_edit5.fet
sed 's/'1:2='//g' SLNS3-SE_sig_edit4.fet > SLNS3-SE_sig_edit5.fet
sed 's/'1:2='//g' AR2-SE_sig_edit4.fet > AR2-SE_sig_edit5.fet
sed 's/'1:2='//g' AR3-SE_sig_edit4.fet > AR3-SE_sig_edit5.fet
sed 's/'1:2='//g' AR4-SE_sig_edit4.fet > AR4-SE_sig_edit5.fet
sed 's/'1:2='//g' AR5-SE_sig_edit4.fet > AR5-SE_sig_edit5.fet

date
exit
