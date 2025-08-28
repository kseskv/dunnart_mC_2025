#!/bin/bash -e

###---Make bigtable
#concatenate chr and pos columns to make a unique identifier - to be used to combine different files into a bigtable
#get methylation values from 0 to 1; get coverage (C+T)
for input in *_CpG.bedGraph; do
    echo $input
        ext='_C_T_cov.bedGraph'
        output=$(echo $input | rev | cut -c 10- | rev)$ext #n+1 - removes n last characters from the string
    echo $output
    awk '{OFS="\t"} {$7 = $5 + $6}1 {print $1,$2,$3,$4,$5,$6,$7}' $input  > $output
 done


#concatenate chr and pos columns to make a unique identifier - to be used to combine different files into a bigtable
  for input in *_C_T_cov.bedGraph; do
    echo $input
        ext='C_cov_ID.bedGraph'
        output=$(echo $input | rev | cut -c 17- | rev)$ext #n+1 - removes n last characters from the string
    echo $output
    awk '{print $1"_"$2"\t"$5"\t"$7}' < $input  > $output
 done


#remove the first line and add the header
for input in *_C_cov_ID.bedGraph; do
vim -u NONE +'1d' +wq! $input
done


#add filename to .C and .cov
for i in *_C_cov_ID.bedGraph; do 
    echo $i; s=`echo $i| awk -F'[._]' '{print $1"_"$2"_"$3}'`; 
    echo -e "ID\t${s}.C\t${s}.cov" > ${i/.bedGraph/.txt}; 
    cat ${i/.bedGraph/.txt} $i >> ${i/.bedGraph/.txt}
done

#-if there is more than one sample - sort ID column and merge based on it
# Sort ID column
# Re-sort all files for join compatibility (lexicographic order)
for file in *_CpG_C_cov_ID.txt; do
    echo "Sorting $file for join..."
    (head -1 "$file"; tail -n +2 "$file" | sort -k1,1) > "${file%.txt}.sorted.txt"
done

#test the merge
join -t$'\t' dunnart_sperm_r1_R1_paired_sorted_dedup_sorted_filter.svg_CpG_C_cov_ID.sorted.txt dunnart_sperm_r2_R1_paired_sorted_dedup_sorted_filter.svg_CpG_C_cov_ID.sorted.txt | head -5

#merge
files=(*.sorted.txt)
echo "Starting with: ${files[0]}"
cp "${files[0]}" merged.txt

for file in "${files[@]:1}"; do
    echo "Joining with $file..."
    join -t$'\t' merged.txt "$file" > temp.txt && mv temp.txt merged.txt
done

echo "Final result:"
head merged.txt
echo "Total lines: $(wc -l < merged.txt)"
echo "Total columns: $(head -1 merged.txt | tr '\t' '\n' | wc -l)"



# header
awk '{OFS="\t"}NR==1{print "chr", "start", $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25}' merged.txt > header.tsv
# data with the header !!! DON'T FORGET TO CHANGE THE NUMBER OF COLUMNS
awk '{OFS="\t"}NR>1{split($1,a,"_"); print a[1], a[2], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13,  $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25}' merged.txt >> header.tsv


mv header.tsv bigtable_unsorted.tsv
sort -k1,1 -k2,2n bigtable_unsorted.tsv -o bigtable.tsv

# check
echo "Your merged file is: $file1"
head bigtable.tsv| column -t



