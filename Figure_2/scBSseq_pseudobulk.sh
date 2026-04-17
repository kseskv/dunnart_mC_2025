#!/bin/bash

# get list of chromosomes in the input files
CHROMS=$(awk '{print $1}' *.bismark.cov | sort | uniq)

# directory for per-chromosome outputs 
mkdir -p chrom_cov_files

# process each chromosome individually
for chr in $CHROMS; do
    echo "Processing chromosome: $chr"

    awk -v chr="$chr" '
    FNR==1 {file++}
    $1 == chr {
        key=$1 FS $2 FS $3;
        cov5[key,file]=$5;
        cov6[key,file]=$6;
        keys[key]++;
    }
    END {
        for (k in keys) {
            sum5=0; sum6=0;
            for (i=1; i<=file; i++) {
                sum5 += cov5[k,i];
                sum6 += cov6[k,i];
            }
            ratio = (sum5 + sum6) > 0 ? sum5 / (sum5 + sum6) : 0;
            split(k, a, FS);
            print a[1] "\t" a[2] "\t" a[3] "\t" ratio "\t" sum5 "\t" sum6;
        }
    }' *.bismark.cov | sort -k2,2n > "chrom_cov_files/scBSseq_E9_pseudobulk.${chr}.bismark.cov"
done

# merge and sort all chromosome files
echo "Merging and sorting final file..."
cat chrom_cov_files/scBSseq_E9_pseudobulk.*.bismark.cov | sort -k1,1 -k2,2n > scBSseq_E9_pseudobulk.bismark.cov

echo "Done."
