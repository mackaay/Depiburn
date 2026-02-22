#!/bin/bash

kmercountexact.sh -Xmx16g in=hg38.fa k=4 out=hg38_4mer_counts.txt


INPUT="hg38_4mer_counts.txt"
OUTPUT="human_genome_4mer_expected.csv"

# 1. Reformat: Put motif and count on the same line, remove the '>'
# This turns: 
# >21643152
# GCTT
# Into:
# GCTT 21643152
awk '/^>/{count=substr($0,2); next} {print $0, count}' $INPUT > temp_counts.txt

# 2. Calculate the Grand Total of all k-mers
TOTAL=$(awk '{sum+=$2} END {print sum}' temp_counts.txt)

# 3. Create the CSV with Header
echo "Motif,Count,Expected_Freq" > $OUTPUT

# 4. Calculate Relative Frequency and append to CSV
awk -v total="$TOTAL" 'BEGIN{OFS=","} {print $1, $2, $2/total}' temp_counts.txt >> $OUTPUT

# Cleanup
rm temp_counts.txt

echo "Done! Created $OUTPUT"
echo "Total k-mers processed: $TOTAL"
