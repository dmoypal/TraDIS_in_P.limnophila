#!/bin/bash

nonessential=$(grep "non-essential" essentiality.tsv | wc -l)
essential=$(grep essential essentiality.tsv | grep -v non | wc -l)
unclear=$(grep unclear essentiality.tsv | wc -l)
total=$((essential + unclear + nonessential))

echo "Total genes and pseudogenes: $total"
echo "Essential: $essential"
echo "Non-essential: $nonessential"
echo "Unclear: $unclear"

essential_percentage=$(printf "%.3f" $(echo "scale=3; $essential * 100 / $total" | bc))
unclear_percentage=$(printf "%.3f" $(echo "scale=3; $unclear * 100 / $total" | bc))
nonessential_percentage=$(printf "%.3f" $(echo "scale=3; $nonessential * 100 / $total" | bc))

echo "Essential pecentage: $essential_percentage%"
echo "Unclear percentage: $unclear_percentage%"
echo "Non-essential percentage: $nonessential_percentage%"

grep essential essentiality.tsv | grep -v non | cut -f1 > list_essential
grep non-essential essentiality.tsv | 
cut -f1 > list_non-essential
grep unclear essentiality.tsv | cut -f1 > list_unclear
