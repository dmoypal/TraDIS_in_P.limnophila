echo "Total number of reads:" 
wc -l all_data.fastq | awk '{print $1/4}'

echo "Clean reads:"
wc -l all_data.only_trans.trimmer.no_adapt.fastq | awk '{print $1/4}'

echo "Unique insertion in the main chromosome:"
main=$(grep CP001744 insertions.bed | wc -l)
echo $main

echo "Unique insertion in the plasmid:"
plasmid=$(grep CP001745 insertions.bed | wc -l)
echo $plasmid

echo "Reads in CDS in main:"
cdsmain=$(cut -f2,5,10 with_trans.tnseq.tsv | grep gene | grep CP001744 | cut -f3 | awk '{sum+=$1} END {print sum}')
echo $cdsmain

echo "Read in CDS in plasmid:" 
cdsplasmid=$(cut -f2,5,10 with_trans.tnseq.tsv | grep gene | grep CP001745 |  cut -f3 | awk '{sum+=$1} END {print sum}')
echo $cdsplasmid

echo "Reads intergenic main:"
intergenicmain=$(($main - $cdsmain))
echo $intergenicmain 

echo "Reads intergenic plasmid:"
intergenicplasmid=$(($plasmid - $cdsplasmid))
echo $intergenicplasmid


mainlength=$(awk '$1=="CP001744.1"{print $2}' ./genome.sizes)
plasmidlength=$(awk '$1=="CP001745.1"{print $2}' ./genome.sizes)

echo "Density of insertions in the main in bp:"
densitymain=$(echo "scale=2; $mainlength / $main" | bc)
echo $densitymain

echo "Density of insertions in the plasmid in bp:"
densityplasmid=$(echo "scale=2; $plasmidlength / $plasmid" | bc)
echo $densityplasmid

