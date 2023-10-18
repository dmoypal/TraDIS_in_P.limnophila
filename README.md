# Essential gene complement of *Planctopirus limnophila* from the bacterial phylum *Planctomycetes* 


In this repository you can find the code used for processing the raw sequence data, and for predicting essential genes using TraDIS in *P. limnophila*, based on Goodall, Emily C A et al. “The Essential Genome of Escherichia coli K-12.” mBio vol. 9,1 e02096-17. 20 Feb. 2018, doi:10.1128/mBio.02096-17. 

#How to use the code
Make sure all the necessary programs and packages are installed. 
Then, download the raw data [(here)]([https://doi.org/10.6084/m9.figshare.24249346]) and run in the following order: 
1. insertions_pipeline.sh
2. data_processing.R
3. stats_insertions.sh
4. tradis.R
5. stats_essentiality.sh
6. domain_essential.py

