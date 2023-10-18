# Essential gene complement of *Planctopirus limnophila* from the bacterial phylum *Planctomycetes* 


In this repository you can find the code used for processing the raw sequence data, and for predicting essential genes using TraDIS in *P. limnophila*, based on Goodall, E.C.A. et al. "The Essential Genome of Escherichia coli K-12." mBio 9.1 (2018): e02096-17.

## How to use the code
Make sure all the necessary software and packages are installed. 
Then, download the [raw data](https://doi.org/10.6084/m9.figshare.24249346) and run in the following order: 
1. ```insertions_pipeline.sh```
2. data_processing.R
3. stats_insertions.sh
4. tradis.R
5. stats_essentiality.sh
6. domain_essential.py

## Credits

This project has been developed thanks to the effort and contribution of the following individuals:

- **David Moyano Palazuelo** - [Link to GitHub Profile](https://github.com/dmoypal)
- **Ildefonso Cases**

We appreciate everyone who has dedicated their time and expertise to improve this project and make it possible.
