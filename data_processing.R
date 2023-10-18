# Author: Ildefonso Cases

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")  
#BiocManager::install("GenomicAlignments", force = TRUE)
#BiocManager::install("GenomicFeatures", force = TRUE)
#install.packages("tidyverse")
#install.packages("future")
#install.packages("furrr")

library(tidyverse)
library(future)
library(furrr)
library(GenomicAlignments)
library(GenomicFeatures)



## calculate insertion points and create bed
dir="/Users/davidmoyanopalazuelo/Dropbox/proyectos/prueba_que_los_scripts_de_github_funcionan_bien/"
bam_file="trimmed2_t2_b2_all.sorted.dedup.bam"
aln_nodup1 <- readGAlignments(paste0(dir,bam_file))
insertions1 <- as.data.frame(aln_nodup1) %>% 
  mutate(insertion_point=if_else(strand=='+',start-1L,end-1L),
         insertion_end=insertion_point+1L,
         frame=(insertion_point %% 3 +1) * if_else(strand=='+',1,-1),
         name='read',
         score=0)

insertions <- insertions1 %>%
  dplyr::select(seqnames,insertion_point,insertion_end,name,score) %>%
  unique()

write_tsv(insertions %>% 
            dplyr::select(seqnames,start=insertion_point,end=insertion_end,name,score) %>% 
            mutate(strand="+"),
          paste0(dir,"insertions.bed"),col_names = F)

count <- function(gene){
  overlap <- insertions %>%
    filter(seqnames == gene$chrom,
           insertion_point >= gene$start,
           insertion_point <= gene$end)
  
  nrow(overlap)
}

annot <- read_tsv(paste0(dir,"Plimno.gff3"),
                  comment="#",
                  col_names=c('chrom',"Source","type","start","end","frame","strand","score","info")) %>% 
  mutate(frame=(start %% 3 +1) * if_else(strand=='+',1,-1))

#plan(multiprocess, workers=15)
future::plan(future::multisession, workers = 10)

annot_t <- nest(annot, data=c(chrom,start,end,strand,frame)) 
annot_t <- annot_t %>%
  mutate(count_all=future_map_dbl(data,count))

annot_t <- annot_t %>% unnest(cols = c(data))

insertions_main <- insertions %>% filter(seqnames=='CP001744.1')
insertions_plasmid <- insertions %>% filter(seqnames=='CP001745.1')

annot_t <- annot_t %>% 
  mutate(count_all_per_kb=log((count_all+1)/(end-start)*1000)) 

write_tsv(annot_t,paste0(dir,"with_trans.tnseq.tsv"))
