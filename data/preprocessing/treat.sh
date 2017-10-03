#!/bin/bash

for f in /storage/mathelierarea/processed/arnaudst/TFFM_DAPSeq/data/data_to_treat/*;
do 
    [ -d $f ] && cd "$f" && pwd
    TF=`pwd | awk 'BEGIN{FS="/"}{print $9}'`
    # awk 'BEGIN{OFS="\t"}{print $1,$1,$3} ' chr1-5_GEM_events.narrowPeak > ${TF}.bed 
    mv peaks.bed ${TF}.bed
    mv meme_m1.txt ${TF}.meme
done; 
