# TFFM_plants
This contains all the scripts and the data to create TFFm from the DAP Seq  Data generated by OMalley.
All the TF studied here are edible to jaspar.

######

cd bin # go to bin

./peak_processing_TFFM_arnaud.sh # launch the creation of all the TFFms

\# needs mapping file in ../data/map.csv

\# needs peaks in ../data/narrowPeaks

\# needs pfms in ../data/jaspar_pfm

\# use TFFM_first_order.py which launches the TFFM program in ../lib/TFFM

########

Read the comments in peak_processing_TFFM for more infos 

output files : ../results/processed
