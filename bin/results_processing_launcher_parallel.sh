# the type of algorithm to be used in the peak processing


#usage for TFFM: bash results_processing_launcher_parallel.sh /storage/scratch/marius/ENCODE/bin/peak_processing/peak_processing_TFFM_one.sh /storage/scratch/marius/ENCODE/data/PFM_JASPAR batch_1-41


#usage for NRG: bash results_processing_launcher_parallel.sh /storage/scratch/marius/ENCODE/bin/peak_processing/peak_processing_energy_one.sh /storage/scratch/marius/ENCODE/data/MEME_JASPAR_2016 batch_1-41

script=$1
pcm_folder=$2
batch=""

#main loop over input folders
peaks_folder="/storage/scratch/marius/GEO/results/reprocess_tffm_no_profile"

# cd $peaks_folder/MACS/$batch
cd $peaks_folder
i=0

for d in */
   do
     
     ((i++))

     #limit the number of processes
 #   if [ "$i" -lt "250" ]
 #       then
 #            continue
 #    elif [ "$i" -gt "300" ]
 #        then
 #           break
 #    fi

#    if [ "$i" -gt "40" ]
#        then
#          break
#     fi


    #  time bash $script MACS $peaks_folder/MACS/$batch $pcm_folder /storage/scratch/marius/GEO/bin/GEO_TF_JASPAR_mapping_10_09_2017.tsv 250 $d & 
      time bash $script MACS $peaks_folder $pcm_folder /storage/scratch/marius/GEO/bin/tffm_new_profiles/GEO_TF_JASPAR_mapping_19_09_2017.tsv 250 $d & 
#      echo "bash $script MACS $peaks_folder/MACS/$batch $pcm_folder /storage/scratch/marius/ENCODE/bin/ENCODE_TF_JASPAR_mapping.tsv 500 $d &"
done
