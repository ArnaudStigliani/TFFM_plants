#!/bin/bash
#usage: bash peak_processing_TFFM.sh JAMM /storage/scratch/marius/ENCODE/results/JAMM/batch_1-41 /storage/scratch/marius/ENCODE/data/PFM_JASPAR /storage/scratch/marius/ENCODE/bin/ENCODE_TF_JASPAR_mapping.tsv 500
workspace="/storage/scratch/marius/ENCODE/results"

#parse arguments
data_type=;
if [ -z "$1" ]
    then
        echo "No data type specified. Exiting...";
        exit 1
    else
        data_type=$1
fi

in_dir=
if [ -z "$2" ]
    then
        echo "No input file folder specified. Exiting...";
        exit 1
    else
        in_dir=$2
fi

pfm_folder=;
if [ -z "$3" ]
    then
        echo "No PFM folder specified. Using: /storage/scratch/marius/ENCODE/data/PFM_JASPAR";
        pfm_folder="/storage/scratch/marius/ENCODE/data/PFM_JASPAR"
    else
        pfm_folder=$3
fi

tf_pfm_map=;
if [ -z "$4" ]
    then
        echo "No TF to PFM mapping file specified. Using: /storage/scratch/marius/ENCODE/bin/ENCODE_TF_JASPAR_mapping_11_08_2017.tsv";
        tf_pfm_map="/storage/scratch/marius/ENCODE/bin/ENCODE_TF_JASPAR_mapping_11_08_2017.tsv"
    else
        tf_pfm_map=$4
fi


window_size=;
if [ -z "$5" ]
    then
        echo "No window size specified. Using 500bp";
        window_size=500
    else
        window_size=$5
fi

########## PARAMETER SETTING  ##################
result_folder="/storage/scratch/marius/ENCODE/results/processed_peaks_new_profiles/TFFM/comparison"
# faulty data sets
no_peaks_folder="/storage/scratch/marius/ENCODE/results/processed_peaks_new_profiles/TFFM/no_peaks"
no_mapping_folder="/storage/scratch/marius/ENCODE/results/processed_peaks_new_profiles/TFFM/no_mapping"
no_thresholds_folder="/storage/scratch/marius/ENCODE/results/processed_peaks_new_profiles/TFFM/no_thresholds"
# script locations
tffm_module="/storage/scratch/marius/ENCODE/bin/TFFM_first_order.py"
compute_thresholds="/storage/scratch/marius/ENCODE/bin/compute_thresholds_entropy.R"
get_thresholds="/storage/scratch/marius/ENCODE/bin/get_thresholds.R"

hg_chr_file="/storage/scratch/marius/ENCODE/data/hg38_chromosome_sizes"
hg_gen_file="/storage/mathelierarea/raw/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"

TRAINING_WINDOW_SIZE=50 #the size of the window on which the TFFM should train

############ END OF PARAMETER SETTING #########

#print script arguments
echo "------------------------------------------------------"
echo "Data type: $data_type"
echo "Window size: $window_size"
echo "Input folder: $in_dir"
echo "PFM folder: $pfm_folder"
echo "Workspace: $workspace"
echo "------------------------------------------------------"

mkdir -p $result_folder/$(basename $in_dir)

#counters
no_peaks_counter=0
no_thresholds_counter=0
no_pfm_counter=0
no_mapping_counter=0

#log the data sets with no TF to PFM mapping entry
log_no_mapping=()
log_no_mapping_filename=$no_mapping_folder/data_sets_no_mapping.log;
# read in all entries if log file exists
if [ -e $log_no_mapping_filename ]
  then
     while IFS= read -r line; do log_no_mapping+=($line); done < $log_no_mapping_filename
  else
     mkdir -p $no_mapping_folder
fi

#function to check if an element is NOT in an array
elementNotPresent () {
  local i
  for i in "${@:2}"; do [[ "$i" == "$1" ]] && return 1; done
  return 0
}

#sort the TF to PFM mapping file by data set (change column if not the 11th)
if [ ! -e ${tf_pfm_map}.sorted ]
   then
      sort -t $'\t' -k 11,11 $tf_pfm_map > ${tf_pfm_map}.sorted
fi
tf_pfm_map=${tf_pfm_map}.sorted

# get ready to process
cd $in_dir
i=1
nb_folders="$(ls -l $in_dir | grep -c ^d)"

#main loop over input folders
for d in */
   do
   
   echo "Processing folder $i out of $nb_folders"
   ((i++))

   #create folder structure and output file name
   out_dir=$result_folder/$(basename $in_dir)/${d%/*}/$data_type
   
   #check if the folder already exists and skip it
   if [ -d "$out_dir" ]
        then
	echo "Folder $out_dir already exists in results. Skipping...";
        continue
   fi

   ##### Parse the data set folder name ########
   # Elements are separated by '_' (e.g., ENCSR000ATW_K562_CBX8)
   # 1. Data set name (e.g., ENCSR000ATW)
   # 2. Cell line (e.g., K562)
   # 3. TF name (e.g., CBX8) 
   IN=$d
   folderIN=(${IN//_/ })
   data_set=${folderIN[0]%/*}
   cell_line=${folderIN[1]%/*}
   tf_name=${folderIN[-1]%/*}

   #in case the TF name is composed with some specificity of the data set
   tf_nameIN=(${tf_name//-/ })
   if [ "${#tf_nameIN[@]}" -gt "1" ]
      then
         tf_name=${tf_nameIN[-1]}
   fi

   #to upper case 
   data_set="$(echo "$data_set" | tr '[:lower:]' '[:upper:]')"
   cell_line="$(echo "$cell_line" | tr '[:lower:]' '[:upper:]')"
   tf_name="$(echo "$tf_name" | tr '[:lower:]' '[:upper:]')"

   mkdir -p $out_dir
   out_file=$out_dir/${data_type}_${tf_name}_peaks

   #check if JAMM had input in this folder
   if [ $data_type = "JAMM" ]
      then
      #echo "I'm in JAMM";
          if [ -d "$d/peaks" ]
            then
               peaks_file=$d/peaks/filtered.peaks.narrowPeak
               if [ -e $peaks_file ]
		 then
            	     ############## 1. Add offset to chrStart #################
                     awk -v OFS='\t' '{print $1,$2+$10,$2+$10+1,$4,$5,$6}' $peaks_file > $out_file
                     # hg_chr_file=$hg38_chr
                     # hg_gen_file=$hg38_gen
                 else
                     echo "Folder $d does not contain a filtered peaks file. Check the JAMM log file for more details. Skipping..";
                     move_folder=$no_peaks_folder/$data_type
                     mkdir -p $move_folder
                     cp -r $d $move_folder
                     ((no_peaks_counter++))  
                     rm -R $out_dir
                     #remove the parent directory also if empty 
                     if [ -z "$(ls -A "$(dirname $out_dir)")" ]
                        then
                           rm -R "$(dirname $out_dir)"
                     fi 
                     continue               		
               fi
            else 
              echo "Folder $d does not contain the \"peaks\" subfolder. Check the JAMM log file for more details. Skipping..";
	      move_folder=$no_peaks_folder/$data_type
              mkdir -p $move_folder
              cp -r $d $move_folder
              ((no_peaks_counter++))
              rm -R $out_dir
              #remove the parent directory also if empty 
              if [ -z "$(ls -A "$(dirname $out_dir)")" ]
                  then
                     rm -R "$(dirname $out_dir)"
              fi 
              continue
          fi
  elif [ $data_type = "MACS" ]
     then
       # echo "I'm in MACS";
        peaks_file=$d/${d%/*}_peaks.narrowPeak
        if [ -e $peaks_file ]
           then
              ################ 1. Add offset to chrStart #################
              awk -v OFS='\t' '{print $1,$2+$10,$2+$10+1,$4,$5,$6}' $peaks_file > $out_file
              #hg_chr_file=$hg38_chr
              #hg_gen_file=$hg38_gen 
           else
              echo "Folder $d does not contain a results file. Check the MACS log file for more details. Skipping..";
              move_folder=$no_peaks_folder/$data_type
              mkdir -p $move_folder
              cp -r $d $move_folder
              ((no_peaks_counter++))
              rm -R $out_dir
              #remove the parent directory also if empty 
              if [ -z "$(ls -A "$(dirname $out_dir)")" ]
                 then
                    rm -R "$(dirname $out_dir)"
              fi 
              continue
        fi
  else
     echo "No supported data_type passed as argument. Only JAMM or MACS supported.";
     exit 1
  fi

   ################ 2. Retrieve matching PFM(s) #################
   # Retrieve the PFM(s) identifier(s) for this TF based on mapping to JASPAR profiles
   # Columns of the mapping file should follow this order:
   # 1. TF_NAME: official name of the transcription factor
   # 2. SYNONYM: all synonyms for this TF as found at TF_checkpoint
   # 3. DNA_BINDING: 'yes' if this TF is known to bind the DNA, 'no' if known not to, 'NA' if unknown
   # 4. DATA_SET_TF_NAME: the TF name as the result of parsing the data set folder name (previous step)
   # 5. JASPAR_TF_NAME: the TF name as it appears in the JASPAR database
   # 6. PFM_ID: the JASPAR identifier for the PFM corresponding to this TF
   # 7. VERSION: the version of the PFM for this TF (should always be the last/newest one used)
   # 8. JASPAR_ID: the internal JASPAR database identifier (i.e., primary key)
   # 9. COLLECTION: the name of the JASPAR collection in which this PFM is present (Note: that it is advised to use CORE)
   # 10.CELL_LINE: the cell line for this ChIP-seq data set  
   # 11.DATA_SET: the name of the ChIP-seq data set
   #
   # Note that not all the above columns are mandatory. The minimum columns needed to ensure the correct functioning
   # of the present script are: 1)DATA_SET_TF_NAME (ENCODE or GEO depending on the data source); 2)PFM_ID; 3)VERSION; 4)CELL_LINE and 5)DATA_SET 
   pfms=()
   noPFM=true
   pfm_id=;
   pfm_version=;
   cell_line_map=;
   while IFS=$'\t' read -r -a line
     do
       data_set_map="${line[10]}"
       data_set_map="$(echo "$data_set_map" | tr '[:lower:]' '[:upper:]')"

       # speed up the file parsing a bit
       if [[ "$data_set_map" < "$data_set" ]]
           then
               continue
       elif [[ "$data_set_map" == "$data_set" ]]
           then
               data_set_tf_name="${line[0]}"
               pfm_id="${line[5]}"
               pfm_version="${line[6]}"
               cell_line_map="${line[9]}"
               
               #to upper case
               data_set_tf_name="$(echo "$data_set_tf_name" | tr '[:lower:]' '[:upper:]')"
               cell_line_map="$(echo "$cell_line_map" | tr '[:lower:]' '[:upper:]')"
               data_set_map="$(echo "$data_set_map" | tr '[:lower:]' '[:upper:]')"
       elif [[ "$data_set_map" > "$data_set" ]]
           then
               break
       fi

       # check if we got a match on data set, cell line and TF
       if [ "$data_set" == "$data_set_map" ] && [ "$cell_line" == "$cell_line_map" ] && [ "$tf_name" == "$data_set_tf_name" ]
          then
              #is there a PFM for this TF
              if [ "$pfm_id" == "NA" ]
                 then
                    echo "No PFM for the transcription factor $data_set_tf_name in data set $data_set was found. Skipping..";
                    ((no_pfm_counter++))
                    rm -R $out_dir
                    # remove parent directory also if empty
                    if [ -z "$(ls -A "$(dirname $out_dir)")" ]
                        then
                           rm -R "$(dirname $out_dir)"
                    fi 
                    break
                 else
                    pfms+=($pfm_folder/${pfm_id}.${pfm_version}.pfm)
                    noPFM=false
#                   echo "PFM selected: $pfm_folder/${pfm_id}.${pfm_version}.pfm"
              fi
          else
            echo "Something does not match:"
            echo "Dataset from file name: $data_set" 
            echo "Dataset from map: $data_set_map" 
            echo "Cell line from file name: $cell_line" 
            echo "Cell line from map: $cell_line_map" 
            echo "TF from file name: $tf_name"
            echo "TF from map: $data_set_tf_name" 
            echo "---------------------------------------"
       fi
   done < $tf_pfm_map

   # no mapping for this dataset - log and skip the folder
   if [[ -z "${pfm_id// }" ]]
      then
        echo "The data set $data_set has no entry in the TF to PWM mapping table. Logging and skipping.." 
        ((no_mapping_counter++))
        rm -R $out_dir
        # remove parent directory also if empty
        if [ -z "$(ls -A "$(dirname $out_dir)")" ]
           then
              rm -R "$(dirname $out_dir)"
        fi
        # add to log if not already in the list
        elementNotPresent "${d%/*}" "${log_no_mapping[@]}"
        notMapped=$?
        if [ $notMapped -eq 0 ] 
          then
            log_no_mapping+=(${d%/*})
            continue
        fi
   fi

   # no PFM(s) found - skip to next folder
   if $noPFM
      then
         continue
   fi

  ################# 3. Increase area around peak max with 'window_size' upstream and downstream ################
  /lsc/bedtoolsolder225/bin/bedtools slop -i $out_file -g $hg_chr_file -b $window_size > ${out_file}_${window_size}bp
  # create TFFM training file: increase area around peak max with 50BP
  /lsc/bedtoolsolder225/bin/bedtools slop -i $out_file -g $hg_chr_file -b $TRAINING_WINDOW_SIZE > ${out_file}_${TRAINING_WINDOW_SIZE}bp
  training_file=${out_file}_${TRAINING_WINDOW_SIZE}bp
  out_file=${out_file}_${window_size}bp

  ################# 4. Map back to the genome and generate a FASTA file ##################
  /lsc/bedtoolsolder225/bin/fastaFromBed -fo ${out_file}.fa -fi $hg_gen_file -bed $out_file
  out_file=${out_file}.fa
  /lsc/bedtoolsolder225/bin/fastaFromBed -fo ${training_file}.fa -fi $hg_gen_file -bed $training_file
  training_file=${training_file}.fa

  #launch processing for each of the PWMs found for the TF in the current data set 
  j=0 #keep track of how many PWMs for this TF
  out_file_bkp=$out_file #keep a copy of the output file name created so far
  for pfm in "${pfms[@]}"
    do
       ((j++))
       if [[ ! -z $pfm ]]
         then
            echo "Using the PfM: $(basename $pfm) for folder ${d#/*} with the transcription factor $tf_name";
            out_file=$out_file_bkp 
          else
            continue
       fi

       #retrieve the current PWM identifier to be added to the output file
       pfm_id=$(basename $pfm)
       pfm_id=${pfm_id%.*}
 
      ################# 5. Scoring the sequences and output ###################
      python $tffm_module --train $training_file --input $out_file --pfm $pfm --outdir $out_dir --tf $tf_name
      out_file=${out_file}_${pfm_id}.tffm
      # Output format:
      # 1. Chromosome:start-end (end-start is SLOP (window_size*2)+1
      # 2. Start offset of the TFBS in the sequence from the window start (Note it is one based, so no need to +1)
      # 3. End of TFBS relative to start offset (2nd column - 3rd column gives the TFBS sequence length)
      # 4. Strand
      # 5. Sequence (Note: that it is 1 based, so there is no need to +1 the start offset of the TFBS 
      # 6. TF name - default TFFM
      # 7. HMM state
      # 8. Score - from 0 to 1 (1 is highest)

      ################# 6. Format file for centrimo and visualization tools ####################
      awk -v OFS='\t' '{print $1,$4,$5,$6,(($3-$2)+1),(($2+$3)/2)-'$window_size',$8,$7,$2,$3}' $out_file | awk -v OFS='\t' -F":" '$1=$1' | awk -v OFS='\t' '{sub(/\-/,"\t",$2)};1' > ${out_file}.ext
      # Output format:
      # 1. Chromosome
      # 2. Start of the TFBS on chromosome (chrStart from step 4 + column 7 from step 4)
      # 3. End of the TFBS on chromosome (resulting chrStart + length of TFBS (8th column - 7th column from step 4)
      # 4. Strand
      # 5. Sequence
      # 6. TF name
      # 7. TFBS sequence length (Note that +1 because it is zero based coming out of the TFFM code)
      # 8. Distance to peak max (!!Note: that this is relative to the chrStart and chrEnd info coming from SLOP)
      # 9. Sequence score - based on probability coming out of TFFM
      # 10. State of the HMM
      # 11. Start offset of the TFBS relative to window start
      # 12. End of the TFBS relative to window start

      #################### 7. Run the thresholding R script ##################
      # It will dynamically calculate score and enrichment zone thresholds, and centrality p-value
      # based on the centrimo concept. 
      # NOTE: you do not need all the sequences per peak, only the top scoring one. Discussed twice, 
      # looked in the original centrimo paper (twice), it suffices to pass only top scoring sequences 
      # (one per peak) and the total number of sequences scored (also one per peak) - a.k.a output
      #  of C score with -b is ok.. :]
      #
      # Parameters:
      #  1. Path to script wrapper for threshold calculation
      #  2. Input file containing the scored sequences (formatted file from step 5)
      #  3. Type of score (relative, prob or abs) - same as score.type parameter in compute_thresholds.R
      #  4. The column holding the TFBS sequence length - same as seq.length.column in compute_thresholds.R
      #  5. The column holding the distance to the peak summit - same as distance.column parameter in compute_thresholds.R
      #  6. The column containing the score - same as score.column parameter in compute_thresholds.R
      #  7. The size of the window used in bedtools slop
      #  8. The name of the output file - same as output.file.name parameter in compute_threshold.R
      thresholds=$(Rscript $get_thresholds $compute_thresholds ${out_file}.ext "prob" 7 8 9 $window_size ${out_file}.top)
      
      if [ "$thresholds" != "OK" ]
         then
            echo "No thresholds calculated for $d. Maybe signal is too poor :/"
            move_folder=$no_thresholds_folder/$d/$data_type
            mkdir -p $move_folder
            # more than one PFM for this
            if [ "$j" -gt "1" ] || [ "${#pwms[@]}" -gt "1" ]
                then
                   #move only files generated for the current PWM and not all
                   mv $out_dir/*$pwm_id* $move_folder
                #signal too poor and we reached end of array then remove all files
                if [ "$j" -eq "${#pwms[@]}" ]
                   then
                      mv $out_dir/* $move_folder
                      rm -R $out_dir
                      #remove the parent directory also if empty 
                      if [ -z "$(ls -A "$(dirname $out_dir)")" ]
                         then
                             rm -R "$(dirname $out_dir)"
                      fi 
                fi
            fi
            ((no_thresholds_counter++))
            continue
      fi
   done #end for over PFMs
done #end for over folders

#output data set names with no TF to PWM mapping entry
if [ "$no_mapping_counter" -gt "0" ]

  then
     printf '%b\n' "${log_no_mapping[@]}" > $log_no_mapping_filename
fi

#ouput processing summary
echo "--------------------------------------------------------------------------------------"
echo "Batch: $(basename $in_dir)"
echo "For data type $data_type there were:"
echo "$no_pfm_counter data set(s) with no PWM(s)"
echo "$no_peaks_counter data set(s) with no peaks found during processing"
echo "$no_thresholds_counter data set(s) with poor signal"
echo "$no_mapping_counter data set(s) with no TF to PWM mapping entry"
echo "--------------------------------------------------------------------------------------"
