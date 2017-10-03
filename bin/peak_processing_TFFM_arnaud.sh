#!/bin/bash
#usage: bash peak_processing_TFFM.sh JAMM ../peaks ../pfm ../doc/map.csv 200
workspace="../results"

#parse arguments
data_type=;
if [ -z "$1" ]
    then
        echo "No data type specified. Using OMalley...";
        data_type="OMalley"
    else
        data_type=$1
fi

in_dir=
if [ -z "$2" ]
    then
        echo "No input file folder specified. Using ../data/narrowPeaks ";
        in_dir="../data/narrowPeaks"
    else
        in_dir=$2
fi

pfm_folder=;
if [ -z "$3" ]
    then
        echo "No PFM folder specified. Using: ../data/pfm";
        pfm_folder="../data/jaspar_pfm"
    else
        pfm_folder=$3
fi

tf_pfm_map=;
if [ -z "$4" ]
    then
        echo "No TF to PFM mapping file specified. Using: ../data/map.csv";
        tf_pfm_map="../data/map.csv"
    else
        tf_pfm_map=$4
fi


window_size=;
if [ -z "$5" ]
    then
        echo "No window size specified. Using 250bp";
        window_size=250
    else
        window_size=$5
fi


########## PARAMETER SETTING  ##################
result_folder="../results"
# faulty data sets
no_peaks_folder="../results/no_peak"
no_mapping_folder="../results/no_mapping"
no_thresholds_folder="../results/no_thresholds"
mkdir -p $no_thresholds_folder
# script locations
tffm_module="./TFFM_first_order.py"
compute_thresholds="./compute_thresholds_centrimo.R"
get_thresholds="./get_thresholds.R"

hg_chr_file="../../data/data_arnaud/tair10.size"
hg_gen_file="../../data/data_arnaud/tair10.fas"
nb_processes=40

TRAINING_WINDOW_SIZE=50 #the size of the window on which the TFFM should train

############ END OF PARAMETER SETTING #########

#print script arguments
#echo "------------------------------------------------------"
#echo "Data type: $data_type"
#echo "Window size: $window_size"
#echo "Input folder: $in_dir"
#echo "PFM folder: $pfm_folder"
#echo "Workspace: $workspace"
#echo "------------------------------------------------------"



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

# #function to check if an element is NOT in an array
# elementNotPresent () {
#   local i
#   for i in "${@:2}"; do [[ "$i" == "$1" ]] && return 1; done
#   return 0
# }

#sort the TF to PFM mapping file by data set (change column if not the 11th)
if [ ! -e ${tf_pfm_map}.sorted ]
   then
      sort -t $'\t' -k 11,11 $tf_pfm_map > ${tf_pfm_map}.sorted
fi
tf_pfm_map=${tf_pfm_map}.sorted



# get ready to process
i=1
nb_peak_files="$(ls -l $in_dir | grep -c narrow)"
echo "There are $nb_peak_files peak files"
#Create a folder to receive the results
out_dir=$result_folder/processed
mkdir -p $out_dir


#retrieving TF names
tf_names=$(ls $in_dir | sed 's/.narrowPeak//')
# into an array
tf_name=(${tf_names// /})


j=0


if [ $data_type = "OMalley" ]
then
    # echo "I'm in OMalley";
    peaks_file=$in_dir/*.narrowPeak
    #into an array
    peaks_file=(${peaks_file// /})
    for k in ${!peaks_file[@]}
    do
	((j++))
	if [ $j == $nb_processes ]
	then 
	    wait
	    echo "new processes"
	    j=0
	fi
	(
	if [ -e ${peaks_file[$k]} ]
	then
	    mkdir -p $out_dir/${tf_name[$k]}
            ################ 1. find the middle of the peak  #################
                   awk -v OFS='\t' '{print $1,$2+100,$2+100+1,$4,$5,$6}' ${peaks_file[$k]} > $out_dir/${tf_name[$k]}/${tf_name[$k]}_peaks
	    # echo $out_dir/${tf_name[$k]} 
	fi
	)&
    done
fi




wait


j=0
sed '1d' $tf_pfm_map | while IFS=$'\t' read -r -a line
do
	tf_map=${line[3]}
	pfm_id_map=${line[5]}
	pfm_version_map=${line[6]}
	out_file=${out_dir}/${tf_map}/${tf_map}_peaks
	pfm_id=$pfm_folder/${pfm_id_map}.${pfm_version_map}.pfm
	tf_name=${line[0]}
	((j++))
	while [ $j -ge $nb_processes ]
	do
	    j=$(pgrep -u $USER python | wc -l)
	    sleep 5 
	done
	echo "new tf : $tf_name"
	(
	# check if the peak file exists, else go to next line
	if [[ ! -e $out_file ]]
	then
	    echo "The  peak_file $out_file does not exist"
	    continue
	fi

	# check if the matrix exists, else go to next line
	if [[ ! -e $pfm_id ]]
	then
	    echo "The pfm $pfm_id does not exist"
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

	### pipeline TFFM ######
	
	python $tffm_module --train $training_file --input $out_file --pfm $pfm_id --outdir $out_dir/${tf_map} --tf $tf_name 
	basename_pfm=$(basename $pfm_id .pfm)
	out_file=${out_file}_${basename_pfm}.tffm
	# Output format:
	# 1. Chromosome:start-end (end-start is SLOP (window_size*2)+1
	# 2. Start offset of the TFBS in the sequence from the window start (Note it is one based, so no need to +1)
	# 3. End of TFBS relative to start offset (2nd column - 3rd column gives the TFBS sequence length)
	# 4. Strand
	# 5. Sequence (Note: that it is 1 based, so there is no need to +1 the start offset of the TFBS 
	# 6. TF name - default TFFM
	# 7. HMM state
	# 8. Score - from 0 to 1 (1 is highest)
	
	awk -v OFS='\t' '{print $1,$4,$5,$6,(($3-$2)+1),(($2+$3)/2)-'$window_size',$8,$7,$2,$3}' $out_file | awk -v OFS='\t' -F":" '$1=$1' | awk -v OFS='\t' '{sub(/\-/,"\t",$2)};1' > ${out_file}.ext
	
	#and for the detailed one
	awk -v OFS='\t' '{print $1,$4,$5,$6,(($3-$2)+1),(($2+$3)/2)-'$window_size',$8,$7,$2,$3}' ${out_file}.detailed | awk -v OFS='\t' -F":" '$1=$1' | awk -v OFS='\t' '{sub(/\-/,"\t",$2)};1' > ${out_file}.detailed.ext

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
	
	#and for the detailed one
	thresholds=$(Rscript $get_thresholds $compute_thresholds ${out_file}.detailed.ext "prob" 7 8 9 $window_size ${out_file}.detailed.top)
	
	if [ "$thresholds" != "OK" ]
	then
            echo "No thresholds calculated for $tf_map. Maybe signal is too poor :/"
	    echo "dir to mv is ${out_dir}/${tf_map}/"
	    mv ${out_dir}/${tf_map}/ $no_thresholds_folder
	fi
	) 2>${out_dir}/${tf_map}/err.log &  
done 
