
#called by the results_processing.sh script

args = commandArgs(trailingOnly = T)
if (length(args) != 0){
  source(args[1])
  
  #------------------------------------------------------------------------
  #     Arguments: 
  # 1. Path to the threshold calculation script (compute_thresholds.R)
  # 2. Input data file - containing the extended score file (.score.ext resulting from the result_processing.sh script)
  # 3. Type of score to be taken into consideration (relative or absolute)
  # 4. The column holding the distance to the peak max
  # 5. The column holding the score information
  # 6. Name of the output file
  #-----------------------------------------------------------------------

    #debug print
#    print("The threshold calculation script:")
#    print(args[1])
#    print("Data file:")
#    print(args[2])
#    print("Score type:")
#    print(args[3])
#    print("Sequence length column:")
#    print(args[4]) 
#    print("Distance column:")
#    print(args[5])
#    print("Score column:")
#    print(args[6])
#    print("Window size:")
#    print(args[7])
#    print("Output file name:")
#    print(args[8])
 
    #call the threshold calculation
#    suppressMessages(require("rJava" , lib.loc="/storage/scratch/marius/ENCODE/src/"))
    compute_thresholds(data.file = args[2], score.type = args[3], seq.length.column = args[4], distance.column = args[5], score.column = args[6], window.size = args[7], output.file.name = args[8])

}  
