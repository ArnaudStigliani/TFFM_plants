import sys
sys.path.append("../lib/TFFM")
import os
import tffm_module
from constants import TFFM_KIND

from Bio import motifs 
from argparse import ArgumentParser


# argument parser
parser = ArgumentParser()
parser = ArgumentParser(description='TFFM arguments')
parser.add_argument('--train', dest='train', type=str,
                    help='the path to the input FASTA file for training the TFFM')
parser.add_argument('--input', dest='data', type=str,
                    help='the path to the input FASTA file to apply the TFFM on')
parser.add_argument('--pfm', dest='pfm', type=str,                    
                    help='the path to the PFM to be used')
parser.add_argument('--outdir', dest='out', type=str,
		    help='the path to the output folder')
parser.add_argument('--tf', dest='tf', type=str,
		    help='the name of the transcription factor')		
arguments = parser.parse_args()

with open(arguments.pfm) as handle:
	motif = motifs.read(handle, "pfm") 

output_file=arguments.out + "/" + os.path.basename(arguments.data) + "_" + os.path.splitext(os.path.basename(arguments.pfm))[0]
#print(output_file)

#will use the same data set for training and applying
tffm_first_order = tffm_module.tffm_from_motif(motif, TFFM_KIND.FIRST_ORDER, arguments.tf)
tffm_first_order.train(arguments.train)
#output tffm model
tffm_first_order.write(output_file + "_tffm.xml")

out = open(output_file + "_summary_logo.svg", "w")
tffm_first_order.print_summary_logo(out)
out.close()
out = open(output_file + "_dense_logo.svg", "w")
tffm_first_order.print_dense_logo(out)
out.close()

# the detailed model
tffm_detailed = tffm_module.tffm_from_motif(motif, TFFM_KIND.DETAILED, arguments.tf)
tffm_detailed.train(arguments.train)
tffm_detailed.write(output_file + "tffm_detailed.xml")
out = open(output_file + "_detailed_summary_logo.svg", "w")
tffm_detailed.print_summary_logo(out)
out.close()
out = open(output_file + "_detailed_dense_logo.svg", "w")
tffm_detailed.print_dense_logo(out)
out.close()




out = open(output_file + ".tffm", "w")
for hit in tffm_first_order.scan_sequences(arguments.data, only_best=True):
    #print only hits and not lines for which the chains in the HMM got stuck
    hit = str(hit)
    if not "None" in hit:     
        out.write(hit + '\n')
out.close()

out = open(output_file + ".tffm.detailed", "w")
for hit in tffm_detailed.scan_sequences(arguments.data, only_best=True):
    #print only hits and not lines for which the chains in the HMM got stuck
    hit = str(hit)
    if not "None" in hit:
        out.write(hit + '\n')
out.close()

