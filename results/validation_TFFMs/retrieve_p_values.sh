#!/bin/bash

path=$(pwd)

cd ../processed
paste <(ls */* | grep tffm.ext.centrimo.pval | sed 's/\/.*//'  ) <(ls */* | grep tffm.ext.centrimo.pval | xargs -I{} cat {}  ) | sort -k2,2n > $path/TF_pval.csv

exit 0 
