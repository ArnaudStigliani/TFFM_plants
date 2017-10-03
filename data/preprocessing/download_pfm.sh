#!/bin/bash

mkdir -p jaspar_pfm
cd jaspar_pfm
while read -r line
do 
    wget http://hfaistos.uio.no:8000/download/JASPAR_CORE/pfm/individual/$line.pfm
done <../list_pfm

sed -i -e '1d' -e 's/[ACGT]\s*\[\(.*\)]/\1/' * 
cd ..
mv jaspar_pfm ..

exit 0
