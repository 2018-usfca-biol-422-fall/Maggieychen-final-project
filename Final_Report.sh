#!/bin/bash

# shellcheck source=/dev/null
source /home/.profile

# This is script to download the sequence data set
# and save them into a proper folders.

# Maggie Chen
# November 26, 2018
# ychen254@dons.usfca.edu


# Downloads the zip file
# wget is used to downloda the data file from link 
# and provide the download state
 echo "Downloading the zip file to data/raw_data directory...."
 wget https://www.mothur.org/MicrobiomeBiomarkerCRC/Zackular_EDRN_fastq_files.gz.tar

# unzip the files and delete the original zip file
# tar to unzip files in the formate of gz.tar
# -x extrac the unziped files
# -v make tar talk a lot, verbose output
# -f must be the last and right before the file name
# tells tar the name and the path of the file
# -C for change the directory to save the unziped files
tar -xvf Zackular_EDRN_fastq_files.gz.tar -C data/raw_data/

# Due to it is a paired-end sequenceing
# this project will only deal with R1 sequqences
# thus, delet the files has R2 in the name rather tan R1
for file in *_R2_*
do 
   rm "data/raw_data/$file" 
done


