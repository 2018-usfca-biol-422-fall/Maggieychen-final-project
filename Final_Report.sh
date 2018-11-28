#!/bin/bash

source /home/.profile

# This is script to download the sequence data set
# and save them into a proper folders.

# Maggie Chen
# November 26, 2018
# ychen254@dons.usfca.edu


# Downloads the zip file
 echo "Downloading the zip file to data/raw_data directory...."
 curl -L https://www.mothur.org/MicrobiomeBiomarkerCRC/Zackular_EDRN_fastq_files.gz.tar data/raw_data/fasta_archive.zip



# In order to use this script to get information about fasta files
# please unzip them, thus they are in the form of .fasta
# create a directory to save the unziped data 
unzip data//raw_data/fasta_archive.zip -d data/raw_data
rm data//raw_data/fasta_archive.zip

