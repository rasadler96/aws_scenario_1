#!/bin/bash

# installing the required packages  

apt-get update

apt install -y curl wget gcc libssl-dev make unzip openjdk-8-jdk git-all bwa samtools virtualenv python3.6-dev default-libmysqlclient-dev mysql-server

cd /opt 
mkdir software 
cd software 

# Installing GATK toolkit (Picard tools are bundled in with this) 
wget https://github.com/broadinstitute/gatk/releases/download/4.1.6.0/gatk-4.1.6.0.zip
unzip gatk-4.1.6.0.zip
rm -r gatk-4.1.6.0.zip

# Installing aws-cli
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install

# Pulling the code into the instance 
cd /home/ubuntu
git clone https://github.com/Becky-Sadler/aws-ec2-pipeline.git
cd aws-ec2-pipeline
mkdir run-directory reference
mv FH.bed reference 

# Creating virtual environment

virtualenv -q -p python3 pipeline_env
pipeline_env/bin/pip install -r aws-requirements.txt

# Sorting out reference sequence
cd reference

wget https://1000genomes.s3.amazonaws.com/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz 

# creating .dict
/opt/software/gatk-4.1.6.0/gatk CreateSequenceDictionary -R hs37d5.fa

# creating .fa.fai
samtools faidx hs37d5.fa

# index for use with BWA

bwa index hs37d5.fa

# stop instance

shutdown now -h


