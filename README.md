# AWS Use Case 1 - EC2

This code was used to run the bioinformatics pipeline found at https://github.com/Becky-Sadler/aws-ec2-pipeline on AWS using EC2. This code was created for the completion of an MSc project.

## Prerequisites
AWS programmatic access keys are required. These should be stored within a config file, in the format in the template config file.

Both an input and output S3 bucket should be created. The FASTQ files to be run should be present in an input bucket within AWS S3. 

## Installation
1. Clone the repo
> git clone https://github.com/Becky-Sadler/DECIPHER_upload

2. Install requirements.txt in a virtualenv of your choice
> pip3 install -r requirements.txt

## Usage

### create_ami.py

#### Input:


#### To run:

> python create_ami.py


#### Output:


### run_pipeline.py

#### Input:


#### To run:

> python run_pipeline.py

#### Output:


#### Future Work
