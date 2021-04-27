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

This script can be used to create the an Amazon Machine Image (AMI) that has all the required dependencies to run the bioinformatics pipeline. This is achieved by launcing an instance in AWS using a bash script on launch to configure the instance which is then used to create the AMI within AWS. Furthermore, the script also creates a keypair for EC2 access, sets up a security group for the instances and creates an EC2 IAM role that allows EC2 instances to access S3 services. 

#### Input:
- Config file for programmatic access to AWS.
- The bash_script.sh that is used as a bootstrap script when launching the instance. 

#### To run:

> python create_ami.py

#### Output:
- A config file (ec2_config.yml) with the keypair name, the security group ID, the AMI ID and the IAM role name 
- Keypair for EC2 instance access
- An AMI that can be used to launch instances with the required software and packages for running the pipeline 

### run_pipeline.py

#### Input:
- Config file for programmatic access to AWS
- ec2_key.pem: Key to allow instance access for troubleshooting purposes, created by create_ami.py
- ec2_config.yml: This config file, created by create_ami.py, contains the relevant AWS information to run the instance that will be used to run the pipeline, such as the AMI that should be used to launch the instance

#### To run:

> python run_pipeline.py

The S3 bucket names will need to be changed within this script (lines ) if run using different S3 buckets. 

#### Output:
- All pipeline outputs can be found in the output S3 bucket. 

#### Future Work
1. Look into how to cost optimise pipelines run this way - what changes may be helpful for this? 
