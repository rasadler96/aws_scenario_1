# AWS Use Case 1 - EC2

This code was used to run the bioinformatics pipeline found at https://github.com/Becky-Sadler/aws-ec2-pipeline on AWS using EC2. This code was created for the completion of an MSc project.

## Prerequisites

AWS programmatic access keys are required. These should be stored within a config file, called config.yml, in the format in the template config file (template_config.txt).

Both an input and output S3 bucket should be created within AWS S3. The FASTQ files to be processed should be present in the input bucket. 

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

This script is used to launch an instance and run the bioinformatics pipelines for the fastq files present in the input bucket. The instance is launched using the information held within the ec2_config.yml file created by the create_ami.py script. The python library paramiko is then used to pass 3 commands to the instance via SSH.

1. Pulling the input files from the input S3 bucket into the instance
2. Running the pipeline
3. Pushing all the results files (including log files) to the output S3 bucket, excluding the FASTQ input files. 

#### Input:

- Config file for programmatic access to AWS
- ec2_key.pem: Key to allow instance access for troubleshooting purposes, created by create_ami.py
- ec2_config.yml: This config file, created by create_ami.py, contains the relevant AWS information to run the instance that will be used to run the pipeline, such as the AMI that should be used to launch the instance

#### To run:

> python run_pipeline.py

The S3 bucket names will need to be manually changed within this script (lines 170 & 172) to run with other buckets. 

#### Output:
- All pipeline outputs can be found in the output S3 bucket. This includes both the final VCF files, any intermediate files and the log files for each process. 

#### Future Work
1. Look into how to cost optimise pipelines run this way - what changes may be helpful for this? 
