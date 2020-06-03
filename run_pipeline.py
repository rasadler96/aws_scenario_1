import boto3 
import yaml 
import botocore 
import os 
import json 
import paramiko
import time

# Config for credentials and ec2 information
amazon_config = yaml.safe_load(open("config.yml"))

amazon_details = amazon_config['amazon'] 
access_key = amazon_details['aws_access_key_id']
secret_key = amazon_details['aws_secret_access_key']
default_region = amazon_details['aws_default_region']

ec2_config = yaml.safe_load(open("ec2_config.yml"))

ec2_details = ec2_config['ec2_information']
ami_id = ec2_details['ami_ID']
keypair_name = ec2_details['keypair_name']
security_group_id = ec2_details['security_group_ID']
role_name = ec2_details['iam_role_name']

# Creating session (configures credentials and default region)
session = boto3.Session(
	aws_access_key_id = access_key,
	aws_secret_access_key = secret_key,
	region_name = default_region
)

# Creating ec2 resource and client for session
ec2_resource = session.resource('ec2', region_name=default_region)
ec2_client = session.client('ec2', region_name=default_region)

# Function to check the state of the ami image. 
def image_state(amiID):
	try: 
		response = ec2_client.describe_images(ImageIds=[amiID])
	except botocore.exceptions.ClientError as e: 
		print(e)
	else: 
		state = response['Images'][0]['State']
		return state

def create_instances(**kwargs):
	try:
		response = ec2_client.run_instances(**kwargs)
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		instance_ID = response['Instances'][0]['InstanceId']
		print('Instance created (%s)'%instance_ID)
		return instance_ID

# Waiter definition **kwargs must be in the form of a dictionary. 
def add_waiter(waiter_type, **kwargs):
	try:
		waiter = ec2_client.get_waiter(waiter_type)
		waiter.wait(**kwargs)
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		print(waiter_type)

# Function to get the public DNS of the instance

def get_public_ip(instance_id):
	try:
		response = ec2_client.describe_instances(InstanceIds=[instance_id])
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		public_DNS = response['Reservations'][0]['Instances'][0]['PublicIpAddress']
		return public_DNS

def run_pipeline(keypair_name, public_ip, commands):
	key = paramiko.RSAKey.from_private_key_file('%s.pem'%keypair_name)

	ssh = paramiko.SSHClient()

	ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	
	ssh.connect(hostname=public_ip, username="ubuntu", pkey=key)

	# executing list of commands within server
	print("Starting execution")
	for command in commands:
		print("Executing command: " + command)
		stdin , stdout, stderr = ssh.exec_command(command)
		exit_status = stdout.channel.recv_exit_status()  
		if exit_status == 0:
			print(stdout.readlines())
			print(stderr.readlines())
		else:
			print("Error", exit_status)
			print('Command failed to execute: ' + command)

	ssh.close()

def terminate_instance(instanceID):
	try:
		ec2_resource.instances.filter(InstanceIds=[instanceID]).terminate()
	except botocore.exceptions.ClientError as e:
		print(e)
	else:
		print('Instance terminated')

# Running script from here (Above are the functions)

# Defining variables for instance
instance_details = {'BlockDeviceMappings' : [
    {
        'DeviceName' : '/dev/sda1',
        'Ebs': {
            'DeleteOnTermination': True,
            'VolumeSize': 15,
            'VolumeType': 'gp2',
            'Encrypted': False
        },
    },
],
'ImageId' : ami_id,
'InstanceType' : 't2.2xlarge',
'KeyName' : keypair_name,
'MinCount' : 1,
'MaxCount' : 1,
'SecurityGroupIds' : [
    security_group_id,
],
'IamInstanceProfile': {
	'Name' : role_name,
}}


# Checking AMI state
state = image_state(ami_id) 

# Making sure that the AMI is available (okay) before launching an instance and running the pipeline. 
if state == 'available':
	print('AMI is available')
	instance_id = create_instances(**instance_details)
	waiter_run = {'InstanceIds' : [
    instance_id,
],
'WaiterConfig' : {
    'Delay': 20,
    'MaxAttempts': 100
}}
	add_waiter('instance_running', **waiter_run)
	
	# Getting public DNS name
	public_ip = get_public_ip(instance_id)

	# Need to ensure that checks are done (can add a waiter for this)
	waiter_status_ok = {'InstanceIds' : [
    instance_id,
],
'WaiterConfig' : {
    'Delay': 20,
    'MaxAttempts': 100
}}
	# waiter to make sure the instance is ok to ssh into 
	add_waiter('instance_status_ok', **waiter_status_ok)

	print('SSH into instance')

	# command list
	commands = [
		"sudo aws s3 sync s3://giab-fastq-files/run-directory /home/ubuntu/aws-ec2-pipeline/run-directory",
		"cd aws-ec2-pipeline ; sudo /home/ubuntu/aws-ec2-pipeline/pipeline_env/bin/python aws-pipeline.py --input /home/ubuntu/aws-ec2-pipeline/run-directory/ -j 2 --verbose 5",
		"aws s3 sync ~/aws-ec2-pipeline/run-directory/ s3://pipeline-output-files/ --exclude='*.fastq.gz'",
	]

	# SSH into instance and run the commands
	run_pipeline(keypair_name, public_ip, commands)

	print('Pipeline run and files transferred to pipeline-output-files S3 bucket')

else:
	print('There is a problem with the selected AMI - state is "' + str(state) + '"') 

# Terminate the instance as no longer required. 
terminate_instance(instance_id)

