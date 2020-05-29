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

# Creating session (configures credentials and default region)
session = boto3.Session(
	aws_access_key_id = access_key,
	aws_secret_access_key = secret_key,
	region_name = default_region
)

# Creating ec2 resource and client for session
ec2_resource = session.resource('ec2', region_name=default_region)
ec2_client = session.client('ec2', region_name=default_region)

# Creating IAM client for session
iam_client = session.client('iam')

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

def create_iam_role(**kwargs):
	try:
		response = iam_client.create_role(**kwargs)
	except botocore.exceptions.ClientError as e:
		print(e)
	else:
		role_name = response['Role']['RoleName']
		print('IAM role created')
		return role_name

def add_policy(policy_arn, role_name):
	try:
		iam_client.attach_role_policy(
    	PolicyArn=policy_arn,
    	RoleName=role_name
		)
	except botocore.exceptions.ClientError as e:
		print(e)
	else:
		print('%s policy added to %s' %(policy_arn, role_name))

def create_instance_profile(instance_profile_name):
	try:
		iam_client.create_instance_profile(
			InstanceProfileName = instance_profile_name
		)
	except botocore.exceptions.ClientError as e:
		print(e)
	else:
		print('Instance profile %s created'%instance_profile_name)

def add_role_to_instance_profile(instance_profile_name, role_name):
	try:
		iam_client.add_role_to_instance_profile(
			InstanceProfileName= instance_profile_name,
			RoleName= role_name
		)
	except botocore.exceptions.ClientError as e:
		print(e)
	else:
		print('Role added to instance profile')

# Running script from here (Above are the functions)

# Creating an IAM role for EC2 to access S3 

# First need to create a trust permission (giving EC2 to ability to take on the role created)
'''ec2_role_access = {
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "",
      "Effect": "Allow",
      "Principal": {
        "Service": "ec2.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}

# Defining role details 
role_details = {
'RoleName':'EC2_S3_Access',
'AssumeRolePolicyDocument' : json.dumps(ec2_role_access),
'Description':'Role to give EC2 access to S3',
'MaxSessionDuration' : 43200}

# Creating the role 
role_name = create_iam_role(**role_details)

# Adding the S3 full access policy (EC2 can fully access S3)
add_policy('arn:aws:iam::aws:policy/AmazonS3FullAccess', role_name)


# Create instance profile and add role (There is a delay with this so for testing I have created role and instance profile once and haed coded the name for creating instance (so I can test the remaining part without having to keep deleting the roles etc... Maybe in future move this step to when ami is created))

create_instance_profile(role_name)
add_role_to_instance_profile(role_name, role_name)'''


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
'InstanceType' : 't2.micro',
'KeyName' : keypair_name,
'MinCount' : 1,
'MaxCount' : 1,
'SecurityGroupIds' : [
    security_group_id,
],
'IamInstanceProfile': {
	'Name' : 'EC2_S3_Access',
}}


# Checking AMI state
state = image_state('ami-01dcf48509ea3e1ff') 

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
	time.sleep(120)

	key = paramiko.RSAKey.from_private_key_file('%s.pem'%keypair_name)

	ssh = paramiko.SSHClient()

	ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	
	ssh.connect(hostname=public_ip, username="ubuntu", pkey=key)

	stdin, stdout, stderr = ssh.exec_command('sudo aws s3 sync s3://giab-fastq-files /home/ubuntu/aws-ec2-pipeline/run-directory')
	exit_status = stdout.channel.recv_exit_status()  

	if exit_status == 0:
		print("Fastq files pulled into run-directory")
	else:
		print("Error", exit_status)

	stdin, stdout, stderr = ssh.exec_command('cd /home/ubuntu/aws-ec2-pipeline ; sudo virtualenv -p python3 pipeline_env ; source pipeline_env/bin/activate ; python aws-pipeline.py --input /home/ubuntu/aws-ec2-pipeline/run_directory/ -j 2 --verbose 5')
	exit_status = stdout.channel.recv_exit_status()  
	
	if exit_status == 0:
		print("Processes done?")
	else:
		print("Error", exit_status)

	ssh.close()

else:
	print('There is a problem with the selected AMI - state is "' + str(state) + '"') 

