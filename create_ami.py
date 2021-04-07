import boto3
import yaml
import os
import botocore
import json

# Function to create key pair to be used to launch instance - input = key pair name (minus the .pem), output = keypair.pem with correct permissions. 
def create_keypair(name_of_keypair):
	key_file = open('%s.pem'%name_of_keypair,'w')
	try:
		key = ec2_resource.create_key_pair(KeyName=name_of_keypair)
		key_pair_contents = str(key.key_material)
		key_file.write(key_pair_contents)
		os.system('chmod 400 %s.pem'%name_of_keypair)
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		print('Key pair %s.pem sucessfully created'%name_of_keypair)
		return(name_of_keypair)

def create_security_group(description, name):
	try:
		response = ec2_client.create_security_group(
			Description = description,
			GroupName = name
		)
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		print('Security group sucessfully created')	
		sg_id = response['GroupId']
		return(sg_id)

def create_sg_rule(groupid, ipPermissions):
	try:	
		response = ec2_client.authorize_security_group_ingress(
	    GroupId= groupid,
	    #GroupName='string',
	    IpPermissions= ipPermissions
	)
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		print('Security group rule added: %s'%ipPermissions)	

def create_instances(**kwargs):
	try:
		response = ec2_client.run_instances(**kwargs)
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		instance_ID = response['Instances'][0]['InstanceId']
		print('Instance created (%s)'%instance_ID)
		return instance_ID	

def get_user_data(file_name):
    f = open(file_name, 'r')
    user_data = f.read()
    return user_data

# Waiter definition **kwargs must be in the form of a dictionary. 
def add_waiter(waiter_type, **kwargs):
	try:
		waiter = ec2_client.get_waiter(waiter_type)
		waiter.wait(**kwargs)
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		print(waiter_type)

# Create AMI function **kwargs must be in the form of a dictionary. 
def create_ami(**kwargs):
	try:
		response = ec2_client.create_image(**kwargs)
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		ami_ID = response['ImageId'] 
		print('AMI created: %s'%ami_ID)
		return ami_ID

# Create IAM role
def create_iam_role(**kwargs):
	try:
		response = iam_client.create_role(**kwargs)
	except botocore.exceptions.ClientError as e:
		print(e)
	else:
		role_name = response['Role']['RoleName']
		print('IAM role created')
		return role_name

# Add policy to IAM rule
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

# Create instance profile (Needed to attach role to an instance)
def create_instance_profile(instance_profile_name):
	try:
		iam_client.create_instance_profile(
			InstanceProfileName = instance_profile_name
		)
	except botocore.exceptions.ClientError as e:
		print(e)
	else:
		print('Instance profile %s created'%instance_profile_name)

# Add role to instance profile
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

# clean up function that terminates the instance and puts the required ec2 values (key pair name, security group ID and AMI id) into a ec2 config file
def clean_up(instanceID, keypairName, sgID, amiID, rolename):
	try:
		ec2_resource.instances.filter(InstanceIds=[instanceID]).terminate()
	except botocore.exceptions.ClientError as e:
		print(e)
	else: 
		print('Instance %s terminated'%instanceID)
		#config here 
		data = {'ec2_information': {'keypair_name': str(keypairName), 'security_group_ID': str(sgID), 'ami_ID' : str(amiID), 'iam_role_name' : str(rolename)}}
		config_file = open('ec2_config.yml', 'w')
		yaml.dump(data, config_file)
		print('ec2_config file created')


# Pulling amazon credentials from the config file. 
amazon_config = yaml.safe_load(open("config.yml"))

amazon_details = amazon_config['amazon']
 
access_key = amazon_details['aws_access_key_id']
secret_key = amazon_details['aws_secret_access_key']
default_region = amazon_details['aws_default_region']

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

# Creating AMI to run pipeline 

# Create ec2 key pair
key_name = create_keypair('ec2_key')

# Create EC2 security group
security_group_id = create_security_group('Security group for EC2 Scenario 1', 'EC2 group')

# Defining a security group rule - this allows SSH access to the instance 
ipPermissions =[
        {
            'FromPort': 22,
            'IpProtocol': 'tcp',
            'IpRanges': [
                {
                    'CidrIp': '0.0.0.0/0',
                    'Description': 'SSH access',
                },
            ],
            'ToPort': 22,
        }
    ]

# Adding rule to the security group 
create_sg_rule(security_group_id, ipPermissions)

# Pulling the bootstrap script into a string to pass to the make instance function. 
bootstrap_script = get_user_data('bash_script.sh')

# Creating variable with all the instance details - must be in the form of a dictionary. AMI_ID used is the official Ubuntu ami
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
'ImageId' : 'ami-0eb89db7593b5d434',
'InstanceType' : 't2.2xlarge',
'KeyName' : key_name,
'MinCount' : 1,
'MaxCount' : 1,
'SecurityGroupIds' : [
    security_group_id,
],
'UserData' : bootstrap_script}

# Creating instance
instance_id = create_instances(**instance_details)

# Specifying waiter details - longer delay (10 minutes) for the stopping waiter compared to the running(20 seconds)

waiter_run = {'InstanceIds' : [
    instance_id,
],
'WaiterConfig' : {
    'Delay': 20,
    'MaxAttempts': 100
}}

waiter_stop = {'InstanceIds' : [
    instance_id,
],
'WaiterConfig' : {
    'Delay': 600,
    'MaxAttempts': 100
}}

# Don't forget to add ** infront of kwargs argument! 
add_waiter('instance_running', **waiter_run)

add_waiter('instance_stopped', **waiter_stop)

ami_details = {
'BlockDeviceMappings' : [
    {
        'DeviceName': '/dev/sda1',
        'Ebs': {
            'DeleteOnTermination': True,
            'VolumeSize': 15,
            'VolumeType': 'gp2',
        },
    },
],
'Description': 'Ec2 AMI',
'InstanceId' : instance_id,
'Name' : 'EC2_AMI',
}

ami_id = create_ami(**ami_details)

# Creating an IAM role for EC2 to access S3 

# Create a trust permission (giving EC2 to ability to take on the role created)
ec2_role_access = {
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

# Defining role details for attaching to later instances 
role_details = {
'RoleName':'EC2_S3_Access',
'AssumeRolePolicyDocument' : json.dumps(ec2_role_access),
'Description':'Role to give EC2 access to S3',
'MaxSessionDuration' : 43200}

# Creating the role 
role_name = create_iam_role(**role_details)

# Adding the S3 full access policy (EC2 can fully access S3)
add_policy('arn:aws:iam::aws:policy/AmazonS3FullAccess', role_name)

# Create instance profile and add role -> Name of instance profile == same as role name (makes it easier and is how this occurs if done through the console)

create_instance_profile(role_name)
add_role_to_instance_profile(role_name, role_name)

clean_up(instance_id, key_name, security_group_id, ami_id, role_name)



