import boto3
import yaml
import os
import botocore

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

# clean up function that terminates the instance and puts the required ec2 values (key pair name, security group ID and AMI id) into a ec2 config file
def clean_up(instanceID, keypairName, sgID, amiID):
	try:
		ec2_resource.instances.filter(InstanceIds=[instanceID]).terminate()
	except botocore.exceptions.ClientError as e:
		print(e)
	else: 
		print('Instance %s terminated'%instanceID)
		#config here 
		data = {'ec2_information': {'keypair_name': str(keypairName), 'security_group_ID': str(sgID), 'ami_ID' : str(amiID)}}
		config_file = open('ec2_config.yml', 'w')
		yaml.dump(data, config_file)
		print('ec2_config file created')

# Creating AMI to run pipeline 

# Create ec2 kwy pair
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
'Description': 'Test_AMI',
'InstanceId' : instance_id,
'Name' : 'Test_AMI',
}

ami_id = create_ami(**ami_details)

clean_up(instance_id, key_name, security_group_id, ami_id)



