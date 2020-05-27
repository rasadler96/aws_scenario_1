import boto3 
import yaml
import botocore

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

def get_DNS(instance_id):
	try:
		response = ec2_client.describe_instances(InstanceIds=[instance_id])
	except botocore.exceptions.ClientError as e: 
		print(e)
	else:
		public_DNS = response['Reservations'][0]['Instances'][0]['PublicDnsName']
		return public_DNS

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
'InstanceType' : 't2.micro',
'KeyName' : keypair_name,
'MinCount' : 1,
'MaxCount' : 1,
'SecurityGroupIds' : [
    security_group_id,
]}


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
	public_DNS = get_DNS(instance_id)
	


else:
	print('There is a problem with the selected AMI - state is "' + str(state) + '"')

