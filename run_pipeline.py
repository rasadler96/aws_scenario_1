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
		print(state)

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


# Running script from here (Above are the functions)

image_state(ami_id)

if state == 'available':
	print('AMI is available')

else:
	print('There is a problem with the selected AMI - state is "' + state + '"')

