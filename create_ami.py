import boto3
import yaml
import os

config = yaml.safe_load(open("config.yml"))
 
access_key = config['aws_access_key_id']
secret_key = config['aws_secret_access_key']
default_region = config['aws_default_region']

# Creating session (configures credentials and default region)
session = boto3.Session(
	aws_access_key_id = access_key,
	aws_secret_access_key = secret_key,
	region_name = default_region
)

# Function to create key pair to be used to launch instance - input = key pair name (minus the .pem), output = keypair.pem with correct permissions. 
def create_keypair(name_of_keypair):
	ec2 = session.resource('ec2', region_name=default_region)
	key_file = open('%s.pem'%name_of_keypair,'w')
	key_pair = ec2.create_key_pair(KeyName=name_of_keypair, )
	key_pair_contents = str(key_pair.key_material)
	key_file.write(key_pair_contents)
	os.system('chmod 400 %s.pem'%name_of_keypair)







