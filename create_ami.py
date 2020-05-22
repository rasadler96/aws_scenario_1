import boto3
import yaml

config = yaml.safe_load(open("config.yml"))
 
access_key = config['aws_access_key_id']
secret_key = config['aws_secret_access_key']
default_region = config['aws_default_region']

session = boto3.Session(
	aws_access_key_id=access_key
	aws_secret_access_key=secret_key
	aws_default_region=default_region
)