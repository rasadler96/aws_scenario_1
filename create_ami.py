import boto3
import yaml

config = yaml.safe_load(open("config.yml"))
 
access_key = config['aws_access_key_id']
secret_key = config['aws_secret_access_key']
default_region = config['aws_default_region']