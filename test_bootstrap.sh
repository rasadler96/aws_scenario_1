#!/bin/sh

# This is a test bootstrap script - simply pulls in code from github

apt-get update 
apt install -y git-all

cd /home/ubuntu
git clone https://github.com/Becky-Sadler/aws-ec2-pipeline.git 
