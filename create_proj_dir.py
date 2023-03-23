#!/usr/bin/env python3

# copy everything from the template directory to the new project directory
# accept the project numer as a command line argument
# use argparse to parse the command line arguments

import argparse
import os
import shutil

parser = argparse.ArgumentParser(description='Create a new project directory')
parser.add_argument('-n', '--number', help='Project number', required=True)
args = parser.parse_args()

# if args.number is not a number, exit
if not args.number.isdigit():
    print('Project number must be a number')
    exit()

# if args.number is not 2 digits, zero pad it
if len(args.number) == 1:
    args.number = '0' + args.number

# get the location of this script
root = os.path.dirname(os.path.abspath(__file__))
template_dir = 'template'

# create the new project directory
proj_dir = 'proj' + args.number
path_to_proj = os.path.join(root, 'projects', proj_dir)
if os.path.exists(path_to_proj):
    print('Project directory already exists. Force overwrite? (y/n)', end=' ')
    if input() == 'y':
        shutil.rmtree(path_to_proj)
    else:
        exit()

# copy the template directory to the new project directory
shutil.copytree(os.path.join(root, template_dir),
                os.path.join(path_to_proj))

# rename projxx.Rproj to proj{args.number}.Rproj
os.rename(os.path.join(path_to_proj, 'projxx.Rproj'),
          os.path.join(path_to_proj, 'proj' + args.number + '.Rproj'))

print('Project directory created.')
