import json
import argparse

from Preparation.Extract_Chr import Extract_Chr

parser = argparse.ArgumentParser()
parser.add_argument("JSON_file", help="path to the task-description JSON")
args = parser.parse_args()

if len(args.JSON_file) == 0:
    print("error: JSON file path is empty")
    quit()

with open(args.JSON_file) as file:
    instruction = json.load(file)

mode = instruction['mode']

# prepare selected chromosome's sequence


# manual mode
if mode == 'manual':
    pass