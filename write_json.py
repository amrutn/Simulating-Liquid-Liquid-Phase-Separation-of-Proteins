"""
All units are in um and uM
"""

import json
import sys
from copy import deepcopy
import os

arguments = sys.argv
if len(sys.argv) == 1:
    config_path = os.path.relpath("config.json")
else:
    _, config_path = deepcopy(sys.argv)
directory = os.path.split(config_path)[0]
if not os.path.isdir(directory):
	os.mkdir(directory)


config = {"data_dir" : "../sim_data",
"specifications" : [
	{
	"concentrations" : [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 10],
	"sample_volume" : 100,
	"sample_diffusive_const" : 5 * 0.01,
	"Number of timesteps" : 100,
	"Molecular Radius" : 0.01,
	"Min Droplet Volume" : 0.00015
	}, 
	{
	"concentrations" : [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 10, 12.8, 25.6],
	"sample_volume" : 20,
	"sample_diffusive_const" : 5 * 0.01, 
	"Number of timesteps" : 30,
	"Molecular Radius" : 0.01,
	"Min Droplet Volume" : 0.00015
	}
	]
	}



json_string = json.dumps(config)


with open(config_path, "w") as config_json:
	config_json.write(json_string)