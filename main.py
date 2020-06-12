import sys
import os
from copy import deepcopy
from model import Sample, Aggregate
import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt
import json
import csv

def main():

	#Path to config file
	_, config_path = deepcopy(sys.argv)
	with open(config_path, "r") as f:
		config = json.load(f)
		data_dir = config["data_dir"]
		if not os.path.isdir(data_dir):
			os.mkdir(data_dir)
		specs = config["specifications"]
	run_num = 0
	for run_vals in specs:
		concentrations = run_vals["concentrations"]
		sample_volume = run_vals["sample_volume"]
		sample_diffusive_const = run_vals["sample_diffusive_const"]
		num_timesteps = run_vals["Number of timesteps"]
		num_droplets = []

		means = []
		for c in concentrations:
			print("Started run for concentration " + str(c) + "uM", flush = True)
			sample = Sample(sample_volume, sample_diffusive_const, c)
			aggs = sample.simulate(num_timesteps)
			volumes = []
			for agg in aggs:
				if agg.is_droplet():
					volumes.append(agg.volume())
			if len(volumes) == 0:
				num_droplets.append(0)
				means.append(0)
			else:
				nobs, minmax, mean, variance, skewness, kurtosis = sp.describe(volumes)
				num_droplets.append(nobs)
				means.append(mean)
			print("Finished run for concentration " + str(c) + "uM", flush = True)
		run_num += 1

		run_dir = os.path.join(data_dir, "run" + str(run_num))
		os.mkdir(run_dir)
		print(run_dir)
		with open(os.path.join(run_dir, "volumes_dat.csv"), 'w', newline = '') as csvfile:
			writer = csv.writer(csvfile, delimiter = ',')
			writer.writerow(volumes)

		with open(os.path.join(run_dir, "config_copy.txt"), 'w', newline = '') as f:
			f.write(json.dumps(config))

		plt.rc('font', family='serif')
		fig = plt.figure(figsize=(10, 10))
		ax = fig.add_subplot(1, 1, 1)
		for item in (ax.get_xticklabels() + ax.get_yticklabels()):
			item.set_fontsize(20)
		for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
			item.set_fontsize(30)
		ax.plot(concentrations, num_droplets, "sb:")
		ax.set_title("Number of Droplets")
		ax.set_xlabel("Concentration (uM)")
		ax.set_ylabel("Number of Droplets")
		plt.savefig(os.path.join(run_dir, 'num_droplets.png'), bbox_inches='tight')

		fig = plt.figure(figsize=(10, 10))
		ax = fig.add_subplot(1, 1, 1)
		for item in (ax.get_xticklabels() + ax.get_yticklabels()):
			item.set_fontsize(20)
		for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
			item.set_fontsize(30)
		ax.plot(concentrations, means, "sb:")
		ax.set_title("Mean Droplet Volume")
		ax.set_xlabel("Concentration (uM)")
		ax.set_ylabel("Mean Droplet Volume (um^3)")
		plt.savefig(os.path.join(run_dir, 'mean_droplet_volume.png'), bbox_inches='tight')




if __name__ == "__main__": 
    main()
