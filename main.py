from model import Sample, Aggregate
import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt

def main():
	concentrations = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8]
	num_droplets = []
	means = []
	for c in concentrations:
		print("Started run for concentration " + str(c), flush = True)
		sample = Sample(1, 5 * 0.01, c)
		aggs = sample.simulate(10)
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
		print("Finished run for concentration " + str(c), flush = True)


	plt.plot(concentrations, num_droplets)
	plt.title("Number of Droplets")
	plt.xlabel("Concentration")
	plt.ylabel("Number of Droplets")
	plt.show()
	plt.close()

	plt.plot(concentrations, means)
	plt.title("Mean Droplet Area")
	plt.xlabel("Concentration")
	plt.ylabel("Mean Droplet Area")
	plt.show()
	plt.close()


if __name__ == "__main__": 
    main()
