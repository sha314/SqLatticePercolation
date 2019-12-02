import numpy as np
import matplotlib.pyplot as plt
import json
import glob

signature="SitePercolation_ps_v10_lattice_"
# filename="SitePercolation_ps_v10_lattice_L100_2019-12-02_192131"
# filename="SitePercolation_ps_v10_lattice_L200_2019-12-02_192700"
# filename="SitePercolation_ps_v10_lattice_L400_2019-12-02_193551"
# filename += ".txt"

files = glob.glob(signature + "*txt")

for file in files:
	with open(file) as f:
		line = f.readline()
		head = json.loads(line[1:])
		L = int(head['length'])
		pass
	# if L > 5:
		# continue
	x, y, clr = np.loadtxt(file, unpack=True, delimiter=',')

	fig = plt.figure(figsize=(10, 10))
	
	index = clr == 0
	plt.scatter(x[index], y[index], c='k', marker="s", s=2, label="Spanning", edgecolors='none')

	index = clr > 0
	plt.scatter(x[index], y[index], c=clr[index], marker="s", s=2,  alpha=0.4)


	plt.legend()
	# plt.show()
	plt.savefig(file.split('.')[0] + "figure.eps")
