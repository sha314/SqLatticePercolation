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
	# if L > 200:
	# 	continue
	x, y, clr = np.loadtxt(file, unpack=True, delimiter=',')

	fig, ax = plt.subplots(figsize=(10, 10), frameon=False, dpi=300)
	

	# ax = fig.add_axes([0, 0, L, L])
	ax.set_position([0.02, 0.02, 0.98, 0.98])
	fig.patch.set_visible(False)
	ax.axis('off')

	index = clr == 0
	ax.scatter(x[index], y[index], c='k', marker="s", s=3, edgecolors='none')

	index = clr > 0
	ax.scatter(x[index], y[index], c=clr[index], marker="s", s=3,  alpha=0.7)


	# plt.legend()
	# plt.show()
	# plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='on')
	fig.savefig(file.split('.')[0] + "figure")
