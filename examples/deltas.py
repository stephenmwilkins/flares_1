
import os
import sys
import numpy as np
from scipy import stats

import matplotlib as mpl
import matplotlib.cm as cm
import h5py

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import flares_1.analysis.plt as fplt

deltas = np.array([0.969639,0.918132,0.851838,0.849271,0.845644,0.842128,0.841291,0.83945,0.838891,0.832753,0.830465,0.829349,0.827842,0.824159,0.821425,0.820476,0.616236,0.616012,0.430745,0.430689,0.266515,0.266571,0.121315,0.121147,-0.007368,-0.007424,-0.121207,-0.121319,-0.222044,-0.222156,-0.311441,-0.311329,-0.066017,-0.066185,-0.00748,-0.007424,0.055076,0.054909,-0.47874,-0.433818])




xlims = [28, 30.5]
ylims = [-2.75,-1.0]

fig, ax = fplt.single()


mu = -0.0032
sigma = 0.0529
x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
ax.fill_between(x, x*0.0, stats.norm.pdf(x, mu, sigma), color='0.9', zorder = 1)


ax.hist(np.log10(deltas+1),bins=28, range=[-0.35, 0.35], ls='-',color='k', alpha=1.0, lw=1, histtype='step', zorder = 2)

# ax.set_ylim(ylims)
# ax.set_xlim(xlims)
# ax.set_ylabel()
ax.set_xlabel(r'$\rm\log_{10}(1+\delta)$')
ax.grid(True)
fig.savefig('figs/delta.pdf')
