# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pandas as pd

"""
displays temprature distribution
"""

plt.rcParams["font.size"] = 10  #default 12
plt.rcParams["font.family"] = "Helvetica"

#file name
dirname ="AIP/AIP_data/2.3.v2_real_case/2.3.2.vp600s/"
fname1 = "sigvCX,ION(Tn1.0E-1).dat"

#data loading
T  = np.loadtxt("../../../../"+dirname+fname1, usecols=0, skiprows=1)
sigv_cx = np.loadtxt("../../../../"+dirname+fname1, usecols=1, skiprows=1)
sigv_iz = np.loadtxt("../../../../"+dirname+fname1, usecols=2, skiprows=1)

#figure
fig = plt.figure(figsize=(2.8, 2.8))

axA = fig.add_subplot(1,1,1)
axA.plot(T,sigv_cx, '-', label="charge exchange")
axA.plot(T,sigv_iz, '--', label="ionization/excitation")
axA.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axA.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
axA.tick_params(which='both',axis='both', direction='in', labelbottom=True)
plt.yscale('log')
plt.xscale('log')
axA.legend()
plt.legend(fontsize=10,loc='upper left',bbox_to_anchor=(0,1))
axA.set(xlim=[10,120],
        ylim=[1E-14,1E-13],
        xlabel="T (eV)",
        ylabel="$\langle$ $\sigma v$ $\\rangle$ (m$^{3}$/s)",
        title="")

plt.savefig('../graphs/fig1.svg',bbox_inches="tight")
#plt.savefig('../2.3.2.n.img_pdf/sigv_cx,iz.pdf',bbox_inches="tight")
plt.show()
