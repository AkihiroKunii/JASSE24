# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pandas as pd

"""
displays temprature distribution
"""
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 10  #default 12

#file name
fname1 = "ext_I_ratio.dat"
fname2 = "aip.time.dat"

#directory name
dirname ="../../../../AIP/AIP_data/2.3.v2_real_case/2.3.2.vp600s/"
dirname1="2.3.2.1.(0,0)/pypost/data/"
dirname2="2.3.2.2.(0,1)/pypost/data/"
dirname3="2.3.2.3.(1,0)/pypost/data/"
dirname4="2.3.2.4.(1,1)/pypost/data/"

#data loading
ratio_1 = np.loadtxt(dirname+dirname1+fname1)
ratio_2 = np.loadtxt(dirname+dirname2+fname1)
ratio_3 = np.loadtxt(dirname+dirname3+fname1)
ratio_4 = np.loadtxt(dirname+dirname4+fname1)
t       = np.loadtxt(dirname+"/2.3.2.4.(1,1)/"+fname2, usecols=1, skiprows=1)

#range
rlim1 = [0, 13.5]
rlim2 = [-2.8, 2.8]
s_throat1 = 2.815
s_throatA = -2.815
width = 0.5
width_main = 1.5

#figure
fig = plt.figure(figsize=(5.2, 2))

axB = fig.add_subplot(1,1,1) 
axB.plot(t*1E3,ratio_1, '-',  linewidth = width_main, label="(cx$_{\\rm{off}}$,iz$_{\\rm{off}}$)")
axB.plot(t*1E3,ratio_3, '--', color="tab:green", linewidth = width_main, label="(cx$_{\\rm{on}}$,iz$_{\\rm{off}}$)")
axB.hlines(ratio_4[-1], 0, 21, "black", linestyle="dashdot", linewidth=width)
axB.vlines(0.6, 0, 2, "black", linestyle="dashdot", linewidth=width)
axB.tick_params(which='both',axis='both', direction='in', labelbottom=True)
axB.legend()
plt.legend(fontsize=10,loc='upper right')
axB.set(xlim=[0,10],
        ylim=[1,1.06],
        xlabel="time (ms)",
        ylabel="|$\Gamma_{\\rm{p}_{\\rm{i}}}$+$\Gamma_{\\rm{p}_{\\rm{e}}}$|/|$\Gamma_{\\rm{mirr}}$|",
        title="")

plt.savefig('../graphs/fig5.svg',bbox_inches="tight")
#plt.savefig('../2.3.2.n.img_pdf/flux_ratio_cx,iz.pdf',bbox_inches="tight")
plt.show()
