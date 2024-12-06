# coding: utf-8
import re
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
fname1 = "aip.var.001200.dat"
fname2 = "aip.var.020000.dat"
fname3 = "output_AB.dat"
fname4 = "aip.pellet.dat"

#directory name
dirname ="AIP/AIP_data/2.3.v2_real_case/2.3.2.vp600s/"
dirname1="2.3.2.1.(0,0)/"
dirname2="2.3.2.2.(0,1)/"
dirname3="2.3.2.3.(1,0)/"
dirname4="2.3.2.4.(1,1)/"

#data loading
data_1 = np.loadtxt("../../../../"+dirname+dirname1+fname1, dtype='float',encoding="utf-8_sig",\
                  converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)
data_2 = np.loadtxt("../../../../"+dirname+dirname2+fname1, dtype='float',encoding="utf-8_sig",\
                  converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)
data_3 = np.loadtxt("../../../../"+dirname+dirname3+fname1, dtype='float',encoding="utf-8_sig",\
                  converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)
data_4 = np.loadtxt("../../../../"+dirname+dirname4+fname1, dtype='float',encoding="utf-8_sig",\
                  converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)
data_std = np.loadtxt("../../../../"+dirname+dirname4+fname2, dtype='float',encoding="utf-8_sig",\
                  converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)

Ti_para_1       = data_1[:,2]
Ti_para_2       = data_2[:,2]
Ti_para_3       = data_3[:,2]
Ti_para_4       = data_4[:,2]
Ti_para_std     = data_std[:,2]
Ti_perp_1       = data_1[:,3]
Ti_perp_2       = data_2[:,3]
Ti_perp_3       = data_3[:,3]
Ti_perp_4       = data_4[:,3]
Ti_perp_std     = data_std[:,3]
s1 = np.loadtxt("../../../../"+dirname+dirname1+fname3, usecols=1)[2::2]
B = np.loadtxt("../../../../"+dirname+dirname1+fname3, usecols=4)[2::2]

dt = 1E-8
diag_period = 50

#range
rlim1 = [0, 2.8]
rlim2 = [-2.8, 2.8]
s_throat1 = 2.815
s_throatA = -2.815
width = 0.5
width_main = 1.5

ratio = 6

aniso_Ti_1, aniso_Ti_2, aniso_Ti_3, aniso_Ti_4, aniso_Ti_std = Ti_perp_1/Ti_para_1, Ti_perp_2/Ti_para_2, Ti_perp_3/Ti_para_3, Ti_perp_4/Ti_para_4, Ti_perp_std/Ti_para_std

#figure (variables)
fig = plt.figure(figsize=(4.8, 2))
plt.subplots_adjust(wspace=0.3)

axB = fig.add_subplot(1,1,1) 
axB.plot(s1,aniso_Ti_1, '-',  linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{off}}$)")
# axB.plot(s1,n2*10E-19, ':', linewidth = width_main , label="(cx$_{\\rm{off}}$, iz$_{\\rm{on}}$)")
axB.plot(s1,aniso_Ti_3, '--', linewidth = width_main , color="tab:green", label="(cx$_{\\rm{on}}$, iz$_{\\rm{off}}$)")
# axB.plot(s1,n4*10E-19, '-.', linewidth = width_main , label="(cx$_{\\rm{on}}$, iz$_{\\rm{on}}$)")
axB.plot(s1,aniso_Ti_std, '-', linewidth = width_main , color="black", label="steady")
axB.tick_params(which='both',axis='both', direction='in', labelbottom=True)
axB.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axB.legend()
plt.legend(fontsize=10,loc='upper left', bbox_to_anchor=(1,1))
axB.set(xlim=rlim1,
        ylim=[1.5,3],
        xlabel="s (m)",
        ylabel="$T_{\\rm{i}\perp} / T_{\\rm{i}\parallel}$",
        title="")

plt.savefig('../graphs/fig6.svg',bbox_inches="tight")
#plt.savefig('../2.3.2.n.img_pdf/JASSE24_fig2.pdf',bbox_inches="tight")
plt.show()
