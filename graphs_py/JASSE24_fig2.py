# coding: utf-8
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MaxNLocator

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


n1      = data_1[:,0]
n2      = data_2[:,0]
n3      = data_3[:,0]
n4      = data_4[:,0]
n_std   = data_std[:,0]
V1      = data_1[:,1]
V2      = data_2[:,1]
V3      = data_3[:,1]
V4      = data_4[:,1]
V_std   = data_std[:,1]
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
Te_1    = data_1[:,4]
Te_2    = data_2[:,4]
Te_3    = data_3[:,4]
Te_4    = data_4[:,4]
Te_std  = data_std[:,4]
s1 = np.loadtxt("../../../../"+dirname+dirname1+fname3, usecols=1)[2::2]
B = np.loadtxt("../../../../"+dirname+dirname1+fname3, usecols=4)[2::2]

dt = 1E-8
diag_period = 50

n1_new = n1*10E-19

#range
rlim1 = [0, 13.5]
rlim2 = [-2.8, 2.8]
s_throat1 = 2.815
s_throatA = -2.815
width = 0.5
width_main = 1.5

labelpad_value = 15  # y軸ラベルと軸の距離
ylabel_position = -0.15  # y軸ラベルの位置（x座標で調整）

ratio = 6

#figure (variables)
fig = plt.figure(figsize=(5.3, 3.5))
plt.subplots_adjust(wspace=0.3)

axA = fig.add_subplot(3,2,1) 
axA.plot(s1,B, '-', linewidth = width_main , color="black", label="steady")
axA.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axA.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axA.text(11,3,"(a)")
axA.set(xlim=rlim1,
        ylim=[0,4],
        xlabel="",
        title="")
axA.set_ylabel("B (T)", labelpad=labelpad_value)
axA.yaxis.set_label_coords(ylabel_position, 0.5)

axB = fig.add_subplot(3,2,3) 
axB.plot(s1,n1*10E-19, '-',  linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{off}}$)")
axB.plot(s1,n2*10E-19, ':', linewidth = width_main , label="(cx$_{\\rm{off}}$, iz$_{\\rm{on}}$)")
axB.plot(s1,n3*10E-19, '--', linewidth = width_main , label="(cx$_{\\rm{on}}$, iz$_{\\rm{off}}$)")
axB.plot(s1,n4*10E-19, '-.', linewidth = width_main , label="(cx$_{\\rm{on}}$, iz$_{\\rm{on}}$)")
axB.plot(s1,n_std*10E-19, '-', linewidth = width_main , color="black", label="steady")
axB.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axB.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axB.text(11,3,"(b)")
# axB.yaxis.set_major_locator(MaxNLocator(integer=True))
axB.set(xlim=rlim1,
        ylim=[0,5],
        xlabel="",
        title="")
axB.set_ylabel("n (10$^{18}$ /m$^3$)", labelpad=labelpad_value)
axB.yaxis.set_label_coords(ylabel_position, 0.5)

axC = fig.add_subplot(3,2,5) 
axC.plot(s1,V1*1E-5, '-', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{off}}$)")
axC.plot(s1,V2*1E-5, ':', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{on}}$)")
axC.plot(s1,V3*1E-5, '--', linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{off}}$)")
axC.plot(s1,V4*1E-5, '-.', linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{on}}$)")
axC.plot(s1,V_std*1E-5, '-', linewidth = width_main, color="black", label="steady")
axC.tick_params(which='both',axis='both', direction='in', labelbottom=True)
axC.vlines(2.8, 0, 100, "black", linestyle="dashdot", linewidth=width)
axC.text(11,0.5,"(c)")
axC.set(xlim=rlim1,
        ylim=[0,3],
        xlabel="s (m)",
        ylabel="",
        title="")
axC.set_ylabel("V (10$^{5}$m/s)", labelpad=labelpad_value)
axC.yaxis.set_label_coords(ylabel_position, 0.5)

axD = fig.add_subplot(3,2,2) 
axD.plot(s1,Ti_perp_1, '-', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{off}}$)")
axD.plot(s1,Ti_perp_2, ':', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{on}}$)")
axD.plot(s1,Ti_perp_3, '--', linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{off}}$)")
axD.plot(s1,Ti_perp_4, '-.', linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{on}}$)")
axD.plot(s1,Ti_perp_std, '-', linewidth = width_main, color="black", label="steady")
axD.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axD.vlines(2.8, 0, 150, "black", linestyle="dashdot", linewidth=width)
axD.text(11,100,"(d)")
axD.set(xlim=rlim1,
        ylim=[0,Ti_perp_std.max()*1.1],
        xlabel="",
        ylabel="",
        title="")
axD.set_ylabel("T$_{\\rm{i}\\perp}$ (eV)", labelpad=labelpad_value)
axD.yaxis.set_label_coords(ylabel_position, 0.5)
        
axE = fig.add_subplot(3,2,4) 
axE.plot(s1,Ti_para_1, '-', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{off}}$)")
axE.plot(s1,Ti_para_2, ':', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{on}}$)")
axE.plot(s1,Ti_para_3, '--', linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{off}}$)")
axE.plot(s1,Ti_para_4, '-.', linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{on}}$)")
axE.plot(s1,Ti_para_std, '-', linewidth = width_main, color="black", label="steady")
axE.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axE.vlines(2.8, 0, 100, "black", linestyle="dashdot", linewidth=width)
axE.text(11,60,"(e)")
axE.set(xlim=rlim1,
        ylim=[0,Ti_para_std.max()*1.1],
        xlabel="",
        ylabel="",
        title="")
axE.set_ylabel("T$_{\\rm{i}\\parallel}$ (eV)", labelpad=labelpad_value)
axE.yaxis.set_label_coords(ylabel_position, 0.5)

axF = fig.add_subplot(3,2,6)
axF.plot(s1,Te_1, '-', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{off}}$)")
axF.plot(s1,Te_2, ':', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{on}}$)")
axF.plot(s1,Te_3, '--', linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{off}}$)")
axF.plot(s1,Te_4, '-.', linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{on}}$)")
axF.plot(s1,Te_std, '-', linewidth = width_main, color="black", label="steady")
axF.vlines(2.8, 0, 100, "black", linestyle="dashdot", linewidth=width)
axF.tick_params(which='both',axis='both', direction='in')
axF.text(11,22,"(f)")
#plt.yscale("log")
axF.legend()
plt.legend(fontsize=10,ncol=2,loc='upper center',bbox_to_anchor=(-0.2,-0.4))
axF.set(xlim=rlim1,
        ylim=[0,Te_std.max()*1.1],
        xlabel="s (m)",
        ylabel="",
        title="")
axF.set_ylabel("T$_{\\rm{e}}$ (eV)", labelpad=labelpad_value)
axF.yaxis.set_label_coords(ylabel_position, 0.5)


plt.savefig('../graphs/fig2.svg',bbox_inches="tight")
#plt.savefig('../2.3.2.n.img_pdf/JASSE24_fig2.pdf',bbox_inches="tight")
plt.show()
