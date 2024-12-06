# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pandas as pd

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 10  #default 12
#plt.style.use('tableau-colorblind10')

#file name
fname1="aip.pellet.dat"
fname2="ext_pnEn.dat"
fname3="ext_gamma.dat"

#directory name
dirname ="AIP/AIP_data/2.3.v2_real_case/2.3.2.vp600s/"
dirname1="2.3.2.1.(0,0)/pypost/"
dirname2="2.3.2.2.(0,1)/pypost/"
dirname3="2.3.2.3.(1,0)/pypost/"
dirname4="2.3.2.4.(1,1)/pypost/"


#data loading
pn1, En1 = np.loadtxt("../../../../"+dirname+dirname1+fname2)
pn2, En2 = np.loadtxt("../../../../"+dirname+dirname2+fname2)
pn3, En3 = np.loadtxt("../../../../"+dirname+dirname3+fname2)
pn4, En4 = np.loadtxt("../../../../"+dirname+dirname4+fname2)
t1, G1 = np.loadtxt("../../../../"+dirname+dirname1+fname3)
t2, G2 = np.loadtxt("../../../../"+dirname+dirname2+fname3)
t3, G3 = np.loadtxt("../../../../"+dirname+dirname3+fname3)
t4, G4 = np.loadtxt("../../../../"+dirname+dirname4+fname3)

#rlim = [-13.5,13.5]
width=0.5
width_main = 1.5
pn1_new = pn1*1E-18
G1_new = G1*1E-21

ratio = 3

#figure (variables)
fig = plt.figure(figsize=(5.1, 1.4*3))

axA = fig.add_subplot(3,1,1)
axA.plot(t1*1E3,pn1*1E-18, '-', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{off}}$)")
axA.plot(t1*1E3,pn2*1E-18, ':', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{on}}$)")
axA.plot(t1*1E3,pn3*1E-18, '--',linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{off}}$)")
axA.plot(t1*1E3,pn4*1E-18, '-.',linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{on}}$)")
plt.vlines(t1[1200]*1E3, pn1_new.min()*0.9,pn1_new.max()*1.1, "black", linestyle="dashdot", linewidth=width)
plt.hlines(pn1_new[-1], 0, t1.max()*1E3, "black", linestyle="dashdot", linewidth=width)
axA.text(8, 4.5, "(a)")
axA.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axA.legend()
plt.legend(fontsize=8,loc='upper left',bbox_to_anchor=(1,1))
axA.set(xlim=[0,10],
        ylim=[3.2,5],
        xlabel="",
        ylabel="stored ion number\n (10$^{18}$)",
        title="")

axB = fig.add_subplot(3,1,3)
axB.plot(t1*1E3,En1*1E0, '-', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{off}}$)")
axB.plot(t1*1E3,En2*1E0, ':', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{on}}$)")
axB.plot(t1*1E3,En3*1E0, '--',linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{off}}$)")
axB.plot(t1*1E3,En4*1E0, '-.',linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{on}}$)")
plt.vlines(t1[1200]*1E3, 50,110, "black", linestyle="dashdot", linewidth=width)
plt.hlines(En4[-1], 0, t1.max()*1E3, "black", linestyle="dashdot", linewidth=width)
axB.text(8, 90, "(c)")
#axB.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
axB.tick_params(which='both',axis='both', direction='in', labelbottom=True)
# axB.legend()
# plt.legend(fontsize=10,loc='lower right')
axB.set(xlim=[0,10],
        ylim=[50,110],
        xlabel="time (ms)",
        ylabel="stored energy\n (J)",
        title="")

axC = fig.add_subplot(3,1,2)
axC.plot(t1*1E3,G1*1E-21, '-', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{off}}$)")
axC.plot(t1*1E3,G2*1E-21, ':', linewidth = width_main, label="(cx$_{\\rm{off}}$, iz$_{\\rm{on}}$)")
axC.plot(t1*1E3,G3*1E-21, '--',linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{off}}$)")
axC.plot(t1*1E3,G4*1E-21, '-.',linewidth = width_main, label="(cx$_{\\rm{on}}$, iz$_{\\rm{on}}$)")
plt.vlines(t1[1200]*1E3, 0,100, "black", linestyle="dashdot", linewidth=width)
plt.hlines(G4[-1]*1E-21, 0, t1.max()*1E3, "black", linestyle="dashdot", linewidth=width)
axC.text(8, 2, "(b)")
#axB.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
axC.tick_params(which='both',axis='both', direction='in', labelbottom=False)
# axC.legend()
# plt.legend(fontsize=10,loc='lower right')
axC.set(xlim=[0,10],
        ylim=[0.5,2.5],
        xlabel="",
        ylabel="$\Gamma$A (10$^{21}$ /s)",
        title="")

plt.savefig('../graphs/fig3.svg',bbox_inches="tight")
#plt.savefig('../2.3.2.n.img_pdf/Eg_cx,iz.pdf',bbox_inches="tight")
fig.align_labels()
plt.show()
