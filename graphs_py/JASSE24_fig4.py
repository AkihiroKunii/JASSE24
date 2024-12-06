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
fname1 = "ext_flux_pi_0.6ms.dat"
fname2 = "ext_flux_pe_0.6ms.dat"
fname3 = "ext_flux_mirr_0.6ms.dat"
fname4 = "ext_flux_S_0.6ms.dat"
fname5 = "ext_flux_dotn_0.6ms.dat"
fname6 = "ext_flux_time_0.6ms.dat"
fname7 = "ext_flux_tot_mod_0.6ms.dat"
fname8 = "output_AB.dat"

#directory name
dirname ="../../../../AIP/AIP_data/2.3.v2_real_case/2.3.2.vp600s/"
dirname1="2.3.2.1.(0,0)/pypost/data/"
dirname2="2.3.2.2.(0,1)/pypost/data/"
dirname3="2.3.2.3.(1,0)/pypost/data/"
dirname4="2.3.2.4.(1,1)/pypost/data/"

#data loading
I_pi_1 = np.loadtxt(dirname+dirname1+fname1)
I_pi_2 = np.loadtxt(dirname+dirname2+fname1)
I_pi_3 = np.loadtxt(dirname+dirname3+fname1)
I_pi_4 = np.loadtxt(dirname+dirname4+fname1)
I_pe_1 = np.loadtxt(dirname+dirname1+fname2)
I_pe_2 = np.loadtxt(dirname+dirname2+fname2)
I_pe_3 = np.loadtxt(dirname+dirname3+fname2)
I_pe_4 = np.loadtxt(dirname+dirname4+fname2)
I_mirr_1 = np.loadtxt(dirname+dirname1+fname3)
I_mirr_2 = np.loadtxt(dirname+dirname2+fname3)
I_mirr_3 = np.loadtxt(dirname+dirname3+fname3)
I_mirr_4 = np.loadtxt(dirname+dirname4+fname3)
I_S_1 = np.loadtxt(dirname+dirname1+fname4)
I_S_2 = np.loadtxt(dirname+dirname2+fname4)
I_S_3 = np.loadtxt(dirname+dirname3+fname4)
I_S_4 = np.loadtxt(dirname+dirname4+fname4)
I_dotn_1 = np.loadtxt(dirname+dirname1+fname5)
I_dotn_2 = np.loadtxt(dirname+dirname2+fname5)
I_dotn_3 = np.loadtxt(dirname+dirname3+fname5)
I_dotn_4 = np.loadtxt(dirname+dirname4+fname5)
I_time_1 = np.loadtxt(dirname+dirname1+fname6)
I_time_2 = np.loadtxt(dirname+dirname2+fname6)
I_time_3 = np.loadtxt(dirname+dirname3+fname6)
I_time_4 = np.loadtxt(dirname+dirname4+fname6)
I_tot_1 = np.loadtxt(dirname+dirname1+fname7)
I_tot_2 = np.loadtxt(dirname+dirname2+fname7)
I_tot_3 = np.loadtxt(dirname+dirname3+fname7)
I_tot_4 = np.loadtxt(dirname+dirname4+fname7)

s = np.loadtxt(dirname+"/2.3.2.1.(0,0)/"+fname8, usecols=1)[2::2]

#range
rlim1 = [2, 2.8]
rlim2 = [-2.8, 2.8]
s_throat1 = 2.815
s_throatA = -2.815
width = 0.5
width_main = 1.5

I_drag_1 = I_S_1+I_dotn_1
I_drag_2 = I_S_2+I_dotn_2
I_drag_3 = I_S_3+I_dotn_3
I_drag_4 = I_S_4+I_dotn_4

I_pi_1_new= I_pi_1*1E-23
I_pe_1_new= I_pe_1*1E-23
I_mirr_1_new= I_mirr_1*1E-24
I_drag_1_new= I_drag_1*1E-23
I_time_1_new= I_time_1*1E-21
I_tot_1_new = I_tot_1*1E-22

#figure
fig = plt.figure(figsize=(5.8, 4))
plt.subplots_adjust(wspace=0.4)#, hspace=0.4)

axA = fig.add_subplot(3,2,1)
axA.plot(s,I_pi_1*1E-24, '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
axA.plot(s,I_pi_3*1E-24, '--', color="tab:green", linewidth = width_main, label="cx$_{\\rm{on}}$")
axA.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axA.text(2.1, 1.3, "(a)")
axA.set(xlim=rlim1,
        ylim=[-1.8,1.8],
        #xlabel="s [m]",
        ylabel="$\Gamma_{ \\rm{p}_{\\rm{i}}}$ \n(10$^{24}$ /m$^2$/s)",
        title="")
        
axB = fig.add_subplot(3,2,3)
axB.plot(s,I_pe_1*1E-24, '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
axB.plot(s,I_pe_3*1E-24, '--', color="tab:green", linewidth = width_main, label="cx$_{\\rm{on}}$")
axB.text(2.1, 1.3, "(b)")
axB.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axB.set(xlim=rlim1,
        ylim=[-1.8,1.8],
        #xlabel="s [m]",
        ylabel="$\Gamma_{\\rm{p}_{\\rm{e}}}$\n (10$^{24}$ /m$^2$/s)",
        title="")

axC = fig.add_subplot(3,2,5)
axC.plot(s,I_mirr_1*1E-24, '-',  linewidth = width_main, label="(cx$_{\\rm{off}}$,iz$_{\\rm{off}}$)")
axC.plot(s,I_mirr_3*1E-24, '--', color="tab:green", linewidth = width_main, label="(cx$_{\\rm{on}}$,iz$_{\\rm{off}}$)")
axC.text(2.1, 1.3, "(c)")
axC.tick_params(which='both',axis='both', direction='in', labelbottom=True)
axC.legend()
plt.legend(fontsize=10,loc='upper left', bbox_to_anchor=(0.4,-0.4), ncol=2)
axC.set(xlim=rlim1,
        ylim=[-1.8,1.8],
        xlabel="s (m)",
        ylabel="$\Gamma_{\\rm{mirror}}$\n (10$^{24}$ /m$^2$/s)",
        title="")

axD = fig.add_subplot(3,2,2)
axD.plot(s,I_drag_1*1E-22, '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
axD.plot(s,I_drag_3*1E-22, '--', color="tab:green", linewidth = width_main, label="cx$_{\\rm{on}}$")
axD.text(2.1, -2.5, "(d)")
axD.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axD.set(xlim=rlim1,
        ylim=[-5,5],
        xlabel="",
        ylabel="$\Gamma_{\\rm{drag}}$\n (10$^{22}$ /m$^2$/s)",
        title="")
        
axE = fig.add_subplot(3,2,4)
axE.plot(s,I_time_1*1E-22, '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
axE.plot(s,I_time_3*1E-22, '--', color="tab:green", linewidth = width_main, label="cx$_{\\rm{on}}$")
axE.text(2.1, -2.5, "(e)")
axE.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axE.set(xlim=rlim1,
        ylim=[-5,5],
        xlabel="",
        ylabel="$\Gamma_{\\rm{transf}}$\n (10$^{22}$ /m$^2$/s)",
        title="")

axF = fig.add_subplot(3,2,6)
axF.plot(s,I_tot_1*1E-22, '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
axF.plot(s,I_tot_3*1E-22, '--', color="tab:green", linewidth = width_main, label="cx$_{\\rm{on}}$")
axF.tick_params(which='both',axis='both', direction='in', labelbottom=True)
axF.text(2.1,-2.5, "(f)")
axF.set(xlim=rlim1,
        ylim=[-5,5],
        xlabel="s (m)",
        ylabel="$\Gamma_{\\rm{sum}}$\n (10$^{22}$ /m$^2$/s)",
        title="")


plt.savefig('../graphs/fig4.svg',bbox_inches="tight")
#plt.savefig('../2.3.2.n.img_pdf/pflux_cx(0.6 ms).pdf',bbox_inches="tight")
plt.show()
