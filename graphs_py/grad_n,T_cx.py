# coding: utf-8
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pandas as pd
from scipy.ndimage import gaussian_filter1d

# display temperature distribution
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 10  # default 12

# File names and directory names
fname1 = "aip.var.001200.dat"
fname2 = "aip.var.020000.dat"
fname4 = "ext_gradient.dat"
fname3 = "output_AB.dat"
dirname ="../../../../AIP/AIP_data/2.3.v2_real_case/2.3.2.vp600s/"
dirname1 = "2.3.2.1.(0,0)/"
dirname2 = "2.3.2.2.(0,1)/"
dirname3 = "2.3.2.3.(1,0)/"
dirname4 = "2.3.2.4.(1,1)/"

# Data loading with error handling for scientific notation
data_1 = np.loadtxt(dirname + dirname1 + fname1, dtype='float', encoding="utf-8_sig", 
                    converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)
data_2 = np.loadtxt(dirname + dirname2 + fname1, dtype='float', encoding="utf-8_sig", 
                    converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)
data_3 = np.loadtxt(dirname + dirname3 + fname1, dtype='float', encoding="utf-8_sig", 
                    converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)
data_4 = np.loadtxt(dirname + dirname4 + fname1, dtype='float', encoding="utf-8_sig", 
                    converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)
data_std = np.loadtxt(dirname + dirname4 + fname2, dtype='float', encoding="utf-8_sig", 
                      converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i), skiprows=1)
data_grad_1 = np.loadtxt(dirname + dirname1 + "pypost/" + fname4, dtype='float', encoding="utf-8_sig", 
                      converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i))
data_grad_3 = np.loadtxt(dirname + dirname3 + "pypost/" + fname4, dtype='float', encoding="utf-8_sig", 
                      converters=lambda i: re.sub(r'(?<=\d)([+-])', r'E\1', i))

# Extracting columns
n1, n2, n3, n4, n_std = data_1[:, 0], data_2[:, 0], data_3[:, 0], data_4[:, 0], data_std[:, 0]
Ti_para_1, Ti_para_2, Ti_para_3, Ti_para_4, Ti_para_std = data_1[:, 2], data_2[:, 2], data_3[:, 2], data_4[:, 2], data_std[:, 2]
Ti_perp_1, Ti_perp_2, Ti_perp_3, Ti_perp_4, Ti_perp_std = data_1[:, 3], data_2[:, 3], data_3[:, 3], data_4[:, 3], data_std[:, 3]
Te_1, Te_2, Te_3, Te_4, Te_std = data_1[:, 4], data_2[:, 4], data_3[:, 4], data_4[:, 4], data_std[:, 4]
grad_n_abl_1, grad_Ti_para_abl_1, grad_Te_abl_1, grad_n_std_1, grad_Ti_para_std_1, grad_Te_std_1, ratio_grad_n1, ratio_grad_Ti_para_1, ratio_grad_Te_1 = data_grad_1[:,0], data_grad_1[:,1], data_grad_1[:,2], data_grad_1[:,3], data_grad_1[:,4], data_grad_1[:,5], data_grad_1[:,6], data_grad_1[:,7], data_grad_1[:,8] 
grad_n_abl_3, grad_Ti_para_abl_3, grad_Te_abl_3, grad_n_std_3, grad_Ti_para_std_3, grad_Te_std_3, ratio_grad_n3, ratio_grad_Ti_para_3, ratio_grad_Te_3 = data_grad_3[:,0], data_grad_3[:,1], data_grad_3[:,2], data_grad_3[:,3], data_grad_3[:,4], data_grad_3[:,5], data_grad_3[:,6], data_grad_3[:,7], data_grad_3[:,8]  

s1 = np.loadtxt(dirname + dirname1 + fname3, usecols=1)[2::2]
B = np.loadtxt(dirname + dirname1 + fname3, usecols=4)[2::2]

# Initialize gradients and ratios with the correct shape
# grad_n1 = np.zeros_like(n_std)
# grad_n2 = np.zeros_like(n_std)
# grad_n3 = np.zeros_like(n_std)
# grad_n4 = np.zeros_like(n_std)
# grad_n_std = np.zeros_like(n_std)
# grad_Ti_para_1 = np.zeros_like(n_std)
# grad_Ti_para_2 = np.zeros_like(n_std)
# grad_Ti_para_3 = np.zeros_like(n_std)
# grad_Ti_para_4 = np.zeros_like(n_std)
# grad_Ti_para_std = np.zeros_like(n_std)
# grad_Ti_perp_1 = np.zeros_like(n_std)
# grad_Ti_perp_2 = np.zeros_like(n_std)
# grad_Ti_perp_3 = np.zeros_like(n_std)
# grad_Ti_perp_4 = np.zeros_like(n_std)
# grad_Ti_perp_std = np.zeros_like(n_std)
# grad_Te_1 = np.zeros_like(n_std)
# grad_Te_2 = np.zeros_like(n_std)
# grad_Te_3 = np.zeros_like(n_std)
# grad_Te_4 = np.zeros_like(n_std)
# grad_Te_std = np.zeros_like(n_std)
diff_T_1 = np.zeros_like(n_std)
diff_T_2 = np.zeros_like(n_std)
diff_T_3 = np.zeros_like(n_std)
diff_T_4 = np.zeros_like(n_std)
diff_T_std = np.zeros_like(n_std)

#smoothing
sm_Ti_para_1, sm_Ti_para_2, sm_Ti_para_3, sm_Ti_para_4, sm_Ti_para_std = gaussian_filter1d(Ti_para_1, sigma=1), gaussian_filter1d(Ti_para_2, sigma=1), gaussian_filter1d(Ti_para_3, sigma=1), gaussian_filter1d(Ti_para_4, sigma=1), gaussian_filter1d(Ti_para_std, sigma=1) 
sm_Te_1, sm_Te_2, sm_Te_3, sm_Te_4, sm_Te_std = gaussian_filter1d(Te_1, sigma=10), gaussian_filter1d(Te_2, sigma=10), gaussian_filter1d(Te_3, sigma=10), gaussian_filter1d(Te_4, sigma=10), gaussian_filter1d(Te_std, sigma=10) 

nmax = len(n_std)

# Define gradient calculation function
ds = 1E-2
def gradient(param):
    grad = np.zeros_like(param)
    for i in range( len(param)):
        if i == nmax-1:
            grad[i] = (param[0] - param[nmax-2]) / (2*ds)
        elif i == 0:
            grad[i] = (param[1] - param[nmax-1]) / (2*ds)
        else:
            grad[i] = (param[i+1] - param[i-1]) / (2*ds)
        # if grad[i] <= 1E-5 and grad[i] >= -1E-5:
        #     grad[i] = 0
    return grad

# Calculate gradients for each parameter
# grad_n1, grad_n2, grad_n3, grad_n4, grad_n_std = gradient(n1), gradient(n2), gradient(n3), gradient(n4), gradient(n_std)
# grad_Ti_para_1, grad_Ti_para_2, grad_Ti_para_3, grad_Ti_para_4, grad_Ti_para_std = gradient(Ti_para_1), gradient(Ti_para_2), gradient(Ti_para_3), gradient(Ti_para_4), gradient(Ti_para_std)
# grad_Ti_perp_1, grad_Ti_perp_2, grad_Ti_perp_3, grad_Ti_perp_4, grad_Ti_perp_std = gradient(Ti_perp_1), gradient(Ti_perp_2), gradient(Ti_perp_3), gradient(Ti_perp_4), gradient(Ti_perp_std)
# grad_Te_1, grad_Te_2, grad_Te_3, grad_Te_4, grad_Te_std = gradient(Te_1), gradient(Te_2), gradient(Te_3), gradient(Te_4), gradient(Te_std)
diff_T_1, diff_T_2, diff_T_3, diff_T_4, diff_T_std = Ti_para_1 - Ti_perp_1, Ti_para_2 - Ti_perp_2, Ti_para_3 - Ti_perp_3, Ti_para_4 - Ti_perp_4, Ti_para_std - Ti_perp_std

# grad_Ti_para_1, grad_Ti_para_2, grad_Ti_para_3, grad_Ti_para_4, grad_Ti_para_std = gradient(sm_Ti_para_1), gradient(sm_Ti_para_2), gradient(sm_Ti_para_3), gradient(sm_Ti_para_4), gradient(sm_Ti_para_std)
# grad_Te_1, grad_Te_2, grad_Te_3, grad_Te_4, grad_Te_std = gradient(sm_Te_1), gradient(sm_Te_2), gradient(sm_Te_3), gradient(sm_Te_4), gradient(sm_Te_std)

# Calculate ratios
ratio_n1, ratio_n2, ratio_n3, ratio_n4 = n1 / n_std, n2 / n_std, n3 / n_std, n4 / n_std
ratio_Ti_para_1, ratio_Ti_para_2, ratio_Ti_para_3, ratio_Ti_para_4 = Ti_para_1 / Ti_para_std, Ti_para_2 / Ti_para_std, Ti_para_3 / Ti_para_std, Ti_para_4 / Ti_para_std
ratio_Ti_perp_1, ratio_Ti_perp_2, ratio_Ti_perp_3, ratio_Ti_perp_4 = Ti_perp_1 / Ti_perp_std, Ti_perp_2 / Ti_perp_std, Ti_perp_3 / Ti_perp_std, Ti_perp_4 / Ti_perp_std
ratio_Te_1, ratio_Te_2, ratio_Te_3, ratio_Te_4 = Te_1 / Te_std, Te_2 / Te_std, Te_3 / Te_std, Te_4 / Te_std
#ratio_grad_n1, ratio_grad_n2, ratio_grad_n3, ratio_grad_n4 = grad_n1 / grad_n_std, grad_n2 / grad_n_std, grad_n3 / grad_n_std, grad_n4 / grad_n_std
#ratio_grad_Ti_para_1, ratio_grad_Ti_para_2, ratio_grad_Ti_para_3, ratio_grad_Ti_para_4 = grad_Ti_para_1 / grad_Ti_para_std, grad_Ti_para_2 / grad_Ti_para_std, grad_Ti_para_3 / grad_Ti_para_std, grad_Ti_para_4 / grad_Ti_para_std
#ratio_grad_Te_1, ratio_grad_Te_2, ratio_grad_Te_3, ratio_grad_Te_4 = grad_Te_1 / grad_Te_std, grad_Te_2 / grad_Te_std, grad_Te_3 / grad_Te_std, grad_Te_4 / grad_Te_std
ratio_delta_T_1, ratio_delta_T_2, ratio_delta_T_3, ratio_delta_T_4 = diff_T_1 / diff_T_std, diff_T_2 / diff_T_std, diff_T_3 / diff_T_std, diff_T_4 / diff_T_std

#range
rlim1 = [2.0, 2.8]
rlim2 = [-2.8, 2.8]
s_throat1 = 2.815
s_throatA = -2.815
width = 0.5
width_main = 1.5

#figure
fig = plt.figure(figsize=(5.8, 4.8))
plt.subplots_adjust(wspace=0.43)

axA = fig.add_subplot(4,2,1) 
axA.plot(s1[:3000], ratio_n1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axA.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axA.plot(s1[:3000], ratio_n3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
#axA.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
axA.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axA.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axA.text(2.2,2.1,"(a)")
axA.set(xlim=rlim1,
        ylim=[1,3],
        xlabel="",
        ylabel="$\\frac{n}{n_{\\rm{std}}} $",
        title="")

axB = fig.add_subplot(4,2,3) 
axB.plot(s1[:3000], ratio_Ti_para_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axB.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axB.plot(s1[:3000], ratio_Ti_para_3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
#axB.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
axB.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axB.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axB.text(2.2,0.25,"(b)")
axB.set(xlim=rlim1,
        ylim=[0,1],
        xlabel="",
        ylabel="$\\frac{T_{\\rm{i}\\parallel}}{T_{\\rm{i}\\parallel,\\rm{std}}} $",
        title="")
        
axC= fig.add_subplot(4,2,5) 
axC.plot(s1[:3000], ratio_Te_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axC.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axC.plot(s1[:3000], ratio_Te_3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
#axC.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
axC.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axC.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axC.text(2.2,0.25,"(c)")
axC.set(xlim=rlim1,
        ylim=[0,1],
        xlabel="",
        ylabel="$\\frac{T_{\\rm{e}}}{T_{\\rm{e},\\rm{std}}}$",
        title="")

axD= fig.add_subplot(4,2,7)
axD.plot(s1[:3000], ratio_delta_T_1[:3000], '-',  linewidth = width_main, label="(cx$_{\\rm{off}}$,iz$_{\\rm{off}}$)")
#axD.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axD.plot(s1[:3000], ratio_delta_T_3[:3000], '--', color="tab:green", linewidth = width_main , label="(cx$_{\\rm{on}}$,iz$_{\\rm{off}}$)")
#axD.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
axD.tick_params(which='both',axis='both', direction='in', labelbottom=True)
axD.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axD.text(2.2,0.75,"(d)")
axD.legend()
plt.legend(fontsize=10,loc='lower left',bbox_to_anchor=(1.5,-0.2))
axD.set(xlim=rlim1,
        ylim=[0,1],
        xlabel="s (m)",
        ylabel="$\\frac{(T_{\\rm{i}\\parallel}-T_{\\rm{i}\\perp})}{(T_{\\rm{i}\\parallel,\\rm{std}}-T_{\\rm{i}\\perp,\\rm{std}})}$",
        title="")

axE= fig.add_subplot(4,2,2) 
axE.plot(s1[:3000],ratio_grad_n1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axE.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axE.plot(s1[:3000],ratio_grad_n3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
#axE.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
# axE.plot(s1[:3000],grad_n1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axE.plot(s1[:3000],grad_n3[:3000], '--',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axE.plot(s1[:3000],grad_n_std[:3000], '-', color = 'black', linewidth = width_main, label="std")
axE.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axE.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axE.text(2.2,1.005,"(e)")
axE.set(xlim=rlim1,
        ylim=[0,2],
        xlabel="",
        ylabel="$\\frac{\\nabla n}{\\nabla n_{\\rm{std}}}$",
        title="")

axF= fig.add_subplot(4,2,4) 
axF.plot(s1[:3000],ratio_grad_Ti_para_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axF.plot(s1,ratio_n2, '--', linewidth = width_main, label="")
axF.plot(s1[:3000],ratio_grad_Ti_para_3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
#axF.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
# axF.plot(s1[:3000],grad_Ti_para_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axF.plot(s1[:3000],grad_Ti_para_3[:3000], '--',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axF.plot(s1[:3000],grad_Ti_para_std[:3000], '-', color = 'black', linewidth = width_main, label="std")
axF.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axF.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axF.text(2.2,4,"(f)")
axF.set(xlim=rlim1,
        ylim=[0,20],
        xlabel="",
        ylabel="$\\frac{\\nabla T_{\\rm{i}\\parallel}}{\\nabla T_{\\rm{i}\\parallel,\\rm{std}}}$",
        title="")

axG = fig.add_subplot(4,2,6) 
axG.plot(s1[:3000],ratio_grad_Te_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axG.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axG.plot(s1[:3000],ratio_grad_Te_3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
#axG.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
# axG.plot(s1[:3000],grad_Te_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axG.plot(s1[:3000],grad_Te_3[:3000], '--',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axG.plot(s1[:3000],grad_Te_std[:3000], '-', color = 'black', linewidth = width_main, label="std")
axG.tick_params(which='both',axis='both', direction='in', labelbottom=True)
axG.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axG.text(2.2,1.005,"(g)")
axG.set(xlim=rlim1,
        ylim=[1.5,2.5],
        xlabel="s (m)",
        ylabel="$\\frac{\\nabla T_{\\rm{e}}}{\\nabla T_{\\rm{e},\\rm{std}}}$",
        title="")


plt.savefig('../graphs/fig5.svg',bbox_inches="tight")
# plt.savefig('../2.3.2.n.img_pdf/grad_n,T_cx.pdf',bbox_inches="tight")
plt.show()

#figure
fig = plt.figure(figsize=(5.8, 4.8))
plt.subplots_adjust(wspace=0.43)

axA = fig.add_subplot(4,2,1) 
axA.plot(s1[:3000], n1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axA.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axA.plot(s1[:3000], n3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
axA.plot(s1[:3000], n_std[:3000], '-', color="black", linewidth = width_main , label="steady")
#axA.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
axA.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axA.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axA.text(2.2,2.1,"(a)")
axA.set(xlim=rlim1,
        ylim=[1,3],
        xlabel="",
        ylabel="$n$",
        title="")

axB = fig.add_subplot(4,2,3) 
axB.plot(s1[:3000], Ti_para_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axB.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axB.plot(s1[:3000], Ti_para_3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
axB.plot(s1[:3000], Ti_para_std[:3000], '-', color="black", linewidth = width_main , label="steady")
#axB.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
axB.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axB.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axB.text(2.2,0.25,"(b)")
axB.set(xlim=rlim1,
        ylim=[0,1],
        xlabel="",
        ylabel="$T_{\\rm{i}\\parallel}$",
        title="")
        
axC= fig.add_subplot(4,2,5) 
axC.plot(s1[:3000], Te_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axC.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axC.plot(s1[:3000], Te_3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
axC.plot(s1[:3000], Te_std[:3000], '-', color="black", linewidth = width_main , label="steady")
#axC.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
axC.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axC.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axC.text(2.2,0.25,"(c)")
axC.set(xlim=rlim1,
        ylim=[0,1],
        xlabel="",
        ylabel="$T_{\\rm{e}}$",
        title="")

axD= fig.add_subplot(4,2,7)
axD.plot(s1[:3000], diff_T_1[:3000], '-',  linewidth = width_main, label="(cx$_{\\rm{off}}$,iz$_{\\rm{off}}$)")
#axD.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axD.plot(s1[:3000], diff_T_3[:3000], '--', color="tab:green", linewidth = width_main , label="(cx$_{\\rm{on}}$,iz$_{\\rm{off}}$)")
axD.plot(s1[:3000], diff_T_std[:3000], '-', color="black", linewidth = width_main , label="steady")
#axD.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
axD.tick_params(which='both',axis='both', direction='in', labelbottom=True)
axD.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axD.text(2.2,0.75,"(d)")
axD.legend()
plt.legend(fontsize=10,loc='lower left',bbox_to_anchor=(1.5,-0.2))
axD.set(xlim=rlim1,
        ylim=[0,1],
        xlabel="s (m)",
        ylabel="$(T_{\\rm{i}\\parallel}-T_{\\rm{i}\\perp})}$",
        title="")

axE= fig.add_subplot(4,2,2) 
axE.plot(s1[:3000],grad_n_abl_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axE.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axE.plot(s1[:3000],grad_n_abl_3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
axE.plot(s1[:3000], grad_n_std_1[:3000], '-', color="black", linewidth = width_main , label="steady")
#axE.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
# axE.plot(s1[:3000],grad_n1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axE.plot(s1[:3000],grad_n3[:3000], '--',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axE.plot(s1[:3000],grad_n_std[:3000], '-', color = 'black', linewidth = width_main, label="std")
axE.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axE.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
axE.text(2.2,1.005,"(e)")
axE.set(xlim=rlim1,
        # ylim=[0,1E-3],
        xlabel="",
        ylabel="$\\nabla n$",
        title="")

axF= fig.add_subplot(4,2,4) 
axF.plot(s1[:3000],grad_Ti_para_abl_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axF.plot(s1,ratio_n2, '--', linewidth = width_main, label="")
axF.plot(s1[:3000],grad_Ti_para_abl_3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
axF.plot(s1[:3000],grad_Ti_para_std_1[:3000], '-', color="black", linewidth = width_main , label="steady")
#axF.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
# axF.plot(s1[:3000],grad_Ti_para_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axF.plot(s1[:3000],grad_Ti_para_3[:3000], '--',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axF.plot(s1[:3000],grad_Ti_para_std[:3000], '-', color = 'black', linewidth = width_main, label="std")
axF.tick_params(which='both',axis='both', direction='in', labelbottom=False)
axF.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
# axF.text(2.2,4,"(f)")
axF.set(xlim=rlim1,
        ylim=[-5E-1,5E-1],
        xlabel="",
        ylabel="$\\nabla T_{\\rm{i}\\parallel}$",
        title="")

axG = fig.add_subplot(4,2,6) 
axG.plot(s1[:3000],grad_Te_abl_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
#axG.plot(s1,ratio_n2, '--', linewidth = width_main , label="")
axG.plot(s1[:3000],grad_Te_abl_3[:3000], '--', color="tab:green", linewidth = width_main , label="cx$_{\\rm{on}}$")
axG.plot(s1[:3000],grad_Te_std_1[:3000], '-', color="black", linewidth = width_main , label="steady")
#axG.plot(s1,ratio_n4, '-.', linewidth = width_main , label="")
# axG.plot(s1[:3000],grad_Te_1[:3000], '-',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axG.plot(s1[:3000],grad_Te_3[:3000], '--',  linewidth = width_main, label="cx$_{\\rm{off}}$")
# axG.plot(s1[:3000],grad_Te_std[:3000], '-', color = 'black', linewidth = width_main, label="std")
axG.tick_params(which='both',axis='both', direction='in', labelbottom=True)
axG.vlines(2.8, 0, 10, "black", linestyle="dashdot", linewidth=width)
# axG.text(2.2,1.005,"(g)")
axG.set(xlim=rlim1,
        ylim=[-1E-1,1E-1],
        xlabel="s (m)",
        ylabel="$\\nabla T_{\\rm{e}}$",
        title="")


plt.savefig('../graphs/fig5_abs_balue.svg',bbox_inches="tight")
# plt.savefig('../2.3.2.n.img_pdf/grad_n,T_cx.pdf',bbox_inches="tight")
plt.show()
