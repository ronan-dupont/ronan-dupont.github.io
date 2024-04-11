from matplotlib.pylab import *

import os

import subprocess

import pandas as pd

from random import random

def rd(a, b):
    return (b-a) * random() + a

def gaussian(x, mu, sig):
    return (1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2))

"""
theta = 0
lims = 15
sig = 8
x = linspace (2*( theta - lims) , 2 * (theta + lims), 1000)
plot(x, gaussian( x, theta, sig))
show()
run_N_refdif(0, 0, 0, 0, 0, 10)
"""

def run_refdif(psi, h0, H0, T_0, direct):

    psi_refdif = h0-psi

    path_folder = "/home/ronan/Bureau/TP_REFDIF/refdif-gcl/"
    script_folder = "/home/ronan/Bureau/TP_REFDIF/Python_scripts/"
    
    mr = len(psi[0,:])
    nr = len(psi[:,0])
    
    dxr, dyr = 1 , 1



    freqs= T_0
    nwavs= 1
    iff1= 1
    iff2= 1 
    iff3= 1
    amp = H0/2


    sl = "\n"
    f = open(path_folder + "datas.input","w")
    f.write("mr= "+str(mr)+ sl)
    f.write("nr= "+str(nr)+ sl)
    f.write("dxr= "+str(dxr)+ sl)
    f.write("dyr= "+str(dyr)+ sl)
    f.write("freqs= "+str(freqs)+ sl)
    f.write("dir= "+str(direct)+ sl)
    f.write("nwavs= "+str(nwavs)+ sl)
    f.write("iff1= "+str(iff1)+ sl)
    f.write("iff2= "+str(iff2)+ sl)
    f.write("iff3= "+str(iff3)+ sl)
    f.write("amp= "+str(amp)+ sl)
    f.close()


    np.savetxt(path_folder+"refdat.dat", psi_refdif.T, fmt="%20.10f", delimiter='')
    #plt.show()

    os.chdir(path_folder)
    subprocess.run('./make-indat1')
    subprocess.run('./refdif2')
    
    #from PT_main import main_PT
    #height = main_PT()
    df = pd.read_csv(path_folder+"datas.input", sep='=', index_col=0, header=None)
    x = np.arange(0, df.loc['mr'].values) * df.loc['dxr'].values
    y = np.arange(0, df.loc['nr'].values) * df.loc['dyr'].values
    X, Y = np.meshgrid(x, y)
    height = np.loadtxt(path_folder+"height.dat").T  
    
    """

    fig, ax = plt.subplots(1,1, dpi=300)
    ax.set_title("Wave height in meters")
    scb = ax.pcolor(X, Y, height, cmap="autumn_r")
    fig.colorbar(scb)
    ax.set_aspect(1)
    show()
    """
    
    return height


def run_N_refdif(psi, h0, H0, T_0, direct, Ntheta):
    mu = direct
    lims = 10
    sig = 5 # sigma grand => Ouverture de gaussienne grande
    theta_i = array([ rd( mu - lims, mu + lims ) for i in range(Ntheta)])
    p_i = gaussian( theta_i, mu, sig)
    S_p = p_i.sum()
    p_i /= S_p # Normalise, now Sum p_i = 1
    Htot = zeros_like(psi)
    #print(Htot)
    for i in range(Ntheta):
        theta = theta_i[i]
        p = p_i[i]
        H = run_refdif(psi, h0, H0, T_0, theta)
        Htot += p * H
    return Htot
