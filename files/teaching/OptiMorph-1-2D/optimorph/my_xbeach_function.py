from matplotlib.pylab import *
from scipy import io
import os
import pandas as pd
from scipy import interpolate
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import shutil
import stat
from zipfile import ZipFile
import subprocess
import xarray as xr

def str_replace(n):
    string = "{}"
    for i in range(n-1):
       string+=' {}'
    return string


def run_xbeach(H0, T0, psi, x, h0, gamma):

# -------------------------------------- INITIALISATION ------------------------------------
    filename_read = "params_ref.txt"
    filename_write = "params.txt"
    bathy_name = "bathy_psi.dep"
    output_name = 'xboutput.nc'
    x_name = "x.grd"
    y_name = "y.grd"

    ntot = len(psi)
    
    fichier_read = open(filename_read, "r")
    fichier_write = open(filename_write, "w")
    
    savetxt(bathy_name,psi - h0,newline = '  ')
    savetxt(x_name, x,newline = '  ')
    savetxt(y_name, zeros(ntot), newline = '  ')
    

    dx = x[1] - x[0]
    i = 0
    save_i = 0

# -------------------------------------- CREATE INIT FILE -----------------------------------------

    for line in fichier_read:
          if (line[0:7])=="depfile":
                split_line = line.split()
                split_line[2] = bathy_name
                fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
                i += 1
          elif (line[0:2])=="nx":
                split_line = line.split()
                split_line[2] = ntot-1
                fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
                i += 1
          elif (line[0:5])=="xfile":
                split_line = line.split()
                split_line[2] = x_name
                fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
                i += 1
          elif (line[0:5])=="yfile":
                split_line = line.split()
                split_line[2] = y_name
                fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
                i += 1
          elif (line[0:4])=="Hrms":
                split_line = line.split()
                split_line[2] = H0
                fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
                i += 1
          elif (line[0:4])=="Trep":
                split_line = line.split()
                split_line[2] = T0
                fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
                i += 1
          elif (line[0:6])=="gammax":
                split_line = line.split()
                split_line[2] = gamma
                fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
                i += 1
          else:
                fichier_write.write(line)
                i += 1
                
    fichier_write.close()   

    # -------------------------------------- RUN CALCULATION -----------------------------------------------------

    subprocess.run(["xbeach",filename_write])
    print("calculation launched and done")
    print("------------------------------------")
    print("Post-traitement")
# -------------------------------------- POST CALCULATION -----------------------------------------------------
    print("------------------------------------")

    Read_nc= xr.load_dataset(output_name)
    H = array(Read_nc.isel(globaltime=-1).H.data[0])
    k = array(Read_nc.isel(globaltime=-1).k.data[0])
    u = array(Read_nc.isel(globaltime=-1).u.data[0])
    taubx = array(Read_nc.isel(globaltime=-1).taubx.data[0])
    H[H<0] = 0
    k[k<0] = 0
    print("Post-traitement Done")
    return H, k, u, taubx
