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

def str_replace(n):
    string = "{}"
    for i in range(n-1):
       string+=' {}'
    return string


def run_swan(H0, T0, psi, x, h0, gamma):

# -------------------------------------- INITIALISATION ------------------------------------
    filename_read = "maupiti1D_1m.swn"
    filename_write = "INPUT.swn"
    bathy_name = "psi.dat"
    filename_save_HSIG = "swan_output_HSIG.dat"
    filename_save_T0 ="swan_output_T0.dat"
    fichier_read = open(filename_read, "r")
    fichier_write = open(filename_write, "w")
    savetxt(bathy_name,psi - h0,newline = '  ')
    ntot = len(psi)
    dx = x[1] - x[0]
    i = 0
    save_i = 0

# -------------------------------------- CREATE INIT FILE -----------------------------------------

    for line in fichier_read:
       if (line[0:5])=="CGRID":
          split_line = line.split()
          split_line[5] = ntot
          split_line[7] = int((ntot-1)//dx)
          fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
          i += 1

       elif (line[0:5])=="INPGR":
          split_line = line.split()
          split_line[6] =int(ntot//dx-1)
          split_line[8] = int(dx)  
          fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
          i += 1

       elif (line[0:14])=="READINP BOTTOM":
          split_line = line.split()
          split_line[3] = "'"+bathy_name+"'"
          fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
          i += 1

       elif (line[0:19])=="BOUNDSPEC SIDE West":
          print('ok')
          split_line = line.split()
          split_line[5] = H0
          split_line[6] = T0
          fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
          i += 1
       elif (line[0:5])=="CURVE":
          split_line = line.split()
          split_line[4] = int(ntot//dx-1)
          split_line[5] = int(ntot//dx-1)
          fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
          i += 1
       elif (line[0:8])=="BREAKING":
          split_line = line.split()
          split_line[3] = gamma
          fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
          i += 1
       elif (line[0:5])=="TABLE" and i!=save_i+1 :
          split_line = line.split()
          split_line[3] = "'"+filename_save_HSIG+"'"
          fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
          save_i = i
          i += 1
       elif (line[0:5])=="TABLE" and i==save_i+1:
          print("boucle")
          split_line = line.split()
          split_line[3] = "'"+filename_save_T0+"'"
          fichier_write.write(str_replace(len(split_line)).format(*split_line)+"\n")
          save_i = i
          i += 1

          """


       elif (line[0:7])=="COMPUTE":
          split_line = line.split()
          split_line[2] = dt
          split_line[4] = float(tf)
          fichier_write.write("{} {} {} {} {:08.8} {}  {}  {}  {}".format(*split_line)+"\n")
          """
       else:
          fichier_write.write(line)
          i += 1
          
    fichier_write.close()      

    # -------------------------------------- RUN CALCULATION -----------------------------------------------------

    subprocess.run(["swanrun","-input",filename_write])
    print("calculation launched and done")
    print("------------------------------------")
    print("Post-traitement")
# -------------------------------------- POST CALCULATION -----------------------------------------------------
    print("------------------------------------")

    H = array(pd.read_csv(filename_save_HSIG,header=6,names=['Hsig']).Hsig)
    T0 = array(pd.read_csv(filename_save_T0,header=6,names=['Tm01']).Tm01)
    H[H<0] = 0
    T0[T0<0] = 0
    print("Post-traitement Done")
    return H,T0


"""
File = io.loadmat(filename_save)
figure(figsize=(12,8))
X= File["Xp"][0]
HS= File["Hsig"]
plot(X,HS[0],"o",markersize=2,label="$H_s$ no interpolate",color="red")
Bathy= File["Botlev"]
text(1,0.5,"MWL",c="blue")
plot([0,L+1],[0,0],color="blue",linestyle="-.")
xlabel("Length [m]")
ylabel("Wave hight and Seabed $\psi$ [m]")
title("$H_s$ with the configuration of $T_0$="+str(T0)+" $s$ and $H_0$="+str(H0)+ " $m$")



# -------------------------------------- INTERPOLATION TO GRID -----------------------------------------------------

X_interp = arange (0,L+1,dx_interp)
f = interpolate.interp1d(X, HS)
HS_interp = f (X_interp)[0]
#print(X_interp)
#print(HS_interp)
plot(X_interp,HS_interp,label="$H_s$ interpolate",color="darkviolet")
plot(X,-Bathy[0],color="black",label="$\psi$")
legend()

# -------------------------------------- SAVE FIG -----------------------------------------------------

savefig("test.png")
print("testok")
"""
