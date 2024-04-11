import random
import os
import subprocess
import warnings
from pathlib import Path
from typing import Union, List
import logging
import time as tm
import yaml
import argparse


from scipy.special import erf
from scipy import interpolate
import pandas as pd
from matplotlib.pylab import * #FIXME
import matplotlib.pyplot as plt 

from optimorph.my_swan_function import run_swan
from optimorph.my_xbeach_function import run_xbeach
from optimorph.my_refdif_function import run_refdif, run_N_refdif
from optimorph.bathy_types import load_bathy, load_bathy_2D
from optimorph.functionnal_types import load_functionnal
from optimorph.utils_functions import (k_guo, calc_H, shoreline, shoaling_window_RONAN, slopecst, 
                                        normalise, vdot2, sandstock, create_forcing, norme, rogner,
                                        trapz_2D,func_to_plot_2D,func_to_plot_2D_flat,H_shoaling_2D,Dif_Hadamard_2D,limiteur_2D,vec_n)

warnings.filterwarnings("ignore")


######
# Params #
######

# ----- plot params
def main(output_dir: Union[str, Path], config_file: Union[str, Path], psi_datas: Union[str, Path]):
    """
    Parameters
    -----------
    output_dir : Path like
        work directory
    config_file : Path like
        yaml config file
    psi_datas : path like
        folder with all psi datas
    """
    t1 = tm.time()
    if isinstance(output_dir, str):
        output_dir = Path(output_dir)

    # ----- load params
    with Path(config_file).open("r") as f:
        param = yaml.load(f, Loader=yaml.Loader)
    logging.info(f'param {param}')
    
    # FIXME
    dir_name: str = param.get("dirname", "simu_without_geotube")
    two_dimension: bool = param.get("two_dimension", {}).get("state", False)
    debug: bool = param.get("debug", False)
    makeGifs: bool = param.get("makeGifs", True)
    figname: str = param.get("figname", 'Simulation of storm in linear bathymetry with geotube')
    output_dir = output_dir / dir_name
    logging.info(f'work dir {output_dir}')
    #print(output_dir)

    if output_dir.exists():
        shutil.rmtree(output_dir, ignore_errors=True)
        output_dir.mkdir(parents=True)
    else:
        output_dir.mkdir(parents=True)

    #Parametres de simulation a post-traiter
    fig_debug, axs_debug = plt.subplots(6, 1, figsize=(12, 10))
    fig_multi_plot, axs_multi_plot = plt.subplots(3, 2, figsize=(17, 8),dpi=300)
    fig_multi_plot.subplots_adjust(top=0.8)
    
    fig_psiH = plt.figure(figsize=(10,5),dpi=150)
    ax = fig_psiH.gca()
    theta0 = 0

    leg: [List[str]] = ['$\psi_0$ the initial bathymetry [m]','$H(x)$ the wave height from SWAN [m]','$\psi$ the final bathymetry [m]',]
    leg_hydro: [List[str]] = ["Shoaling","SWAN","XBeach"]
    laabelsize = 18
    legendsize = 16
    width = 2

    xlb1 = 'Distance from deep sea [m]'
    ylb1 = "Height $[m]$"

    # ----- hydro-morpho params
    g = 9.81
    T0: float = param.get("T0", 6)
    Hmax: float = param.get("Hmax", 2)
    h0: float = param.get("h0", 7)
    h0_ref = h0
    C0 = g / ( 2 * pi) * T0
    sigma0 = 2 * pi / T0
    gamma: float = param.get("gamma", 0.55)
    seuil = 1e-15

    niter: int = param.get("n_iteration", 100)
    ifre : int = param.get("ifre", 50)
    time = linspace(0, niter, niter + 1)
    dynamic: bool = param.get("dynamic", False)

    
    bathy_type: int = param.get("bathy_type", 0)
    logging.info(f"bathy_type: {bathy_type}")

    nwater: int = param.get("nwater", 600)
    nsand: int = param.get("nsand", 40)
    # ----- geometric params
    geotube: bool = param.get("geotube", {}).get("state", False)
    geo_pos: int = param.get("geotube", {}).get('position_x', 0)
    geo_pos_y: int = param.get("geotube", {}).get('position_y', 0)
    geo_length: int = param.get("geotube", {}).get('length', 0)
    geo_height: int = param.get("geotube", {}).get('height', 0)
    slope_max = param.get("slope_max", 7)
    mobility: float = param.get("mobility", 0.002)
    id_cost: int= param.get("id_cost_fct", 1)
    id_hydro: int = param.get("hydro_mode", 1)
    nwater_x: int = param.get("two_dimension", {}).get('nwater_x',250)
    L_x: int = param.get("two_dimension", {}).get('L_x',300)
    L_y: int = param.get("two_dimension", {}).get('L_y',20)
    direct: float = param.get("two_dimension", {}).get("theta",0)
    
    # ----- maree params
    coef_maree: float= param.get("coef_maree", 60)
    u_maree: float = param.get("u_maree", 6.1)
    maree_duration: float = param.get("maree_duration", 12.5)
    ma = coef_maree * u_maree / 100 # marnage
    total_simu_time = maree_duration * 6 # heures
    vect_t = linspace(0, total_simu_time, niter)
    h0_maree = ma / 2 * sin( 2 * pi * vect_t / maree_duration) 
    #plot(vect_t, h0_maree)
    #savefig("maree")
    
    depl_max=0.05    #max depl par iter
    depl_total_max_absolu=1  #max cumul en m
    max_loss_depth=0.5  # % perte profondeur max

    minidepth=0.02   # si<minidepth:  pas d'interaction
    h_pondera = 0.2 # ponderation from 0 to h_pondera

    histJ=[]
    sandStock = []
    filenames_psiH = []
    filenames_multi = []
    
    
    if dynamic:
        H0 = create_forcing(time, max(time)/2, 2*10**2,Hmax, -1e-4)
    else:
        H0 = ones(len(time)) * Hmax
    
    if not(two_dimension):
        
        

        
        xm = - geo_length // 2
        xM = geo_length // 2 + 1
        ntot = nwater + nsand
        x = linspace(0, ntot-1, ntot)
        xstep = x[1] - x[0]

        x, h0, ntot, psi, psif, bool_H0, HRMS = load_bathy(bathy_type, x, h0, nwater, ntot, seuil, str(psi_datas))
        h = psi - h0

        if not(bool_H0):
            if dynamic==0:
                H0 = ones(len(time))*Hmax
            else:
                H0 = create_forcing(time, max(time)/2, 2*10**2/2,Hmax, -1e-4)
                H0 = H0 / max(H0) * Hmax

        # ----- add geotube
        

        if geotube:
            xgeo = arange(xm, xM)
            add_geo = zeros(ntot)
            add_geo[ geo_pos + xm: geo_pos + xM ] = -(4*geo_height/geo_length**2)*xgeo**2+geo_height
            psi += add_geo

        psi0 = psi.copy()

        
        Mslope = slope_max * ones(ntot)

        # ----- hydro calculation
        k = k_guo(h,sigma0,g)
        H,Hpsi = calc_H(H0[0], k,h0-psi,gamma)

        # Last hydro model params
        dwin = 5*xstep
        Nwin = int(floor(dwin/xstep))
        #print(Nwin)


        ######################
        ####### Verification ########
        ######################
        
        rho0=mobility*ones(ntot)
        if geotube:
            rho0 [ geo_pos + xm : geo_pos + xM ] = 0
            Mslope  [ geo_pos + xm-2 : geo_pos + xM +2] = 1

        

        """
        if id_hydro==1: #SWAN
            #subprocess.run(["module","load","use.own"])
            #subprocess.run(["module","load","swan"])
        """

        leg[1] = "$H(x)$ the wave height from "+leg_hydro[id_hydro]+ " [m]"


        sand0 = trapz(psi0,x)

        for itera in range(niter):
            h0 = h0_ref + h0_maree[itera]
            h_t = h0-psi
            if id_hydro == 0:
                xS_t, xB_t, nbS_t, nbB_t = shoreline(h_t, x) # shoreline calculation
                #gamma = 0.55
                k_t, C_t, Cg_t, H_t, xB_t, nbB_t, H_S, tab_xB_start, tab_nbB_start, tab_xB_end, tab_nbB_end, Hw, tab_ibreak = shoaling_window_RONAN(nbS_t, nbB_t, h_t, sigma0, C0, theta0, H0[itera], x, gamma, Nwin, g) 
                C_t = C0*tanh(k_t*h_t)
                Cg_t = 0.5*C_t*(1+2*k_t*h_t/sinh(2*k_t*h_t))
            elif id_hydro == 1:
                #gamma = 0.45
                H_t, T0_t = run_swan(H0[itera], T0, psi, x, h0, gamma)   # hydro calculation
                H_t[isnan(H_t)] = 0 # put NaN to 0
                T0_t[abs(T0_t)<1e-10] = 0 # put 0 to NaN
                sigma0 = 2*pi/T0_t
                k_t = k_guo(h_t,sigma0,g)
                C_t = C0*tanh(k_t*h_t)
                Cg_t = 0.5*C_t*(1+2*k_t*h_t/sinh(2*k_t*h_t))
            elif id_hydro == 2:       
                H_t, k_t, u_t,tau_t = run_xbeach(H0[itera], T0, psi, x, h0, gamma)   # hydro calculation
                H_t[isnan(H_t)] = 0 # put NaN to 0
                tau_t = abs(tau_t)
                u_t[h_t<0] = 0
                tau_t[h_t<0] = 0 
                C_t = C0*tanh(k_t*h_t)
                Cg_t = 0.5*C_t*(1+2*k_t*h_t/sinh(2*k_t*h_t))
            else:
                raise ValueError(f"id_hydro {id_hydro} shall be in {0, 1, 2}")

            k_t[k_t==inf] = nan
            #print(max(k_t))
            dj, cost = load_functionnal(psi,H_t,id_cost)

            dj[0:2] = 0 #2*dj[1]-dj[2]

            histJ.append(cost)
            savedj = dj.copy()

            ponderatau = ones(ntot)
            #ponderatau = (1/tau_t)**1
            
            ponderatau[H_t<=h_pondera] *= - H_t[H_t<=h_pondera]
            #ponderatau/=max(ponderatau)
            #ponderatau[ponderatau<0]=0
            
            #ponderatau[h_t<=minidepth] = 0

            dumping=1/cosh(k_t*h_t)
            dumping*=(1-exp(-50*(x/ntot)**2))    #pas d'interaction si tres profond
            dumping[h_t<=minidepth]=0            #peu d'eau : pas de deplacement 
            dumping/=max(dumping)            #dumping borne par 1
            dumping[dumping<0]=0                #tjs >=0 : entre 0 et 1
            dj*=-rho0*dumping*ponderatau                  #mobilite avant projection   
            
            for i in range(1,ntot-1):    #max deplacement per iteration  avt project
                dj[i]=max(min(dj[i],depl_max),-depl_max)
                dymax=max(0,-psi[i]+psi0[i]+depl_total_max_absolu)
                dymin=min(0,-psi[i]+psi0[i]-depl_total_max_absolu)
                dj[i]=max(min(dj[i],dymax),dymin)
                dj[i]=min(dj[i],psi0[i]-psi[i]+max_loss_depth*(max(psi0)-psi0[i]))
            
            psi1 = psi + dj            # pente max implicite avt project
            psi1= slopecst(psi1, Mslope, rho0, x) 
            dj= psi1-psi

            dj[0]=0                  #bords a zero
            dj[ntot-1]=0
            dj[h_t<=minidepth]=0        #peu d'eau : pas de deplacement

            diff=trapz(psi-psi0,x)      #contrainte conservation 
            dc = psi*diff
            dc[h_t<=minidepth]=0
            dc = normalise(dc)
            dj[h_t>=minidepth] = dj[h_t>=minidepth] - vdot2(dj, dc)*dc[h_t>=minidepth]   #project avant descente 

            dj[0]=0                  #bords a zero
            dj[ntot-1]=0
            dj[h_t<=minidepth]=0        #peu d'eau : pas de deplacement 
                    
            d = dj
            
            psi += d
                
            psi[h_t>=minidepth]*=sand0/trapz(psi,x)      #contrainte conservation a posteriori
                    
            #%%
            sandStock.append(sandstock(psi, psi0, x)) # save sandstock
            if( itera%ifre == 0 or itera==niter-1 ):
            
                
                logging.info(f'iter: {itera} sandstock: {sandStock[-1]: 0.6f}') # print sandstock
            
                ##############
                # Plot psi H #
                ##############
                
                fig_psiH.clear()
                ax = fig_psiH.gca()
                
                ax.plot(x,H_t[H_t>=0]+h0_maree[itera],label="$H_S$ "+leg_hydro[id_hydro], color = "dodgerblue")
                ax.plot([0,x[-1]],[h0_maree[itera],h0_maree[itera]], color = 'cornflowerblue',linestyle=":",linewidth = 2, label ="$WL$")
                ax.plot([0,x[-1]],[0,0], color = 'purple',linestyle="-",linewidth = 3, alpha= 0.1)
                ax.plot(x,psi0-h0+h0_maree[itera], linestyle="--", color = (0.7, 0.7, 0.7) ,label= r"Initial seabed" )
                ax.plot(x,psi-h0+h0_maree[itera], label= r"Final seabed (num)", color = 'sienna')
                ax.plot(x,psif+h0_maree[itera], label= r"Final seabed (exp)",linestyle="--", color = 'red')
                ax.set_ylabel("Seabed and Wave Height",fontsize=laabelsize)
                ax.set_xlabel(xlb1,fontsize=laabelsize)
                ax.set_title('Seabed and wave evolution with parameters:' + " \n H0 = "+str(round(H0[itera],2))+" m - T0 = "+str(T0)+' s - h0 = '+str(h0_ref)+' m' ,fontsize=18, weight='bold')
                fig_psiH.text(0.75, 0.48,"Iteration = "+str(itera)+" / " + str(niter-1),size=14,rotation=0,verticalalignment='center', bbox=dict(facecolor='none', edgecolor="black", linewidth=1, boxstyle='square'))
                
                
                
                
                
                
                ax.scatter(HRMS[0]-18.180,HRMS[1],color="blue")         
                ax.legend(fontsize=legendsize,loc="best")
                xlabel(xlb1,fontsize=laabelsize)
                ax.axis([0, max(x), -h0_ref, Hmax*1.3+max(h0_maree)])
                
                ax.grid(which='major', linestyle=':',linewidth=1)
                ax.grid(which='minor', linestyle=':',linewidth=0.3)
                ax.tick_params(which='both',direction="in")



            
                ###################
                # Plot Multi Vars #
                ###################
                
                axs_multi_plot[0,0].clear()
                axs_multi_plot[1,0].clear()
                axs_multi_plot[2,0].clear()
                axs_multi_plot[0,1].clear()
                axs_multi_plot[1,1].clear()
                axs_multi_plot[2,1].clear()
                
                axs_multi_plot[0,0].grid("on")
                axs_multi_plot[1,0].grid("on")
                axs_multi_plot[2,0].grid("on")
                axs_multi_plot[0,1].grid("on")
                axs_multi_plot[1,1].grid("on")
                axs_multi_plot[2,1].grid("on")
                
                axs_multi_plot[0,0].grid(which='major', linestyle=':',linewidth=1)
                axs_multi_plot[1,0].grid(which='major', linestyle=':',linewidth=1)
                axs_multi_plot[2,0].grid(which='major', linestyle=':',linewidth=1)
                axs_multi_plot[0,1].grid(which='major', linestyle=':',linewidth=1)
                axs_multi_plot[1,1].grid(which='major', linestyle=':',linewidth=1)
                axs_multi_plot[2,1].grid(which='major', linestyle=':',linewidth=1)
                
                axs_multi_plot[0,0].grid(which='minor', linestyle=':',linewidth=0.3)
                axs_multi_plot[1,0].grid(which='minor', linestyle=':',linewidth=0.3)
                axs_multi_plot[2,0].grid(which='minor', linestyle=':',linewidth=0.3)
                axs_multi_plot[0,1].grid(which='minor', linestyle=':',linewidth=0.3)
                axs_multi_plot[1,1].grid(which='minor', linestyle=':',linewidth=0.3)
                axs_multi_plot[2,1].grid(which='minor', linestyle=':',linewidth=0.3)
                
                axs_multi_plot[0,0].tick_params(which='both',direction="in")
                axs_multi_plot[1,0].tick_params(which='both',direction="in")
                axs_multi_plot[2,0].tick_params(which='both',direction="in")
                axs_multi_plot[0,1].tick_params(which='both',direction="in")
                axs_multi_plot[1,1].tick_params(which='both',direction="in")
                axs_multi_plot[2,1].tick_params(which='both',direction="in")
                
                fig_multi_plot.suptitle(figname + " \n Parameters: H0 = "+str(round(H0[itera],2))+" m - T0 = "+str(T0)+' s - h0 = '+str(h0_ref)+' m' '\n Iteration = '+str(itera)+" / " + str(niter-1) ,fontsize=18, weight='bold', y=0.98)
                
                
                # Figure 1 : Wave forcing over time
                        
                axs_multi_plot[0,0].plot(time, H0, label= r'Wave height at $x=0$', color = 'dodgerblue') 
                axs_multi_plot[0,0].scatter(itera,[H0[itera]], color = 'royalblue', label = 'Current forcing wave height')
                
                axs_multi_plot[0,0].set_xlabel("Iteration",fontsize = 'small'  )
                axs_multi_plot[0,0].set_ylabel("Wave height at $x$ = 0",fontsize = 'small')
                axs_multi_plot[0,0].legend(loc="best")
                axs_multi_plot[0,0].axis([0, max(time), 0, max(H0)*1.3])
                
                # Figure 2 : Wave characteritics over the cross-shore profile


                maxk = max(k_t)
                maxC = max(C_t)
                maxCg = max(Cg_t)

                colork = "maroon"
                colorC = "mediumseagreen"
                colorCg = "lightgreen"
                colorspeed = "green"
                
                ax2 = axs_multi_plot[1,0].twinx()
                
                
                axs_multi_plot[1,0].plot(x, k_t, color=colork, label = r"$k$")	
                axs_multi_plot[1,0].axis([0, max(x), 0, maxk*1.1])
                axs_multi_plot[1,0].set_ylabel(r"Wave number ($m^{-1}$)", fontsize = 'small')
                axs_multi_plot[1,0].tick_params(axis='y', labelcolor=colork)
                axs_multi_plot[1,0].plot(nan, color=colorC, label= r"$C$")
                axs_multi_plot[1,0].plot(nan, color=colorCg, label= r"$C_g$")
                axs_multi_plot[1,0].legend(loc="upper right")
                
                axs_multi_plot[1,0].set_xlabel("Distance from deep sea [m]", fontsize = 'small')
                
                
                
                ax2.set_ylabel(r"Velocity ($m.s^{-1}$)", fontsize = 'small')
                ax2.plot(x, C_t, color=colorC)	
                ax2.plot(x, Cg_t, color=colorCg)
                ax2.axis([0, max(x), 0, max(maxC,maxCg)*1.15])
                ax2.tick_params(axis='y', labelcolor=colorspeed)        
                
                # Figure 3 : Sand mobility over the cross-shore profile
                rhomax = max(rho0)
                color1 = 'olivedrab'
                color2 = 'darkred'
                
                axs_multi_plot[2,0].plot(x, rho0, label= r"Sand mobility", color = color1)
                axs_multi_plot[2,0].plot(nan, color=color2, label= r"Maximal slope")

                axs_multi_plot[2,0].set_xlabel("Distance from deep sea [m]", fontsize = 'small')
                axs_multi_plot[2,0].set_ylabel("Sand mobility ($m.s.kg^{-1}$)", fontsize = 'small')
                
                axs_multi_plot[2,0].tick_params(axis='y', labelcolor=color1)
                axs_multi_plot[2,0].legend(loc="best")
                axs_multi_plot[2,0].axis([0, max(x), 0, max(rho0)*1.1])
                
                ax3 = axs_multi_plot[2,0].twinx()
                ax3.set_ylabel("Maximal seabed slope", fontsize = 'small')
                ax3.plot(x, Mslope, color=color2, label= r"Maximal slope")
                ax3.set_yscale('log')
                ax3.tick_params(axis='y', labelcolor=color2)
                
                # Figure 4 :  Wave height over the cross-shore profile
                        
                axs_multi_plot[0,1].plot(x, H_t, label = r"$Hs$ "+leg_hydro[id_hydro],color = 'dodgerblue')

                axs_multi_plot[0,1].set_xlabel("Distance from deep sea [m]", fontsize = 'small')
                axs_multi_plot[0,1].set_ylabel("Wave height (m)", fontsize = 'small')
                axs_multi_plot[0,1].legend(loc="best")
                axs_multi_plot[0,1].axis([0, max(x), 0, max(H0)*1.3])

                # Figure 5 : Gradient calculation
                
                axs_multi_plot[1,1].plot(x,savedj, color = "red", label="$d$ descent direction")

                axs_multi_plot[1,1].set_xlabel("Distance from deep sea [m]", fontsize = 'small')
                axs_multi_plot[1,1].set_ylabel("$d$ descent direction", fontsize = 'small')
                axs_multi_plot[1,1].legend(loc="best")
                amplitude_d = abs(max(savedj) - min(savedj))
                axs_multi_plot[1,1].axis([0, max(x), min(savedj)-amplitude_d*0.1, max(savedj)+amplitude_d*0.1])
                
                # Figure 6 : Seabed / Bedrock / Mean Water Level over the cross-shore profile 
                
                axs_multi_plot[2,1].plot(x, psi0-h0+h0_maree[itera], linestyle="--", color = (0.7, 0.7, 0.7) ,label= r"Initial seabed" )
                axs_multi_plot[2,1].plot(x, psi-h0+h0_maree[itera], label= r"Final seabed", color = 'sienna')    
                axs_multi_plot[2,1].plot([0,x[-1]],[h0_maree[itera],h0_maree[itera]], label=r'$WL$', linestyle=":", color = 'cornflowerblue',linewidth = 2)
                axs_multi_plot[2,1].plot([0,x[-1]],[0,0], color = 'purple',linestyle="-",linewidth = 3, alpha= 0.1)

            

                axs_multi_plot[2,1].set_xlabel("Distance from deep sea [m]", fontsize = 'small')
                axs_multi_plot[2,1].set_ylabel("Seabed elevation (m)", fontsize = 'small')
                axs_multi_plot[2,1].legend(loc="best")
                axs_multi_plot[2,1].axis([0, max(x), -h0_ref,  Hmax*1.3+max(h0_maree)])
                
                
                                
                
                
            
                ##############
                # Plot debug #
                ##############
                
                if debug:
                    axs_debug[0].clear()
                    axs_debug[1].clear()
                    axs_debug[2].clear()
                    axs_debug[3].clear()
                    axs_debug[4].clear()
                    axs_debug[5].clear()
                    
                    axs_debug[0].set_title('Itération = '+str(itera)+" - H0 = "+str(round(H0[itera],2))+" m - T0 = "+str(T0)+' s - h0 = '+str(h0)+' m')
                    

                    axs_debug[0].plot(dumping)
                    axs_debug[1].plot(savedj)
                    axs_debug[2].plot(d)
                    axs_debug[3].plot(psi-h0,color="black",label="Final Model")
                    axs_debug[4].plot(ponderatau)
                    axs_debug[5].plot(ponderatau*dumping)
                            
                    
                    axs_debug[0].grid("on")
                    axs_debug[1].grid("on")
                    axs_debug[2].grid("on")
                    axs_debug[3].grid("on")
                    axs_debug[4].grid("on")
                    axs_debug[5].grid("on")
                    
                    axs_debug[0].set_ylabel("fonction attenuation verticale")
                    axs_debug[1].set_ylabel("$dj$")
                    axs_debug[2].set_ylabel("$d$")
                    axs_debug[3].set_ylabel("$H$ and $\psi$")
                    axs_debug[4].set_ylabel("ponderaH")
                    axs_debug[5].set_ylabel("Upsilon *ponderaH")

                    
                    axs_debug[3].plot(H_t[H_t>=0],label="H "+leg_hydro[id_hydro])
                    axs_debug[3].plot([0,x[-1]],[0,0],color="blue",linestyle="dashed")
                    axs_debug[3].plot(psi0-h0,color="black",linestyle="dashed",label="Initial")
                    axs_debug[3].plot(psif-h0,color="red")
                    axs_debug[3].legend()
                    
                    fig_debug.tight_layout()
                    fig_debug.savefig(output_dir / f"debug_iter_{itera}.png")





                #axs_debug[3].axis([x[psi<h0-20][-1], x[-1], psi[psi<h0-20][-1], max(H_t)*1.2+h0])
                            
                
                fig_psiH.tight_layout()
                fig_psiH.savefig(output_dir / f"psiH_iter_{itera}.png")
                filenames_psiH.append(output_dir / f"psiH_iter_{itera}.png")
                
                #fig_multi_plot.tight_layout(rect=[0, 0, 1, 0.9])
                fig_multi_plot.tight_layout()
                fig_multi_plot.savefig(output_dir / f"multi-plot_iter_{itera}.png")
                ax2.clear()
                ax2.axes.get_yaxis().set_visible(False)
                filenames_multi.append(output_dir / f"multi-plot_iter_{itera}.png")
                
                savetxt(output_dir / f'psi_{itera}.dat', psi)
                savetxt(output_dir / f'H_{itera}.dat', H_t)
                
            #    print('iter=',itera)
                
            #%%
            
        # plot Jhist
        figure()
        plot(histJ,color="green")
        xlabel("Iterations",weight='bold',fontsize = 14)
        ylabel('Jhist',weight='bold',fontsize = 14)
        title("Historical values of J with respect to iterations",weight='bold',fontsize = 16)
        grid(which='major', linestyle=':',linewidth=1)      
        grid(which='minor', linestyle=':',linewidth=0.3)
        tick_params(which='both',direction="in")
        
        tight_layout()
        savefig(output_dir / f'Jhist={itera}.png')
        savetxt(output_dir / "histJ.txt", histJ)
        
        # plot SandStock
        figure()
        plot(sandStock,color="red")
        xlabel("Iterations",weight='bold',fontsize = 14)
        ylabel('Sandstock conservation',weight='bold',fontsize = 14)
        title("SandStock conservation with respect to iterations",weight='bold',fontsize = 15)
        grid(which='major', linestyle=':',linewidth=1)      
        grid(which='minor', linestyle=':',linewidth=0.3)
        tick_params(which='both',direction="in")
        yscale('log')
        #print(sandStock)
        
        tight_layout()

        savefig(output_dir / f'sandStock={itera}.png')
        savetxt(output_dir / "sandStock.txt", sandStock)
        t2 = tm.time()
        logging.info('----------------------------------------------------') # 
        logging.info('----------------- Calculation done -----------------')
        logging.info('--------------- Computing time= '+str(round(t2-t1))+" s ---------------") 
        logging.info('----------------------------------------------------') #
        if makeGifs:
            import imageio
            from PIL import Image
            images_psiH = []
            images_multi = []
            for filename in filenames_psiH:
                images_psiH.append(imageio.imread(filename))
            for filename in filenames_multi:
                images_multi.append(imageio.imread(filename))
            imageio.mimsave(str(output_dir) + "/psiH"+'.gif', images_psiH)
            imageio.mimsave(str(output_dir) + "/multi_plot"+'.gif', images_multi)
            
            
    else:
        
        

        
        print('launching calculation in 2D')


        fig = figure(figsize=(17,10))
        fig.tight_layout(rect=[0.4, 0.4, 0.6, 0.6])
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        
        fig_multi_plot, axs_multi_plot = plt.subplots(1, 2, figsize=(17, 8),dpi=300)
        
        print(L_x,L_y)
        
         
        x = arange(0, L_x, 1)
        y = arange(0, L_y, 1)
        
        A = geo_height
        a = 0.05
        x_pos= geo_pos
        y_pos = geo_pos_y
        #print(x_pos,y_pos)
        #bathy_type = 1
        nwater = nwater_x
        #h0 = 7
        
        X, Y, h0, n_i, n_j, psi = load_bathy_2D(bathy_type, x, y, h0, nwater_x, seuil, psi_datas)
        
        #X, Y = np.meshgrid(x,y)
        gauss = A*exp(-a*((X-x_pos)**2+(Y-y_pos)**2))
        to_block= X[gauss>0.1]
#        x_block_inf, x_block_sup = int(to_block[0]), int(to_block[-1])+1
        #psi = 0.01*X+0.05*Y
        
        if geotube:
            psi +=  gauss
        psi0 = psi.copy()
        Nwin = 1
        
        #
        
        #
        Ntheta = 20
        multi_h = False

        if id_hydro == 3:
            if multi_h == True:
                H = run_N_refdif(psi, h0, H0[0], T0, direct, Ntheta)
            else:
                H = run_refdif(psi, h0, H0[0], T0, direct)
        else:
            H, K = H_shoaling_2D(psi, h0, x, sigma0, C0, H0[0], gamma, Nwin, g)

        
        
        
        func_to_plot_2D(psi, H, X, Y, h0, 0, niter, ax)
        #show()
        #exit()
        Hpsi = Dif_Hadamard_2D (psi, H, X , Y)
        
        rho0=mobility*ones(( n_j , n_i ))
        
        histJ=[]
        sand0= trapz_2D(x, y, psi)
        time = linspace (0,niter,niter+1)
        
        for itera in range(niter):
            h_t = h0-psi
            #H_t, k_t = H_shoaling_2D(psi, h0, x, sigma0, C0, H0[itera], gamma, Nwin, g)
            #H_t = run_refdif(psi, h0, H0[0], T0, direct)
            if id_hydro == 3:
                if multi_h == True:
                    H_t = run_N_refdif(psi, h0, H0[0], T0, direct, Ntheta)
                else:
                    H_t = run_refdif(psi, h0, H0[0], T0, direct)
            else:
                H_t, k_t = H_shoaling_2D(psi, h0, x, sigma0, C0, H0[0], gamma, Nwin, g)

            k_t = k_guo(h_t,sigma0,g)
            
            H_t[isnan(H_t)] = 0 # put NaN to 0
            H_psi_num_N = Dif_Hadamard_2D (psi, H_t, X , Y)
            H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
            cost=0.5*norme(H_t)**2
            dj   = H_t*H_psi_num_N    #grad J

            dj[:,0:2]=0 #2*dj[1]-dj[2]
            dj[0,:]=0
            dj[n_j-1,:]=0

            histJ.append(cost)
            savedj = dj.copy()

            dumping=1/cosh(k_t*h_t)
            dumping*=(1-exp(-50*(X/n_i)**2))     #pas d'interaction si tres profond
            dumping[h_t<=minidepth]=0              #peu d'eau : pas de deplacement 
            dumping/=np.max(dumping)               #dumping borne par 1
            dumping[dumping<0]=0                   #tjs >=0 : entre 0 et 1

            
            dj*=-rho0*dumping                      #mobilite avant projection   
            
            for j in range(1,n_i-1):    #max deplacement per iteration  avt project
                for i in range(1,n_j-1):
                     dj[i,j]=max(min(dj[i,j],depl_max),-depl_max)
                     dymax=max(0,-psi[i,j]+psi0[i,j]+depl_total_max_absolu)
                     dymin=min(0,-psi[i,j]+psi0[i,j]-depl_total_max_absolu)
                     dj[i,j]=max(min(dj[i,j],dymax),dymin)
                     dj[i,j]=min(dj[i,j],psi0[i,j]-psi[i,j]+max_loss_depth*(np.max(psi0)-psi0[i,j]))
            
            psi1 = psi + dj              # pente max implicite avt project
            #psi1= slopecst(psi1, Mslope, rho0, x) 
            dj= psi1-psi

            dj[:,0]=0                      #bords a zero
            dj[:,n_i-1]=0
            dj[h_t<=minidepth]=0         #peu d'eau : pas de deplacement

            diff=trapz_2D(x, y, psi-psi0)       #contrainte conservation 
            dc = psi*diff
            dc[h_t<=minidepth]=0
            dc = normalise(dc)
            dj[h_t>=minidepth] = dj[h_t>=minidepth] - vdot2(dj, dc)*dc[h_t>=minidepth]   #project avant descente 

            dj[:,0]=0                      #bords a zero
            dj[:,n_i-1]=0
            dj[h_t<=minidepth]=0         #peu d'eau : pas de deplacement
                    
            d = dj
            
            psi += d
                
            psi[h_t>=minidepth]*=sand0/trapz_2D(x, y, psi)       #contrainte conservation a posteriori
            sandStock.append(trapz_2D(x, y, psi)/sand0) # save sandstock
            
            if( itera%ifre == 0 or itera==niter-1 ):
                 logging.info(f'iter: {itera} sandstock: {sandStock[-1]: 0.6f}') # print sandstock

                 func_to_plot_2D(psi, H_t, X, Y, h0, itera, niter, ax)
                 fig.tight_layout()
                 fig.savefig(output_dir / f"psiH-2D_iter_{itera}.png")
                 rogner(output_dir / f"psiH-2D_iter_{itera}")
                 filenames_psiH.append(output_dir / f"psiH-2D_iter_{itera}.png")

                 cb1, cb2 = func_to_plot_2D_flat(psi, H_t, X, Y, h0, L_x, L_y, itera, niter, axs_multi_plot)
                 fig_multi_plot.suptitle('Wave height $H$ and Sea Bottom Profile $\psi$ \n at iteration '+str(itera+1)+' / '+str(niter), weight="bold",fontsize = 20)
                 fig_multi_plot.tight_layout()
                 fig_multi_plot.savefig(output_dir / f"multi-plot-2D_iter_{itera}.png")
                 cb1.remove()
                 cb2.remove()
                 #filenames_psiH.append(output_dir / f"psiH-2D_iter_{itera}.png")
                 
                 
                 savetxt(output_dir / f'psi_{itera}.dat', psi)
                 savetxt(output_dir / f'H_{itera}.dat', H_t)
                 #rogner(ttlf)
                 #savefig(ttlf+".svg")
            t2 = tm.time()
            logging.info('Calculation done, computing time= '+str(round(t2-t1))+" s") # 
        if makeGifs:
            import imageio
            from PIL import Image
            images_psiH = []
            for filename in filenames_psiH:
                images_psiH.append(imageio.imread(filename))
            imageio.mimsave(str(output_dir) + "/psiH-2D"+'.gif', images_psiH)
        
        






        
        
        
        







    # -------------- POST TRAITEMENT ----------------------
    """
    ttl1="Wave height $H(x)$ from "+leg_hydro[id_hydro] + " and bathymetric profil after "+str(niter)+" itérations - T0= "+str(T0)+" s - h0 = "+str(h0)+' m'
    ttlfile= "OptiSwan_T0="+str(T0)+"-Hmax="+str(Hmax)+"-Omega="+str(nwater)+'-h0='+str(h0)
    fig=figure(figsize=(14,6))
    F = pd.read_csv(filename_save,header=6,names=['Hsig'])
    H = array(F.Hsig)
    H[H<0] = 0

    plot(x,H+h0,linestyle=lstyle[0],linewidth=width,label=leg[1],color=color[3])
    plot(x,psi0,linestyle=lstyle[3],linewidth=width,label=leg[0],color=color[5])
    plot(x,loadtxt("psi.dat")+h0,linestyle=lstyle[0],linewidth=width,label=leg[2],color="brown")
    plot(x,h0*ones(len(x)),linestyle=lstyle[1],linewidth=width,label="$h_0$",color="blue")



    xlabel(xlb1,fontsize=laabelsize)
    ylabel(ylb1,fontsize=laabelsize)
    title(ttl1, fontsize=titlesize)
    legend(fontsize=legendsize,loc="best")

    grid(visible=True, which='major', linestyle=':',linewidth=1)
    tick_params(which='both',direction="in")	


    savefig(repertoire+'/'+ttlfile+".png",format='png')
    savefig(repertoire+'/'+ttlfile+".svg",format='svg')

    """


    #############################################################################
    #coder en detail pour verif
    # dc=psi.copy()   
    # diff=0
    # for i in range(ntot):
    #    diff+=psi[i]-psi0[i]
    # Ndc=0
    # for i in range(ntot):
    #    dc[i]=max(diff*psi[i],0)
    #    if(h_t[i]<=0):
    #        dc[i]=0
    #    Ndc+=dc[i]**2
    # dc/=max(Ndc,1.e-30)
    # ps=0
    # for i in range(ntot):
    #    ps+=dj[i]*dc[i]
    # dj-=ps*dc
    # dj[0]=0                  #bords à zero
    # dj[ntot-1]=0
    # dj[h_t<=minidepth]=0        #peu d'eau : pas de deplacement 
    # psi = psi1 - dj
    ##############################################################################

if __name__ == "__main__":
    # TODO remove
    import shutil
    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('-i', '--input', dest='config_file', default='user_config.yaml', help='Input configuration file')
    args = parser.parse_args()
    logging.getLogger('matplotlib').setLevel(logging.INFO)
    logging.basicConfig(level=logging.INFO)
    output_dir = Path(__file__).parents[1] / "Results"
    #if output_dir.exists():
        #shutil.rmtree(output_dir, ignore_errors=True)
    main(
        output_dir=Path(__file__).parents[1] / "Results", 
        config_file= Path(__file__).parents[1] / args.config_file, 
        psi_datas = Path(__file__).parents[1] / "datas_psi/"
    )
