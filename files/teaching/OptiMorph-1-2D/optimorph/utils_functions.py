from matplotlib.pylab import *
from scipy import interpolate
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import os as os
from scipy.special import erf
import cv2


def dif_LaxWendroff(alpha, psi):
    dif_ordre_2 = zeros_like(psi)
    dif_centre = zeros_like(psi)
    dif_centre[1:-1] = psi[2:] - psi[0:-2]
    dif_ordre_2[1:-1] = psi[2:] - 2 * psi[1:-1] + psi[0:-2]
    return alpha / 2 * dif_centre - alpha**2/2 * dif_ordre_2
    

def calc_ponderatau(ntot, H_t, h_t, h_pondera, minidepth):
    ponderatau = ones(ntot)
    ponderatau[H_t<=0.2] *=  H_t[H_t<=0.2] #to uncomment
    ponderatau/=max(ponderatau)        
    ponderatau[h_t<=minidepth] = 0
    return ponderatau
    
def calc_dumping(ntot, x, k_t, h_t, minidepth):
    dumping=1/cosh(k_t*h_t)
    dumping*=(1-exp(-50*(x/ntot)**2))    #pas d'interaction si tres profond
    dumping[h_t<=minidepth]=0            #peu d'eau : pas de deplacement 
    dumping/=max(dumping)            #dumping borne par 1
    dumping[dumping<0]=0                #tjs >=0 : entre 0 et 1
    return dumping

def calc_dumping_inv(ntot, x, k_t, h_t, minidepth):
    dumping=1/cosh(k_t*h_t)
    dumping*=(1-exp(-50*(x/ntot)**2))    #pas d'interaction si tres profond
    dumping = dumping[::-1]
    dumping[h_t<=minidepth]=0            #peu d'eau : pas de deplacement 
    dumping/=max(dumping)            #dumping borne par 1
    dumping[dumping<0]=0                #tjs >=0 : entre 0 et 1
    return dumping
    
def dif_centre_neuman(psi, dx):
    psi_x = zeros_like(psi)
    n_x = len(psi)
    for i in range(1,n_x-1):
        psi_x[i] = (psi[i+1] - psi[i-1]) / (2*dx)
    psi_x[0] = psi_x[1]
    psi_x[n_x-1]=psi_x[n_x-2]
    return lim_pente_1(psi_x)

def calc_vect_t(psi,dx):
    t_x = zeros_like(psi)
    t_y = zeros_like(psi)
    dy = zeros_like(psi)
    n_x = len(psi)
    for i in range(1,n_x-1):
        dy[i] = (psi[i+1]-psi[i-1])
        t_x[i] = dx / (dx**2+dy[i]**2)**0.5
        t_y[i] = dy[i] / (dx**2+dy[i]**2)**0.5 
    t_y[0] = t_y[1]
    t_y[n_x-1] = t_y[n_x-2]
    t_x[0] = t_x[1]
    t_x[n_x-1] = t_x[n_x-2] 
    return lim_pente_1(t_x), lim_pente_1(t_y)

def apply_max_displacement(dj, psi, psi0, ntot, depl_max, depl_total_max_absolu, max_loss_depth):
    for i in range(1,ntot-1):    #max deplacement per iteration  avt project
        dj[i]=max(min(dj[i],depl_max),-depl_max)
        dymax=max(0,-psi[i]+psi0[i]+depl_total_max_absolu)
        dymin=min(0,-psi[i]+psi0[i]-depl_total_max_absolu)
        dj[i]=max(min(dj[i],dymax),dymin)
        dj[i]=min(dj[i],psi0[i]-psi[i]+max_loss_depth*(max(psi0)-psi0[i]))
    return dj

def apply_sand_conservation(dj, psi, psi0, x, h_t, minidepth):
    diff=trapz(psi-psi0,x)      #contrainte conservation 
    dc = psi*diff
    dc[h_t<=minidepth]=0
    dc = normalise(dc)
    dj[h_t>=minidepth] = dj[h_t>=minidepth] - vdot2(dj, dc)*dc[h_t>=minidepth]   #project avant descente 
    return dj

def apply_post_sand_conservation(psi, sand0, x, h_t, minidepth):
    psi[h_t>=minidepth] *= sand0/trapz(psi,x)
    return psi

def cut_boundary_minidepth(dj,h_t,minidepth, ntot):
    dj[0]=0                  #bords a zero
    dj[ntot-1]=0
    dj[h_t<=minidepth]=0        #peu d'eau : pas de deplacement
    return dj

def rogner(filename):
    
    img = cv2.imread(str(filename)+".png")
    white = True
    save_i_bas=[]
    save_i_haut=[]
    save_j_bas=[]
    save_j_haut=[]
    for i in range(len(img[:,0])):
        if white:
            for j in range (len(img[0,:])):
                if (img[i,j]!=[255,255,255]).all():
                    save_i_bas.append(i)
                    white = False
                    break
    white = True
    for i in range(len(img[:,0])-1,-1,-1):
        if white:
            for j in range (len(img[0,:])):
                if (img[i,j]!=[255,255,255]).all():
                    save_i_haut.append(i)
                    white = False
                    break
    white = True
    for i in range(len(img[0,:])):
        if white:
            for j in range (len(img[:,0])):
                if (img[j,i]!=[255,255,255]).all():
                    save_j_bas.append(i)
                    white = False
                    break
    white = True
    for i in range(len(img[0,:])-1,-1,-1):
        if white:
            for j in range (len(img[:,0])):
                if (img[j,i]!=[255,255,255]).all():
                    save_j_haut.append(i)
                    white = False
                    break
    cv2.imwrite(str(filename)+".png", cv2.imread(str(filename)+".png")[save_i_bas[0]-1:save_i_haut[0]+1,save_j_bas[0]-1:save_j_haut[0]+1])        
    return

def npdf(x):
	return 1/sqrt(2*pi)*exp(-x**2/2)

def cdf(x) : 
	return 0.5*(1+erf(x/sqrt(2)))
     
def pdf(x, alpha) :
	return 2*npdf(x)*cdf(alpha*x)

def create_forcing(t, loc, scale, hmax, alpha = 0) : 
	r = 2/scale*npdf((t-loc)/scale)*cdf(alpha*(t-loc)/scale)
	M = max(r)+0.0000001
	return hmax/M*r + 0.1

def vdot2(dj,dc):
     return sum([dj[i]*dj[i] for i in range(len(dj))])


def Dif_finite(psi, H1):
     Hp = interpolate.interp1d(psi, H1)
     epsdf=1.e-4
     H_psi_dif = zeros(len(H1))
     H_psi_dif[-1]= (Hp(psi[-1])-Hp(psi[-1]-epsdf))
     H_psi_dif[0]=(Hp(psi[0]+epsdf)-Hp(psi[0]))
     for i in range(1,len(psi)-1):
          if psi[i]+epsdf>max(psi):
              H_psi_dif[i] = (Hp(psi[i])-Hp(psi[i]-epsdf))
          elif psi[i]-epsdf<min(psi):
              H_psi_dif[i] = (Hp(psi[i]+epsdf)-Hp(psi[i]))
          else:
              H_psi_dif[i] = (Hp(psi[i]+epsdf)-Hp(psi[i]-epsdf))/2
     return H_psi_dif/epsdf


def sym_calc(f,x):
     return  array([f.subs(t,var) for var in x])

def Dif_Hadamard_NN (psi, H):

    n=np.size(H)
    dx=1
    H_psi_num=np.zeros(n)
    epsL=1.e-30
    epsdf=1.e-4
    finite_index = []
    for i in range(1,n-1):    
        dH = (H[i+1]-H[i-1])
        dy = (psi[i+1]-psi[i-1])
        t_x = dx / (dx**2+dy**2)**0.5
        t_y = dy / (dx**2+dy**2)**0.5 
        n_x = - t_y
        n_y =   t_x
    
        dHdx=dH/dx
        if(abs(dy)>epsL):
            dHdy=dH/dy
        elif (abs(dy)<=epsL):
            print(i,dy,"degrade to finite differences calling hydro(psi+eps)")
            finite_index.append(i)
            dHdy = 0     
        else:
            dHdy = NaN
        H_psi_num[i] = dHdx * n_x + dHdy * n_y 

    if finite_index!=[]:
         for i in finite_index:
              dH = (H[i+1]-H[i-1])
              dy = (psi[i+1]-psi[i-1])    
              t_x = dx / (dx**2+dy**2)**0.5
              t_y = dy / (dx**2+dy**2)**0.5 
              n_x = - t_y
              n_y =   t_x
              dHdx=dH/dx
              Df = Dif_finite(psi, H)
              H_psi_num[i] = dHdx * n_x + Df[i] * n_y 
         
    H_psi_num[0] = 2*H_psi_num[1] - H_psi_num[2]
    H_psi_num[n-1] = 2*H_psi_num[n-2] - H_psi_num[n-3]
    
    return lim_pente_1(H_psi_num)
    return H_psi_num


def lim_pente_1(y):
     n = len(y)
     ener=1
     for dx in range(1,20):
         y0=y.copy()
         for i in range(dx,n-dx):
             ymin=min(y0[i-dx],y0[i+dx])
             ymax=max(y0[i-dx],y0[i+dx])
             y[i]=max(min(y0[i],ymax),ymin)
     return y

def lim_pente_2(y):
    nx =len(y)
    y0 = zeros(nx)
    ener=1
    for dx in range(1,10):
        y0[0:nx]=y[0:nx]   
        for i in range(dx,nx-dx):
            ymin=min(y0[i-dx],y0[i+dx])
            ymax=max(y0[i-dx],y0[i+dx])
            y[i]=max(min(y0[i],ymax),ymin)
        ener0=ener            
        ener= sqrt((y-y0).dot(     y-y0))
        if(dx==1):
            e0=ener
        ener/=e0
        if(dx>1 and ener>ener0):
            y=y0
            break
    return y










# --------- Functions for hydrodynamic model

def k_guo(h,sigma,g):
     hh=h.copy()
     hh[h<0] = NaN
     x = hh * sigma / sqrt( g*hh )
     beta = 2.4901
     k = x**2/hh * (1 - exp(-x**beta))**(-1/beta)
     k[isnan(k)] = 0
     return k

def Ks(h,k):
     X = k*h
     U = tanh(X)*(1+2*X/sinh(2*X))
     kpsi = k**2/(cosh(X)*sinh(X)+X)
     Xpsi = kpsi*h-k
     Upsi = Xpsi* (2*cosh(X)**2-X*sinh(2*X))/cosh(X)**4
     return (U)**-0.5,-0.5*U**(-1.5)*Upsi

def kkpsi(h,k):
     X = k*h
     kpsi = k**2/(cosh(X)*sinh(X)+X)
     return k,kpsi


def norme(X):
     return sum(X**2)**0.5

def normalise (X):
     N = norme(X)
     if N>1e-20:
          return X/N
     else:
          return X*0

def calc_H(H0,k,h,gamma):
     KS,KSpsi = Ks(h,k)
     H1 = H0 * KS
     Hpsi = H0 * KSpsi
     for i in range(len(h)):
          if H1[i]/h[i]>gamma:
               H1[i] = gamma * h[i]
               Hpsi[i] = -gamma
     return H1,Hpsi

def K_shoaling(n2, C, C0) :
    """
        Method which calculates the shoaling coefficient
        :param n: the ratio between the wave celerity and the group celerity
        :param C : the wave celerity
        :param C0 : wave velocity at x = 0 (Forcing)
        :type n: float64
        :type C: float64
        :type C0: float64
        :return: the wave height
            :rtype: float64 
    """

    if C == 0:
        C = 0.0000000000000000001
    return sqrt(1/(2*n2)*C0/C)
     
def shoreline(h_t, x) :
    """
        Method that calculates the shoreline xS, i.e. the point in which the seabed meets the mean water level.
        :param x: The x-axis array
        :param h : The depth of the water : h = h0 - psi 
        :type x: Array of float64
        :type h: Array of float64
        :return: the value of the shoreline and its position in the x array.
    """

    xS_t = -1                                
    nbS_t = -1            
    nbB_t = -1
    xB_t = -1

    xS_t = -1    
    for i in range(len(x)) :                    
        if h_t[i] <= 0 or isnan(h_t[i]):
            xS_t = x[i]                        
            nbS_t = i            
            nbB_t = nbS_t-1
            xB_t = x[nbB_t]
            break
    if     xS_t == -1 :
        print("PB : No intersection between seabed and h0")
        #sys.exit()
        #break
    
    return xS_t, xB_t, nbS_t , nbB_t

def slopecst(psi, Mslope, mobility, x) :
    """
        Method that applies the slope contraint to the given seabed. If the slope of the seabed is steeper than allowed, the seabed is corrected.
    """
    i0 = 0
    
    psi_forward = copy(psi)
    psi_backward = copy(psi)

    for i in range(i0+1, len(x)) :
        psi_forward[i] = min(max(-Mslope[i]*(x[i] - x[i-1]) + psi_forward[i-1], psi[i]) , Mslope[i]*(x[i] - x[i-1]) + psi_forward[i-1])

    for i in range(len(x)-2, i0, -1) : 
        psi_backward[i] = min(max(-Mslope[i]*(x[i+1] - x[i]) + psi_backward[i+1], psi[i]), Mslope[i]*(x[i+1] - x[i]) + psi_backward[i+1])


    return (psi_backward + psi_forward)/2



def shoaling_window_LC(nbS_t, nbB_t, h_t, sigma0, C0, theta0, current_H0, x, gamma, Nwin, dwin) :
    """
    Function which calculates the shoaling variables accoss the cross shore profile, such as C, Cg, k, and H. This function is used by the incremental hydrodynamic model.
    """
    
    first = True
    first2 = True

    length = len(h_t)        
    k_t = zeros(length)
    C_t = zeros(length)
    Cg_t = zeros(length)
    H_t =  zeros(length)
    H_S =  zeros(length)         # H shoaling
    nbB_t = nbS_t-1
    xB_t = x[nbB_t]
    

    Hw = zeros(length)

    tab_xB_start = []
    tab_xB_end = [x[0]]        # cannot start with breaking ...

    tab_nbB_start = []
    tab_nbB_end = [0]
    
    oldC0 = C0
    oldH0 = current_H0

    w = zeros(Nwin)
    k_t= k_guo(h,sigma0,g)
    for i in range(Nwin):
        w[Nwin-i-1] = exp(log(0.01)/(Nwin-1)**2*i**2)

    ibreak = Nwin

    tab_ibreak = ones(len(x))*Nwin

    lcp = ones(len(x))            # linear combinaison parameter
    ii = 0
    for i in range(Nwin) :
        lcp[i] = x[i]/dwin

    Ks_0 = 1

    for i in range(nbS_t) :
        if i == 0 :
            H0 = current_H0
        else : 
            newHi = 0
            newCi = 0
            sumw = 0
            
            for j in range(min(i, Nwin, ibreak)):
                newHi = newHi + w[min(i, Nwin,ibreak)-j-1]*H_t[i-j-1]
                newCi = newCi + w[min(i, Nwin, ibreak)-j-1]*C_t[i-j-1]
                sumw = sumw +  w[min(i, Nwin, ibreak)-j-1]
            
            H0 = newHi/sumw
            C0 = newCi/sumw

        Hw[i] = H0    
        

        #k_t[i] = wvnum(h_t[i], 0, sigma0,i)
        C_t[i] = C0*tanh(k_t[i]*h_t[i])**(1/2)
        
        Cg_t[i] = 0.5*C_t[i]*(1+2*k_t[i]*h_t[i]/sinh(2*k_t[i]*h_t[i]))
        n_t = Cg_t[i]/C_t[i]

        H_S[i] = oldH0 * K_shoaling(n_t, C_t[i], C0)
        if i == 0 :
            Ks_0 = K_shoaling(n_t, C_t[i], C0)

         # linear combinaison of H0 and Hw when w < dwin, otherwise Hw
        H_t[i] = ((1-lcp[i])*oldH0 + lcp[i]*Hw[i] )* K_shoaling(n_t, C_t[i], C0)


        if H_t[i]/h_t[i] >= gamma :            
            if first == True and i > 1 :
                tab_xB_start = append(tab_xB_start, x[i-1])
                tab_nbB_start = append(tab_nbB_start, i-1)
                first = False
            H_t[i] = gamma*h_t[i]
            first2 = True

            ibreak = 1

        else :
            first = True

            if first2 == True and i >1:
                tab_xB_end = append(tab_xB_end, x[i-1])
                tab_nbB_end = append(tab_nbB_end, i-1)
                first2 = False    
        
        tab_ibreak[i] = ibreak
        ibreak = ibreak+1
        

    if len(tab_xB_start) != 0 :
        xB_t = tab_xB_start[0]
        nbB_t = int(tab_nbB_start[0])
    else :
        nbB_t = int(nbS_t-1)
        xB_t = x[nbB_t]
        tab_xB_start = append(tab_xB_start, x[nbB_t])
        tab_nbB_start = append(tab_nbB_start, nbB_t)

    tab_xB_end = append(tab_xB_end, x[nbS_t])
    tab_nbB_end = append(tab_nbB_end, nbS_t)

    H_t = H_t/Ks_0            # to ensure deep water conditions at x = 0
    #H_t += array([(random.random()-0.5)*0.5 for i in H_t])

    return k_t, C_t, Cg_t, H_t, xB_t, nbB_t, H_S, tab_xB_start, tab_nbB_start, tab_xB_end, tab_nbB_end, Hw, tab_ibreak



def sandstock(psi, psi0, x) :
    """
        Method that returns the current sand-stock, i.e. (quantity of sand) / (original quantity of sand).
    """
    initial_stock = trapz(psi0,x)
    current_stock = trapz(psi,x)

            
    return current_stock/initial_stock



def shoaling_window_RONAN(nbS_t, nbB_t, h_t, sigma0, C0, theta0, current_H0, x, gamma, Nwin, g) :
    """
    Function which calculates the shoaling variables accoss the cross shore profile, such as C, Cg, k, and H. This function is used by the incremental hydrodynamic model.
    """
    
    first = True
    first2 = True

    length = len(h_t)        
    k_t = zeros(length)
    C_t = zeros(length)
    Cg_t = zeros(length)
    H_t =  zeros(length)
    H_S =  zeros(length)         # H shoaling
    nbB_t = nbS_t-1
    xB_t = x[nbB_t]
    

    Hw = zeros(length)

    tab_xB_start = []
    tab_xB_end = []

    tab_nbB_start = []
    tab_nbB_end = []
    
    oldC0 = C0
    oldH0 = current_H0
    Nwin = 1

    w = zeros(Nwin)

    k_t = k_guo(h_t,sigma0,g)
    KS = Ks(h_t,k_t)[0]
    fact = 1/cosh(k_t*h_t)
    saveK=[]
    eps = 1e-12
    #Nwin =10
    for i in range(Nwin):
        if Nwin == 1:
            w[0] = 1
            break
        else:
            w[Nwin-i-1] = exp(log(0.01)/(Nwin-1)**2*i**2)
        #w*=0+1

    ibreak = Nwin
    #exit()

    tab_ibreak = ones(len(x))*Nwin


    for i in range(nbS_t) :
        if i == 0 :
            H0 = current_H0
        else : 
            newHi = 0
            newCi = 0
            sumw = 0
            
            for j in range(min(i, Nwin, ibreak)):
                #print(min(Nwin,ibreak)-j-1,i-j-1,H_t[i-j-1],KS[i])
                newHi = newHi + w[min(i, Nwin, ibreak)-j-1]*H_t[i-j-1]
                newCi = newCi + w[min(i,Nwin, ibreak)-j-1]*C_t[i-j-1]
                sumw = sumw +  w[min(i,Nwin, ibreak)-j-1]
            
            H0 = newHi/sumw
            C0 = newCi/sumw

        Hw[i] = H0
        #print(Hw[i])
        


        C_t[i] = C0*tanh(k_t[i]*h_t[i])
        Cg_t[i] = 0.5*C_t[i]*(1+2*k_t[i]*h_t[i]/sinh(2*k_t[i]*h_t[i]))
        n_t = Cg_t[i]/C_t[i]

        H_t[i] = H0 * K_shoaling(n_t, C_t[i], C0)
        if i>Nwin:
            H_t[i] = H0 * KS[i]/KS[i-Nwin]
        else:
            H_t[i] = H0 * 1
        #print(KS[i]/KS[i-1])
        saveK.append(K_shoaling(n_t, C_t[i], C0) )
        H_S[i] = oldH0 * K_shoaling(n_t, C_t[i], C0)



        if H_t[i]/h_t[i] >= gamma :            
            if first == True and i !=0 :
                tab_xB_start = append(tab_xB_start, x[i-1])
                tab_nbB_start = append(tab_nbB_start, i-1)
                first = False
            H_t[i] = gamma*h_t[i]
            first2 = True

            ibreak = 1

        else :
            first = True

            if first2 == True and i !=0:
                tab_xB_end = append(tab_xB_end, x[i-1])
                tab_nbB_end = append(tab_nbB_end, i-1)
                first2 = False    
        
        tab_ibreak[i] = ibreak
        ibreak = ibreak+1
        

    if len(tab_xB_start) != 0 :
        xB_t = tab_xB_start[0]
        nbB_t = int(tab_nbB_start[0])
    else :
        nbB_t = int(nbS_t-1)
        xB_t = x[nbB_t]
        tab_xB_start = append(tab_xB_start, x[nbB_t])
        tab_nbB_start = append(tab_nbB_start, nbB_t)

    tab_xB_end = append(tab_xB_end, x[nbS_t])
    tab_nbB_end = append(tab_nbB_end, nbS_t)

    curve_coef = -1
    
    f1 = lambda x: x**1*(1+curve_coef*(x-1))
    #f = lambda x: exp(curve_coef*x)
    f2 = lambda x: x**(0.2)
    f3 = lambda x: 0.5+0.5*sin(x*pi+1.5*pi)
    f4 = lambda x: 0.5*tanh(30*x-1.4*pi/2)+0.5
    def f5(x):
        if x<0.09:
            return exp(50*x-5)
        else:
            return -exp(-x*50+3.5)+1
    if len(tab_nbB_end)>=len(tab_nbB_start)+1 :
        #print(tab_nbB_start,tab_nbB_end)
        for i in range(len(tab_xB_start)):
            istart = int(tab_nbB_start[i])+ 1
            istop = int(tab_nbB_end[i+1]) +1
            if (istart>=0 and istart!=istop and istop>=0 ):
                    #print(istart,istop,"wtf")

                    #istop = int(min(istop,len(x)-1))
                    
                    
                    pre_deferlement =int( (istop- istart)*0.1 ) #10%
                    if istart-pre_deferlement>=0:
                            istart -= pre_deferlement
                    Ha = H_t[istart]
                    Hb = H_t[istop]
                    htot = h_t[istart:istop]
                    
                #print(htot)
                    hmin = min(htot)
                    hmax = max(htot)
                    if abs(hmin - hmax)>eps and abs(istop-istart)>eps:
                        #print ( hmax - hmin)
                        for ii in range(istart,istop):
                            X = (ii-istart)/(istop-istart)
                            Y = (hmax - h_t[ii] ) / ( hmax - hmin) # in 0 1
                            H_t[ii] = Ha + (Hb - Ha) * f3(X) * f5(Y) #3 5 pas mal 2/5

    """

    figure(1)

    plot(saveK,"b")
    plot(KS,"r")
    title("KS")
    show()
    """

    return k_t, C_t, Cg_t, H_t, xB_t, nbB_t, H_S, tab_xB_start, tab_nbB_start, tab_xB_end, tab_nbB_end, Hw, tab_ibreak



def vec_n(X,Y,psi,i,j):
    n_i, n_j = shape(psi)
    P0 = array([X[i , j], Y[i , j], psi[i , j]])    
    if i<n_i-1 and j<n_j-1:
        P1 = array([X[i+1 , j], Y[i+1 , j], psi[i+1 , j]])
        P2 = - array([X[i , j+1], Y[i , j+1], psi[i, j+1]])
        t1 = P1 - P0
        t2 = P2 - P0

    else:
        P1 = array([X[i-1 , j], Y[i-1 , j], psi[i-1 , j]])
        P2 = - array([X[i , j-1], Y[i , j-1], psi[i, j-1]])
        t1 = - (P1 - P0)
        t2 = (P2 - P0)
        if j==n_y-1 and i==0:
            t2 = - (P2 - P0)

    t1, t2 = normalise(t1), normalise(t2)
    return normalise(cross(t1,t2))
    
    
def limiteur_2D(z):
     n_i, n_j = shape(z)
     z_lim = zeros(( n_i, n_j ))
     for i in range(n_i):
          z_lim [i, :] = lim_pente_1(z[i, :])
     return z_lim


def Dif_Hadamard_2D (psi, H, X , Y):
    epsL=1.e-30
    epsdf=1.e-8
    dx = X[0,1] - X[0,0]
    dy = Y[1,0] - Y[0,0]
    finite_index = []
    #print(dx, dy)
    n_i, n_j = shape(psi)
    H_psi_num = zeros((n_i,n_j))
    for i in range(1, n_i-1):
        for j in range(1, n_j-1):
            n_x, n_y, n_z = vec_n(X,Y,psi,i,j)
            dHx =(H[i+1,j]-H[i-1,j])  #calculate dHx
            dHy =(H[i,j+1]-H[i,j-1]) #calculate dHy
            
            dzx = (psi[i+1,j]-psi[i-1,j]) # calculate dzx
            dzy = (psi[i,j+1]-psi[i,j-1]) # calculate dz

            dHdx = dHx / dx
            
            dHdy = dHy / dy
            #print(dHdx,dHdy)
            
            dHdz = dHx/dzx + dHy/dzy
            #print("grad",dHdx,dHdy,dHdz)
            #print("n",n_x , n_y  , n_z)
            
            
            if(abs(dzx)>epsL and abs(dzy)>epsL):
                dHdz = (dHx/dzx + dHy/dzy)/2
            elif (abs(dzx)<=epsL or abs(dzy)<=epsL) :
                print(i,dz,"degrade to finite differences calling hydro(psi+eps)")
                finite_index.append(i)
                dHdz = 0
            else:
                dHdz = NaN

            H_psi_num[i,j] = dHdx * n_x + dHdy * n_y  + dHdz * n_z

    if finite_index!=[]:
        for i in finite_index:
            print('banane')
            
    H_psi_num[0,:] = 2*H_psi_num[1,:] - H_psi_num[2,:]
    H_psi_num[:,0] = 2*H_psi_num[:,1] - H_psi_num[:,2]

    H_psi_num[n_i-1,:] = 2*H_psi_num[n_i-2,:] - H_psi_num[n_i-3,:]
    H_psi_num[:,n_j-1] = 2*H_psi_num[:,n_j-2] - H_psi_num[:,n_j-3]

    return limiteur_2D(H_psi_num)


def H_shoaling_2D(psi, h0, x, sigma0, C0, H0, gamma, Nwin, g):
     n_i, n_j = shape(psi)
     H = zeros((n_i,n_j))
     K = zeros((n_i,n_j))
     for i in range(n_i):
          psi_x = psi[i,:]
          h = h0 - psi_x
          
          xS_t, xB_t, nbS_t, nbB_t = shoreline(h, x)
          k_t, C_t, Cg_t, H_t, xB_t, nbB_t, H_S, tab_xB_start, tab_nbB_start, tab_xB_end, tab_nbB_end, Hw, tab_ibreak = shoaling_window_RONAN(nbS_t, nbB_t, h, sigma0, C0, 0, H0, x, gamma, Nwin, g)
          H_t[H_t<=1e-6] = NaN
          H[i,:] = H_t
          K [i,:] = k_t
     return H, K

def func_to_plot_2D(psi, H, X, Y, h0, itera, niter, ax):
     ax.clear()
     ratio= [1*3,40,30*3]
     ax.set_box_aspect((ptp(X)*ratio[0], ptp(Y)*ratio[1], ptp(psi)*ratio[2]))

     #ax.imshow(psi, cmap=cm.get_cmap("viridis").reversed(), origin='lower', extent=[0,L_x,0,L_y],interpolation='spline16', aspect='auto')

     
     ax.plot_surface(X, Y, psi-h0, cmap=cm.get_cmap("copper"), linewidth=0, rstride=1, cstride=1, alpha=0.8)
     ax.plot_surface(X, Y, H , cmap=cm.get_cmap('winter'), linewidth=0, rstride=1, cstride=1, alpha=0.8)
     ax.set_xlabel('x',fontsize = 15)
     ax.set_ylabel('y',fontsize = 15)
     #ax.set_title('Wave height $H$ and Sea Bottom Profile $\psi$ \n at iteration '+str(itera+1)+' / '+str(niter), weight="bold",fontsize = 20)
     ax.set_zlabel('$\psi$')
     ax.view_init(5, 270)

def func_to_plot_2D_flat(psi, H, X, Y, h0, L_x, L_y, itera, niter, axs):



     axs[0].clear()
     axs[1].clear()
     
     labsize = 15
     ttlsize = 20
     
     im0 = axs[0].imshow(H, cmap=cm.get_cmap('winter'), origin='lower', extent=[0,L_x,0,L_y],interpolation='spline16', aspect='auto')
     im1 = axs[1].imshow(psi-h0, cmap=cm.get_cmap("copper"), origin='lower', extent=[0,L_x,0,L_y],interpolation='spline16', aspect='auto')
     
     cb1 = colorbar(im0, ax=axs[0])
     cb2 = colorbar(im1, ax=axs[1])
     
     axs[1].set_xlabel('x',fontsize = labsize)
     axs[1].set_ylabel('y',fontsize = labsize)

     axs[0].set_xlabel('x',fontsize = labsize)
     axs[0].set_ylabel('y',fontsize = labsize)
     
     #axs[0].set_title('Wave height $H$', weight="bold",fontsize = ttlsize)
     #axs[1].set_title('Sea Bottom Profile $\psi$', weight="bold",fontsize = ttlsize)
     return cb1, cb2
     
     
     
     
def trapz_2D(x, y, z):
     n_i, n_j = shape(z)
     I = 0
     Z_x =[]
     for i in range(n_i):
          Z_x. append(trapz(z[ i , : ], x))
     return trapz ( Z_x, y )
