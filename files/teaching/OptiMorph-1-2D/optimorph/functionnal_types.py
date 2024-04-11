from matplotlib.pylab import *
from optimorph.utils_functions import *


def load_functionnal(psi,H_t,id_cost):
    if(id_cost==1):
       H_psi_num_N = Dif_Hadamard_NN(psi,H_t) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2       
       dj   = H_t*H_psi_num_N    #grad J
    elif(id_cost==2):
       H_psi_num_N = 0.5*Dif_Hadamard_NN(psi,H_t**2) # Numerical solution from Hadamard detailed numerical
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2                         
       dj   = H_psi_num_N      
    elif(id_cost==3):       
       H_psi_num_N1 = Dif_Hadamard_NN(psi,H_t) 
       H_psi_num_N1[isnan(H_psi_num_N1)] = 0 # put NaN to 0
       H_psi_num_N2 = 0.5*Dif_Hadamard_NN(psi,H_t**2) # Numerical solution from Hadamard detailed numerical
       H_psi_num_N2[isnan(H_psi_num_N2)] = 0 # put NaN to 0
       H_psi_num_N=0.5*(H_psi_num_N1+H_t*H_psi_num_N2)
       cost=0.5*norme(H_t)**2
       dj   = H_psi_num_N
    elif(id_cost==4): # U^2
       H_psi_num_N = Dif_Hadamard_NN(psi,u_t**2) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2       
       dj   = H_psi_num_N   #grad J
    elif(id_cost==5): #U^2.H
       H_psi_num_N = Dif_Hadamard_NN(psi,u_t**2*H_t) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2       
       dj   = H_psi_num_N   #grad J
    elif(id_cost==6): #U.H^2
       H_psi_num_N = Dif_Hadamard_NN(psi,u_t*H_t**2) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2       
       dj   = H_psi_num_N   #grad J
    elif(id_cost==7): #U^2.H^2
       H_psi_num_N = Dif_Hadamard_NN(psi,u_t**2*H_t**2) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2
       dj   = H_psi_num_N   #grad J
    elif(id_cost==8): #tau.H^2
       H_psi_num_N = Dif_Hadamard_NN(psi,H_t**2*tau_t) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2
       dj   = H_psi_num_N   #grad J
    elif(id_cost==9): #tau^2.H
       H_psi_num_N = Dif_Hadamard_NN(psi,tau_t**2*H_t) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2
       dj   = H_psi_num_N   #grad J
    elif(id_cost==10): #tau.H^2.u
       H_psi_num_N = Dif_Hadamard_NN(psi,H_t**2*u_t*tau_t) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2
       dj   = H_psi_num_N   #grad J
    if(id_cost==11):
       H_psi_num_N = Dif_Hadamard_NN(psi,H_t)
       Tau_psi_num_N = Dif_Hadamard_NN(psi,tau_t) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       H_psi_num_N[isnan(Tau_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2       
       dj   = 9.81*1000/16*H_t*H_psi_num_N + Tau_psi_num_N    #grad J
    else: #H2 
       H_psi_num_N = Dif_Hadamard_NN(psi,H_t) 
       H_psi_num_N[isnan(H_psi_num_N)] = 0 # put NaN to 0
       cost=0.5*norme(H_t)**2       
       dj   = H_t*H_psi_num_N    #grad J
    return dj, cost
