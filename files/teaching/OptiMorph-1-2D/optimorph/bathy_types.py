from pathlib import Path
from matplotlib.pylab import *

def load_bathy(bathy_type, x, h0, nwater, ntot, seuil, psi_datas):
    psi_datas = Path(psi_datas)
    if not psi_datas.exists():
        raise OSError(f"file not found {psi_datas}")
        
    psif = ones(ntot) * NaN
    HRMS = ones(ntot) * NaN
    if bathy_type==0:
        psi = (x*h0/nwater)+seuil
        h = h0 - psi
        plot(psi)

    elif bathy_type==1:
        psi_lim = (x*h0/nwater)+seuil
        last_v = psi_lim[-1]
        psi = (last_v/ 0.99)*(0.99*(x/x[-1])**2) #concave
        h = h0 - psi
    elif bathy_type==2:     
        psi_lim = (x*h0/nwater)+seuil
        last_v = psi_lim[-1]
        psi= 2*psi_lim-(last_v/ 0.99)*(0.99*(x/x[-1])**2) #convexe
        h = h0 - psi
    elif bathy_type==3:     
        psi_lim=(x*h0/nwater)+seuil
        xgauss=5
        bosse_length=100 #épaisseur de la bosse
        bosse_size=1.5*1.1 #taille de la bosse
        bosse_profondeur=3 #profondeur de la bosse

        X=[]
        X=np.linspace(-xgauss,xgauss,bosse_length)
        psi0_f=pdf(X, -1)
        M=max(psi0_f)
        ind=np.where(psi0_f==M)[0][0]
        bosse_xmidle=np.where(psi_lim>(h0-bosse_profondeur))[0][0]

        bosse_xstart=bosse_xmidle-bosse_length//2 #position  x du commencement de la bosse
        bosse_hauteur=h0-bosse_profondeur+psi_lim[bosse_xmidle]
        bosse_xmidle=ind+bosse_xstart
        psi0_f=psi0_f*bosse_size/M
        psi0_f2=[]
        psi0_f2=[0 for i in range(bosse_xstart)]+list(psi0_f)
        psi0_f2 +=[0 for i in range(ntot-len(psi0_f2))]
        psi = array(psi0_f2+psi_lim)
        h = h0 - psi
    elif bathy_type==4:
        psi = loadtxt(psi_datas / 'psi_LIP1B_initial.dat')
        h0 = 4.1
        h = h0 - psi
        #nwater = int(param[3])
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        psif =  loadtxt(psi_datas / 'psi_LIP1B_final.dat')
        HRMS = loadtxt(psi_datas / 'HRMS_LIP1B.dat')
    elif bathy_type==5:     
        psi = loadtxt(psi_datas / 'psi_T06_lineaire.dat')
        h0 = 10
        h = h0 - psi
        ntot = len(psi)
        nwater = ntot - 10
        x = linspace(0,ntot-1,ntot)
    elif bathy_type==6:     
        psi = loadtxt(psi_datas / 'psi_T06_concave.dat')
        h0 = 10
        h = h0 - psi
        ntot = len(psi)
        nwater = ntot - 10
        x = linspace(0,ntot-1,ntot)
    elif bathy_type==7:     
        psi = loadtxt(psi_datas / 'psi_T06_convexe.dat')
        h0 = 10
        h = h0 - psi
        ntot = len(psi)
        nwater = ntot - 10
        x = linspace(0,ntot-1,ntot)
    elif bathy_type==8:     
        psi = loadtxt(psi_datas / 'psi_T016_lineaire.dat')
        h0 = 10
        h = h0 - psi
        ntot = len(psi)
        nwater = ntot - 10
        x = linspace(0,ntot-1,ntot)
    elif bathy_type==9:     
        psi = loadtxt(psi_datas / 'psi_T016_concave.dat')
        h0 = 10
        h = h0 - psi
        ntot = len(psi)
        nwater = ntot - 10
        x = linspace(0,ntot-1,ntot)
    elif bathy_type==10:     
        psi = loadtxt(psi_datas / 'psi_T016_convexe.dat')
        h0 = 10
        h = h0 - psi
        ntot = len(psi)
        nwater = ntot - 10
        x = linspace(0,ntot-1,ntot)
    elif bathy_type==11:     
        psi = loadtxt(psi_datas / 'psi_LIP1C_initial.dat')
        h0 = 4.1
        h = h0 - psi
        #H0 = ones(len(time))*0.6
        #nwater = int(param[3])
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        psif =  loadtxt(psi_datas / 'psi_LIP1C_final.dat') - h0
        HRMS = loadtxt(psi_datas / 'HRMS_LIP1C.dat')

    elif bathy_type==12:     
        psi = loadtxt(psi_datas / 'psi_LIP1C_final_interp.dat')
        h0 = 4.1
        h = h0 - psi
        #H0 = ones(len(time))*0.6
        #nwater = int(param[3])
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        psif =  loadtxt(psi_datas / 'psi_LIP1B_final.dat')

    elif bathy_type==13:     
        psi = loadtxt(psi_datas / 'psi_LIP1B_initial_interp.dat')
        h0 = 4.1
        h = h0 - psi
        #H0 = ones(len(time))*0.6
        #nwater = int(param[3])
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        psif =  loadtxt(psi_datas / 'psi_LIP1C_final.dat')

    elif bathy_type==14:     
        psi = loadtxt(psi_datas / 'psi_LIP1B_final.dat')
        h0 = 4.1
        h = h0 - psi
        #H0 = ones(len(time))*0.6
        #nwater = int(param[3])
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        HRMS = loadtxt(psi_datas / 'HRMS_LIP1B.dat')
        
    elif bathy_type==15:     
        psi = loadtxt(psi_datas / 'psi_LIP1C_final.dat')
        h0 = 4.1
        h = h0 - psi
        #H0 = ones(len(time))*0.6
        #nwater = int(param[3])
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        HRMS = loadtxt(psi_datas / 'HRMS_LIP1C.dat')

    elif bathy_type==16:     
        psi = loadtxt(psi_datas / 'psi_LIP1B_initial_interp_d1.dat')
        h0 = 4.1
        h = h0 - psi
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        psif =  loadtxt(psi_datas / 'psi_LIP1B_final.dat')
        
    elif bathy_type==17:     
        psi = loadtxt(psi_datas / 'psi_LIP1C_initial_interp_d1.dat')
        h0 = 4.1
        h = h0 - psi
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        psif =  loadtxt(psi_datas / 'psi_LIP1C_final.dat')

    elif bathy_type==18:     
        psi = loadtxt(psi_datas / 'bathy_copterI.dat')
        h0 = 0.55
        h = h0 - psi
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        H0 = loadtxt(psi_datas / 'H(t,0)_xb_2.dat')
        return x, h0, ntot, psi, psif, H0, HRMS
    elif bathy_type==19:     
        psi = loadtxt(psi_datas / 'kevin_i.dat')
        h0 = 3.91
        h = h0 - psi
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        psif =  loadtxt(psi_datas / 'kevin_f.dat')
        H0 = 0.53
        return x, h0, ntot, psi, psif, H0, HRMS
    elif bathy_type==20:     
        psi = loadtxt(psi_datas / 'SANDS_initial.dat')
        h0 = abs(min(psi))
        psi += h0
        h = h0 - psi
        ntot = len(psi)
        x = linspace(0,ntot-1,ntot)
        dx = 0.1
        x = arange(0, ntot * dx, 0.1)
        psif =  loadtxt(psi_datas / 'SANDS_final.dat')
        H0 = 0.53
        return x, h0, ntot, psi, psif, H0, HRMS
    else:
        raise ValueError("bathy_type {bathy_type} shall be in [0, 20]")
    return x, h0, ntot, psi, psif, False, HRMS 






def is_point_inside_shape(x, y, shape):
    """
    Check if a point (x, y) is inside a 2D shape defined by a list of points.

    Parameters:
    - x, y: Coordinates of the point to be checked.
    - shape: List of tuples representing the vertices of the shape.

    Returns:
    - True if the point is inside the shape, False otherwise.
    """
    num_vertices = len(shape[:,0])

    # Initialize counter for intersections
    intersections = 0

    # Iterate through each pair of consecutive vertices
    for i in range(num_vertices):
        x1, y1 = shape[i]
        x2, y2 = shape[(i + 1) % num_vertices]

        # Check if the ray from the point intersects the edge
        if ((y1 <= y < y2) or (y2 <= y < y1)) and (x < (x2 - x1) * (y - y1) / (y2 - y1) + x1):
            intersections += 1

    # If the number of intersections is odd, the point is inside the shape
    return intersections % 2 == 1

def struct_2d(mr,nr,x,y,Pts_str,h_str):
    Mat_str = np.zeros((nr,mr))
    for i in range(nr):
        for j in range(mr):
            #print(is_point_inside_shape(x[i,j],y[i,j],Pts_str),"coord?")
            Mat_str[i,j] = h_str * is_point_inside_shape(x[i,j],y[i,j],Pts_str)
        
 
    return Mat_str


def load_bathy_2D(bathy_type, x, y, h0, nwater, seuil, psi_datas):
    psi_datas = Path(psi_datas)
    if not psi_datas.exists():
        raise OSError(f"file not found {psi_datas}")

    if bathy_type==0:
        X, Y = np.meshgrid(x,y)
        
        psi = 0.01*X+0.01*Y
        h = h0 - psi
        n_i = len(x)
        n_j = len(y)

    elif bathy_type==1:
        X, Y = np.meshgrid(x,y)
    
        psi = (X*h0/nwater)+seuil + 0.001 * Y
        #imshow(psi, cmap=cm.get_cmap('winter'), origin='lower',interpolation='spline16', aspect='auto')
        #show()
        
        #show()
        h = h0 - psi

        n_i = len(x)
        n_j = len(y)

    elif bathy_type==2:
        X, Y = np.meshgrid(x,y)
        n_i = len(x)
        n_j = len(y)
        
        psi = (X*h0/nwater)+seuil + 0.0001 * Y
        h = h0 - psi
        
        # Geometric parameters of the structure
        h_str = 1.5 # could be vectorial
        thick_str = 4
        alpha_str = 45
        beta_str = 45
        L1_str = 10
        L2_str = 10
        l1_str = 15
        l2_str = 15
        ystart = 20

        nb_pts = 18
        Pts_str = np.zeros((nb_pts,2))

        Pts_str[0] = [0,0]
        Pts_str[1] = [-L1_str,0]
        Pts_str[2] = Pts_str[1] - l1_str * array([cos(radians(alpha_str)),sin(radians(alpha_str))])
        Pts_str[3] = Pts_str[2] - [thick_str,0]
        Pts_str[4] = Pts_str[1] - [thick_str,0]
        Pts_str[5] = Pts_str[4] - [L2_str,0]
        Pts_str[6] = Pts_str[5] - l2_str * array([cos(radians(beta_str)),sin(radians(beta_str))])
        Pts_str[7] = Pts_str[6] - [thick_str,0]
        Pts_str[8] = Pts_str[5] - [thick_str,0]
        Pts_str[9] = Pts_str[8] + [0,thick_str]
        Pts_str[10] = Pts_str[9] + l2_str*array([-cos(radians(beta_str)),sin(radians(beta_str))])
        Pts_str[11] = Pts_str[10] + [thick_str,0]
        Pts_str[12] = Pts_str[9] + [thick_str,0]
        Pts_str[13] = Pts_str[12] + [L2_str,0]
        Pts_str[14] = Pts_str[13] + l1_str*array([-cos(radians(alpha_str)),sin(radians(alpha_str))])
        Pts_str[15] = Pts_str[14] + [thick_str,0]
        Pts_str[16] = Pts_str[13] + [thick_str,0]
        Pts_str[17] = Pts_str[16] + [L1_str,0]
        
        # shift the structure
        Pts_str[:,0] += 100 # x shifting
        Pts_str[:,1] += n_j//2 # y shifting
        
        
        additional = struct_2d(n_i,n_j,X,Y,Pts_str,h_str)
        
        psi += additional


    elif bathy_type==3:
        X, Y = np.meshgrid(x,y)
        n_i = len(x)
        n_j = len(y)

        psi = (X*h0/nwater)+seuil + 0.001 * Y
        
        psi[:int(n_j * 0),int(n_i * 0.2):int(n_i * 0.3)] += 5
        psi[int(n_j * 0.7):,int(n_i * 0.2):int(n_i * 0.3)] += 5
        
        
        h = h0 - psi




        

    else:
        raise ValueError("bathy_type {bathy_type} shall be in [[0, 3]]")
    return X, Y, h0, n_i, n_j, psi 

