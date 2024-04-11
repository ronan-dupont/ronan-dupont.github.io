from matplotlib.pylab import *
plot(loadtxt('bathy_copterI.txt'))
show()

"""
Psi_I = loadtxt('psi_LIP1C_initial.dat')[:-15]
Psi_F = loadtxt('psi_LIP1C_final.dat')[:-15]
x = linspace(0,len(Psi_F)-1,len(Psi_F))
plot(Psi_I)
plot(Psi_F)
Ii = trapz(Psi_I,x)
If = trapz(Psi_F,x)
print (Ii,If,"Initial integration","final")
print(If/Ii,"Rapport")
print((Ii-If)/If,"Rapport")
show()

"""
