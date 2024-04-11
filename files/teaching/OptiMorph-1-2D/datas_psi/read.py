from matplotlib.pylab import *
from scipy import interpolate

initial = loadtxt('initial.txt', delimiter=';')
final = loadtxt('final.txt', delimiter=';')
Xi, Yi = initial[:,0], initial[:,1]
Xf, Yf = final[:,0], final[:,1]


fi = interpolate.interp1d(Xi,Yi)
ff = interpolate.interp1d(Xf,Yf)

dx = 0.1
print(int(Xi[0]),int(Xi[-1]))
new_x = arange(int(Xi[0]), int(Xi[-1]), dx)

Zi = fi( new_x )
Zf = ff( new_x )


h0i = abs(min(Zi))
h0f = abs(min(Zf))

print('---------')
print("h0 initial ", h0i)
print("h0 final ", h0f)
print('---------')


new_x -= min(new_x)

x0 = 20
ind = abs(new_x - x0)<1e-3
ind = where(ind==True)[0][0]
y0 = Zi[ind]

x1 = new_x[-1]
y1 =Zi[-1]

a = (y0 - y1) / (x0 - x1)
b = (y0 * x1 - y1 * x0) / (x1 - x0)

def f_lin(a, b, x):
    return a*x + b


Zi[:] = f_lin(a, b, new_x[:])

print(abs(min(Zi)))
print(Zi[x0])
savetxt("SANDS_initial.dat",Zi[ind:])
savetxt("SANDS_final.dat",Zf[ind:])


plot(new_x, Zi + h0i)
plot(new_x, Zf + h0f)
print(new_x[-1])
show()
