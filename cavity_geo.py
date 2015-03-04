# Cavity
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib.ticker import MaxNLocator, LinearLocator

# plt.close('all')
## ABCD Matrices
#({{1, 0},{-2/R1, 1}}).({{1, L - z},{0, 1}}).({{1, 0},{-2/R2, 1}}).({{1, z},{0, 1}})
deg = 4
R1nom = 10e-2
R2nom = 25e-2
theta = deg*np.pi/180
IC = .008
rh = 50e-6
fac = 1.01

# X-axis
R1x = R1nom * np.cos(theta)
R2x = R2nom * np.cos(theta)
# Y-axis
R1y = R1nom / np.cos(theta)
R2y = R2nom / np.cos(theta)


fsr = 154.564e6 
cl = 3e8 
lam = 1.070e-6 
L = cl/fsr
# ti = .008
# rh = 50e-6

num = 100
# X-axis
Rx = np.zeros(num)
w0x = np.zeros(num)
w1x = np.zeros(num)
w2x = np.zeros(num)
# Y-axis
Ry = np.zeros(num)
w0y = np.zeros(num)
w1y = np.zeros(num)
w2y = np.zeros(num)

# X-axis
def Ax(z):
    return (-2 * L + R2x + 2 * z) / R2x

def Bx(z):
    return (L * R2x - 2 * L * z + 2 * z ** 2) / R2x

def Cx(z):
    return -2 * (-2 * L + R1x + R2x + 2 * z) / (R1x * R2x)

def Dx(z):
    return (-2 * L * R2x + R1x * R2x + 4 * L * z - 2 * R1x * z - 4 * z ** 2) / (R1x * R2x)

def lim_findx(z):
    return abs(Ax(z) + Dx(z)) - 1.95

# Y-axis
def Ay(z):
    return (-2 * L + R2y + 2 * z) / R2y

def By(z):
    return (L * R2y - 2 * L * z + 2 * z ** 2) / R2y

def Cy(z):
    return -2 * (-2 * L + R1y + R2y + 2 * z) / (R1y * R2y)

def Dy(z):
    return (-2 * L * R2y + R1y * R2y + 4 * L * z - 2 * R1y * z - 4 * z ** 2) / (R1y * R2y)

def lim_findy(z):
    return abs(Ay(z) + Dy(z)) - 1.95

# X-axis
zmaxx = fsolve( lim_findx, R2x+.01)
zminx = fsolve( lim_findx, R1x-.01)
zspanx = np.linspace(zminx, zmaxx, num)

# Y-axis
zmaxy = fsolve( lim_findy, R2y+.01)
zminy = fsolve( lim_findy, R1y-.01)
zspany = np.linspace(zminy, zmaxy, num)

for i in range(num):
    # X
    Rx[i] = 2 * Bx(zspanx[i]) / (Ax(zspanx[i]) - Dx(zspanx[i]))
    w1x[i] = ( 2 * lam * abs(Bx(zspanx[i])) / (np.pi * (4 - (Ax(zspanx[i]) + Dx(zspanx[i])) ** 2) ** .5)) ** .5
    w0x[i] = lam * Rx[i] * w1x[i] / (np.pi **2 * w1x[i] ** 4 + lam **2 * Rx[i] ** 2) ** .5
    qpx = 1 / Rx[i] - 1j * lam / (np.pi * w1x[i] ** 2)
    propx = zspanx[i] - qpx.real / abs(qpx) ** 2
    w2x[i] = w0x[i] * (1 + (lam * propx / np.pi / w0x[i]**2) ** 2) ** .5
    # Y
    Ry[i] = 2 * By(zspany[i]) / (Ay(zspany[i]) - Dy(zspany[i]))
    w1y[i] = ( 2 * lam * abs(By(zspany[i])) / (np.pi * (4 - (Ay(zspany[i]) + Dy(zspany[i])) ** 2) ** .5)) ** .5
    w0y[i] = lam * Ry[i] * w1y[i] / (np.pi **2 * w1y[i] ** 4 + lam **2 * Ry[i] ** 2) ** .5
    qpy = 1 / Ry[i] - 1j * lam / (np.pi * w1y[i] ** 2)
    propy = zspany[i] - qpy.real / abs(qpy) ** 2
    w2y[i] = w0y[i] * (1 + (lam * propy / np.pi / w0y[i]**2) ** 2) ** .5

# plots
fig, ax = plt.subplots()
ax.plot(zspanx*10**2, w0x*10**6, 'k-', label = 'X')
ax.plot(zspany*10**2, w0y*10**6, 'k--', label = 'Y')
ax.set_xlabel('Mirror Separation [cm]')
ax.set_ylabel('$\omega_0$ [$\mu$m]')
ax.set_title('R1 = ' + str(R1x/ np.cos(theta)) + ' cm, R2 =  ' + str(R2x/ np.cos(theta)) + ' cm, deg = ' + str(deg))
ax.yaxis.set_major_locator(LinearLocator(5))
plt.legend()
ax2 = ax.twinx()
ax2.plot(zspanx*10**2, w2x*10**3, 'r-')
ax2.plot(zspany*10**2, w2y*10**3, 'r--')
ax2.set_ylabel('$\omega_{mirror}$ [$\mu$m]', color = 'r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.yaxis.set_major_locator(LinearLocator(5))
ax2.grid()

## Calculation loss and other things
def finesse(loss, ti):
    return 2 * np.pi / (loss + ti)
def buildup(fin, ti):
    return ti * fin**2/np.pi**2
def loss(rh, rspot):
    r = rh/rspot
    return 1 - np.exp(-4 * r ** 2)
def outcouple(rh, rspot, q):
    r = rh/rspot
    return np.exp(-4 * q ** 2 * r ** 2)
# plots
fig, ax = plt.subplots()
ax.plot(zspanx*10**2, finesse(loss(rh,(w2x*w2y)**.5), IC), 'k-',)
ax.set_xlabel('Mirror Separation [cm]')
ax.set_ylabel('Finesse')
ax.set_title('r$_{hole}$ = ' + str(rh*10**6) +  ' $\mu$m \n R1 = ' + str(R1x/ np.cos(theta)) + ' cm, R2 =  ' + str(R2x/ np.cos(theta)) + ' cm, deg = ' + str(deg) + '\n')
ax.yaxis.set_major_locator(LinearLocator(5))
ax.xaxis.set_major_locator(LinearLocator(5))
ax2 = ax.twinx()
ax2.plot(zspanx*10**2, buildup(finesse(loss(rh,(w2x*w2y)**.5), IC), IC), 'r-')
ax2.set_ylabel('Buildup', color = 'r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.yaxis.set_major_locator(LinearLocator(5))
ax2.grid()

fig, ax = plt.subplots()
ax.plot(zspanx*10**2, finesse(loss(rh,(w2x*w2y)**.5), IC), 'k-',)
ax.set_xlabel('Mirror Separation [cm]')
ax.set_ylabel('Finesse')
ax.set_title('r$_{hole}$ = ' + str(rh*10**6) +  ' $\mu$m \n R1 = ' + str(R1x/ np.cos(theta)) + ' cm, R2 =  ' + str(R2x/ np.cos(theta)) + ' cm, deg = ' + str(deg) + '\n')
ax.yaxis.set_major_locator(LinearLocator(5))
ax.xaxis.set_major_locator(LinearLocator(5))


## prop around the cavity.
#start at CM 1 and go to CM2
num = 500
# X-axis
zx = fac*zminx
Rx = 2 * Bx(zx) / (Ax(zx) - Dx(zx))
w1x = ( 2 * lam * abs(Bx(zx)) / (np.pi * (4 - (Ax(zx) + Dx(zx)) ** 2) ** .5)) ** .5
w0x = lam * Rx * w1x / (np.pi **2 * w1x ** 4 + lam **2 * Rx ** 2) ** .5
qpx = (1 / Rx - 1j * lam / (np.pi * w1x ** 2))
propx = zx - qpx.real / abs(qpx) ** 2

z_propx = np.linspace(-1*(zx-propx), propx, num)
w_propx = np.zeros(num)
R_propx = np.zeros(num)

# Y-axis
zy = 1.01*zminx
Ry = 2 * By(zy) / (Ay(zy) - Dy(zy))
w1y = ( 2 * lam * abs(By(zy)) / (np.pi * (4 - (Ay(zy) + Dy(zy)) ** 2) ** .5)) ** .5
w0y = lam * Ry * w1y / (np.pi **2 * w1y ** 4 + lam **2 * Ry ** 2) ** .5
qpy = (1 / Ry - 1j * lam / (np.pi * w1y ** 2))
propy = zy - qpy.real / abs(qpy) ** 2

z_propy = np.linspace(-1*(zy-propy), propy, num)
w_propy = np.zeros(num)
R_propy = np.zeros(num)


for i in range(num):
    # X
    w_propx[i] = w0x * (1 + (lam * z_propx[i]/ np.pi / w0x**2) ** 2) ** .5
    R_propx[i] = z_propx[i] * (1 + (np.pi * w0x**2/lam/z_propx[i])**2)
    # Y
    w_propy[i] = w0y * (1 + (lam * z_propy[i]/ np.pi / w0y**2) ** 2) ** .5
    R_propy[i] = z_propy[i] * (1 + (np.pi * w0y**2/lam/z_propy[i])**2)

# then the rest
# X-axis 
z_totalx = L - zx
qp2x = 1/(1/R_propx[i] -1j * lam / np.pi/w_propx[i]**2)
qp3x = qp2x/(1-2*qp2x/R2x)
distx = 1 * qp3x.real
zrx = qp3x.imag
w2x = (lam/np.pi*qp3x.imag)**.5
w02x = (lam/np.pi*zrx)**.5

z_prop2x = np.linspace(0, z_totalx, num)+distx
w_prop2x = np.zeros(num)
R_prop2x = np.zeros(num)

# Y-axis
z_totaly = L - zy
qp2y = 1/(1/R_propy[i] -1j * lam / np.pi/w_propy[i]**2)
qp3y = qp2y/(1-2*qp2y/R2y)
disty = 1 * qp3y.real
zry = qp3y.imag
w2y = (lam/np.pi*qp3y.imag)**.5
w02y = (lam/np.pi*zry)**.5

z_prop2y = np.linspace(0, z_totaly, num)+disty
w_prop2y = np.zeros(num)
R_prop2y = np.zeros(num)


for i in range(num):
    # X
    w_prop2x[i] = w02x * (1 + (lam * z_prop2x[i]/ np.pi / w02x**2) ** 2) ** .5
    R_prop2x[i] = z_prop2x[i] * (1 + (np.pi * w02x**2/lam/z_prop2x[i])**2)
    # Y
    w_prop2y[i] = w02y * (1 + (lam * z_prop2y[i]/ np.pi / w02y**2) ** 2) ** .5
    R_prop2y[i] = z_prop2y[i] * (1 + (np.pi * w02y**2/lam/z_prop2y[i])**2)


f, ax = plt.subplots(figsize = (15,6))
ax.plot(z_propx*10**2, w_propx*10**3, 'k-', label = 'X' )
ax.plot(z_propy*10**2, w_propy*10**3, 'r-', label = 'Y' )

ax.plot(z_prop2x*10**2 - z_prop2x[0]*10**2 + z_propx[-1]*10**2, w_prop2x*10**3, 'k-' )
ax.plot(z_prop2y*10**2 - z_prop2y[0]*10**2 + z_propy[-1]*10**2, w_prop2y*10**3, 'r-' )

# ax.plot([z_prop[0]*10**2, z_prop[0]*10**2],[0,3], 'r-',label = 'CM 1')
# ax.plot([z_prop[-1]*10**2, z_prop[-1]*10**2],[0,3], 'b-',label = 'CM 2')

plt.legend(loc = 'upper right')
ax.set_xlim((1.5*z_propx[0]*10**2,1.05*max(z_prop2x*10**2 - z_prop2x[0]*10**2 + z_propx[-1]*10**2) ))
ax.set_xlabel('Propagation distance [cm]')
ax.set_ylabel('Beam waist [$\mu$m]')
ax.set_title('CM separation = ' + "%.2f" % np.round( zx * 10**2,2) + ' [cm]' )
ax.text(100, .2, 'Beam waist = (' + "%.1f" % np.round( w0x * 10**6,1) +', '+ "%.1f" % np.round( w0y * 10**6,1)+ ' )[$\mu$m] ')
plt.show()






