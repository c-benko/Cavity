# hole loss
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator

plt.close('all')
rh = 50e-6
ti = .008
q = 13

def loss(rh, rspot):
    r = rh/rspot
    return 1 - np.exp(-4 * r ** 2)
def outcouple(rh, rspot, q):
    r = rh/rspot
    return np.exp(-4 * q ** 2 * r ** 2)

r = np.linspace(.1,3, 100)*10**-3

f, ax = plt.subplots()
ax.semilogy(r*10**3, loss(rh, r))
ax.set_xlabel('Beam Waist [mm]')
ax.set_ylabel('Loss')
ax.set_title('Hole Radius = ' + str(rh*10**6) + ' [$\mu$m]')
plt.show()

f, ax = plt.subplots()
ax.plot(r*10**3, outcouple(rh, r, q))
ax.set_xlabel('Beam Waist [mm]')
ax.set_ylabel('Outcoupled fraction')
ax.set_title('Hole Radius = ' + str(rh*10**6) + ' [$\mu$m], ' + str(q) + '$^{th}$ harmonic.' )
plt.show()



def finesse(loss, ti):
    return 2 * np.pi / (loss + ti)
def buildup(fin, ti):
    return ti * fin**2/np.pi**2

f, ax = plt.subplots()
ax.plot(r*10**3, finesse(loss(rh, r), ti))
ax.set_xlabel('Beam Waist [mm]')
ax.set_ylabel('Finesse')
ax.set_title('Hole Radius = ' + str(rh*10**6) + ' [$\mu$m], Ti = ' + str(ti))
ax2 = ax.twinx()
ax2.plot(r*10**3, buildup(finesse(loss(rh, r), ti),ti), 'r-')
ax2.set_ylabel('Buildup', color = 'r')
ax.yaxis.set_major_locator(LinearLocator(5))
ax2.yaxis.set_major_locator(LinearLocator(5))
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.grid()
plt.show()
