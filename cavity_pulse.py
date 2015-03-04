import numpy as np 
from numpy.fft import * 
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.mlab import movavg
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from matplotlib.ticker import MaxNLocator, LinearLocator
from mirrors import * 
from plotter import plotter

plt.close('all')

############################Calculation######################################
##### units fs and fs**-1 ########
## generate Cavity Mirrors, assume M1 is input coupler
# M1 = mirror('Data/PR_2p4_1040v2.csv')
# M1 = mirror('Data/HighPowerPR_0p8_1070v2.csv')
M1 = mirror('Data/HighPowerPR_1p5_1070v2.csv')
M2 = mirror('Data/HighPowerHR_0p002_1070v2.csv')
M3 = mirror('Data/HighPowerHR_0p002_1070v2.csv')
M4 = mirror('Data/HighPowerHR_0p002_1070v2.csv')

## General paramters 
wavelength = np.linspace(M1.lam_min, M1.lam_max, M1.pts)
frequency = 3e2/(wavelength)
lam0 = M2.center
f0 = 3e2/lam0

##pulse parameters
tau = 120 # pulse duration

## construct transfer function
R = (M1.R() * M2.R() * M3.R() * M4.R()) 
phi = M1.Dispersion() + M2.Dispersion() + M3.Dispersion() + M4.Dispersion()

## remove offset and linear term
z = np.polyfit(wavelength - lam0, phi, 10)
phi -= z[-1] + z[-2] * (wavelength- lam0) + 0*z[-3] * (wavelength- lam0)**2

## Grid Parameters
nt = 2**12  # grid points
t = np.linspace(-40*tau,40*tau, nt)
dt = t[1]-t[0]
freq = fftfreq(nt, d = dt) + f0
wave = (3e2/(freq))

## Initial Pulse Shapes
ut_i = ((2/tau*(np.log(2)/np.pi)**.5)**.5 * np.exp(-2*np.log(2) * t ** 2/tau ** 2))
uf_i = (fft(ut_i))/(2*np.pi)**.5
startt = abs(ut_i) ** 2 / max(abs(ut_i) ** 2)
startf = fftshift(abs(uf_i) ** 2 / max(abs(uf_i) ** 2))

##################Modify with cavity#####################################
phase_shift = interp1d(frequency[::-1], phi, kind = 'linear')
phase_shift_loc = np.where((freq >3e2/wavelength.max()) & (freq < 3e2/wavelength.min()))
## Transfer function
h = (1-M1.R())**.5 / (1-R*np.exp( 1j * phase_shift(freq[phase_shift_loc])))
## Apply
uf_i[phase_shift_loc] = uf_i[phase_shift_loc] * h
uf_f = uf_i
ut_f = (ifft(uf_f))/(2*np.pi)**.5
## Steady State
stopt = abs(ut_f) ** 2 / max(abs(ut_f) ** 2)
stopf = fftshift(abs(uf_f) ** 2 / max(abs(uf_f) ** 2))
waveplot = fftshift(wave)
time_phase = np.angle(ut_f)
## Useful plotting quantities
beta  = fftshift(abs(((1-R)/ (1-R*np.exp(2* 1j * phase_shift(freq[phase_shift_loc])))))**2)
coupling_wavelength = np.sort(3e2/fftshift(freq[phase_shift_loc]))
coupling_phase = fftshift(np.angle(h[np.argsort(coupling_wavelength)]))

##################Plotting############################################
plotter(wavelength, phi, M2.GDD(), t,time_phase, startt, stopt,
    waveplot, startf, stopf, coupling_wavelength, coupling_phase, beta)

