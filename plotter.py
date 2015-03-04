#plotter.py
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, LinearLocator

def plotter(wavelength, phi, GDD, t, time_phase,startt, stopt,
    waveplot, startf, stopf, coupling_wavelength, coupling_phase, beta):

    # IC phase shift
    f, ax =  plt.subplots(2)
    ax[0].plot(wavelength, phi*10**3, 'k-')
    ax[0].set_xlabel('Wavelength [nm]')
    ax[0].set_ylabel('Phase Shift [mrad]')
    ax[0].set_ylim(-100,100)
    # GDD 
    ax[1].plot(wavelength, GDD, 'k-')
    ax[1].yaxis.set_major_locator(LinearLocator(5))
    ax[1].set_xlabel('Wavelength [nm]')
    ax[1].set_ylabel(r'GDD [fs$^{2}$]')
    ax[1].set_ylim(-50,50)
    plt.show()

    f, ax = plt.subplots(2)
    ax[0].plot(t,startt ,'-b', label  = 'Input')
    ax[0].plot(t,stopt ,'-k', label  = 'Circulating')
    ax[0].yaxis.set_major_locator(LinearLocator(5))
    ax2 = ax[0].twinx()
    ax2.plot(t, time_phase, 'r-')
    ax2.yaxis.set_major_locator(LinearLocator(5))
    ax2.set_ylabel('Phase [rad]', color = 'r')
    ax2.grid()
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    ax[0].legend( loc = 'upper left')
    ax[0].set_xlim(-250,250)
    ax[0].set_xlabel('Time [fs]')
    ax[0].set_ylabel('Normalizezd Power')
    
    #Wavelength
    #input
    ax[1].plot(waveplot, startf ,'-b', label = 'Input') 
    ax[1].plot(waveplot, stopf ,'-k', label = 'Circulating') 
    ax[1].yaxis.set_major_locator(LinearLocator(5))
    #phase
    ax2 = ax[1].twinx()
    ax2.plot(coupling_wavelength, coupling_phase,  'r')
    ax2.set_ylabel('Phase [rad]', color = 'r')
    ax2.set_ylim(-4,4)
    ax2.yaxis.set_major_locator(LinearLocator(5))
    ax2.grid()
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    ax[1].legend()
    ax2.set_xlim(1000,1140)
    ax[1].set_xlabel('Wavelength [nm]')
    ax[1].set_ylabel('Spectral Power')
    ax[1].legend()
    plt.show()
    
    # coupling efficiency and phase
    f, ax = plt.subplots(2)
    ax[0].plot(coupling_wavelength, coupling_phase, 'k-')
    ax[0].set_ylabel('Phase [rad]')
    ax[0].set_xlabel('Wavelength [nm]')
    ax[1].plot(coupling_wavelength, beta , 'k-')
    ax[1].set_ylabel('Coupling Efficiency')
    ax[1].set_xlabel('Wavelength [nm]')
    plt.show()
    
