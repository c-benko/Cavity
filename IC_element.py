IC_element.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.mlab import movavg

class IC_element:
    '''
    Data containing detailed dispersion, basic reflectivity.
    Currently, I have no data for any elements. Will add things for
    Brewster plates, gas targets, etc.
    '''
    def __init__(self, element_name = 'Data/HighPowerHR_0p002_1070.csv'):
        self.element_name = element_name
        self.df= pd.read_csv(element_name,  delimiter = ',')
        self.center = self.df['Wavelength2'][0]
        self.lam_min = 1000
        self.lam_max = 1140
        self.pts = 1000
        self.wl  = np.linspace(self.lam_min, self.lam_max,self.pts, dtype = float)

    def R(self):
        '''
        Nominal reflectivity (or loss in this case). 
        The Refelctivity is unifrom over the relevent range
        of the dispersion data. So, we assume it is constant. 
        '''
        return float(self.df['R'][0])

    def GDD(self):
        '''
        Data obtained from material properties.
        '''
        wavelength = movavg(self.df['Wavelength'][:],5)
        GDD = movavg(self.df['GDD'][:],5)
        wl = np.linspace(self.lam_min, self.lam_max,self.pts, dtype = float)
        zGDD = np.poly1d(np.polyfit(wavelength, GDD, 7))
        return zGDD(wl)
    
    def Dispersion(self):
        '''
        Returns the phase shift of the data. No polynomial dependence has
        been removed yet. 
        '''
        wavelength = movavg(self.df['Wavelength'][:],5)
        GDD = movavg(self.df['GDD'][:],5)
        wl = np.linspace(self.lam_min, self.lam_max,self.pts, dtype = float)
        zGDD = np.poly1d(np.polyfit(wavelength, GDD, 7))
        phase = np.cumsum(np.cumsum(zGDD(wl)))*(wl[1]-wl[0]) ** 2 * 9e4 * 4 * np.pi ** 2 / wl ** 4
        return phase

