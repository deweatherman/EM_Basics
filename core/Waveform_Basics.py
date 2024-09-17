#MIT License

#Copyright (c) 2023 Mario Echeverri Bautista



import numpy as np
import scipy
import pandas as pd
        

def uerf(t,ta,tb,td,tf):
    """ Quasi-rectangular unitary pulse function.
    Definition can be found in pag 21 of:
    https://ece-research.unm.edu/summa/notes/SDAN/0041.pdf
    ta-tb = FWHM
    td = 10 to 90 % rise-time
    tf = 90 to 10 % fall-time"""
    arg1 = (t-tb)/(0.55*td)
    arg2 = (t-ta)/(0.55*tf)
    result = 0.5*(scipy.special.erf(arg1) - scipy.special.erf(arg2)) 
    return result
      
        
class waveForm:
    """The waveForm class contains the definition of a waveform; moreover
    basic signal properties in both time and frequency domain are made
    available. The class can be initialized both locally (e.g creating
    a signal and passing time and amplitude arrays) or passing a ".csv"
    with an external signal (csv containing time, amplitude values, comma
    separated)."""
    
    def __init__(self, time, amplitude, 
                 timeMultiplier = 1.0, amplitudeMultiplier = 1.0):  # Basic constructor (from time and amplitude arrays)
        """Baseline constructor; time and amplitude scales 
        are provided via multipliers. Multipliers default to 1.0.
        Inputs are: 
        two numpy arrays with time and amplitude data 
        and
        two floats for defining the scale of the time and amplitude
        respectively."""         
        self.time = time*timeMultiplier
        self.amplitude = amplitude*amplitudeMultiplier  
        if (len(time)!=len(amplitude)):
            raise TypeError("Time and amplitude must have the same lenght.")
        
    def __str__(self):
        return f"{self.amplitude}({self.time})"   
        
    @classmethod
    def from_time_csv(cls, pathTocsv, 
                      timeMultiplier = 1.0, amplitudeMultiplier = 1.0):  # Constructor 1 (from csv)
        """csv constructor; time and amplitude scales 
        are provided via multipliers. Multipliers default to 1.0.
        Inputs are :
        a string with the path to the csv file containing
        the signal definition (comma separated file),
        and
        two floats for defining the scale of the time and amplitude
        respectively."""              
        if type(pathTocsv) != str:
            raise TypeError("Path to signal file must be a string")
            
        array = pd.read_csv(pathTocsv, delimiter=",", on_bad_lines='skip') # Read from comma delimited csv, creates Pandas DataFrame 
        
        if (len(array.columns)!=2):
            raise TypeError("Waveform can only be defined with time \
                             and amplitude (2 arrays); might be extended in the future.")
            
            
        return cls(time = array.iloc[:,0].to_numpy(), 
                   amplitude = array.iloc[:,1].to_numpy(),
                   timeMultiplier = timeMultiplier, 
                   amplitudeMultiplier = amplitudeMultiplier)  # time is refered to "zero"   

    
    #Methods of class waveForm:
    
    def interpolateSignal(self, ninterpolated):  # Interpolates signal to ninterpolated samples
        # TODO: Add tests
        """ """
        time_interp = np.linspace(self.time[0], self.time[-1], ninterpolated)
        amplitude_interp = np.interp(time_interp, self.time, self.amplitude)
        
        return time_interp, amplitude_interp         
        
        
        
 
 
