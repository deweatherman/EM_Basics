#MIT License

#Copyright (c) 2023 Mario Echeverri Bautista



import numpy as np
import scipy
        

def uerf(t,ta,tb,td,tf):
    """ Quasi-rectangular unitary pulse function.
    ta-tb = FWHM
    td = 10 to 90 % rise-time
    tf = 90 to 10 % fall-time"""
    arg1 = (t-tb)/(0.55*td)
    arg2 = (t-ta)/(0.55*tf)
    result = 0.5*(scipy.special.erf(arg1) - scipy.special.erf(arg2)) 
    return result
    
def spectrumRealSignal(n, fft_E):
    """spectrumRealSignal returns an array with the two-sided spectrum reconstructed
    from the one-sided version (for real signals in the time domain).
    The ordering and components are assumed as from:
    https://numpy.org/doc/stable/reference/routines.fft.html#implementation-details"""

    if(n%2==0):
        # works for even "n"
        fft_E_fminus = np.flip(np.conjugate(fft_E[1:int(n/2)]))  #  A[n/2+1:] negative freqs in decreasing neg. freq. 
        fft_all = np.concatenate((np.append(fft_E[0:int(n/2)], fft_E[int(n/2)]), fft_E_fminus))  #  cat(A[0], A[1:n/2]), A[1:n/2] pos. freqs.  
        #print(fft_all.shape)
    else:
        # Not tested
        print('Warning: if n odd, this might not work yet')
        fft_E_fminus = np.flip(np.conjugate(fft_E[1:int((n+1)/2)])) 
        fft_all = np.concatenate((fft_E[0:int((n-1)/2)+1], fft_E_fminus)) 
        
    return fft_all    
        
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
            
        array = np.genfromtxt(pathTocsv, delimiter=",")   # Read from comma delimited csv
        
        if (array.shape[1]!=2):
            raise TypeError("Waveform can only be defined with time \
                             and amplitude (2 arrays); might be extended in the future.")
            
        if (len(array[:,0]) != len(array[:,1])):
            raise TypeError("Length of time and amplitude arrays must be the same.")
            
        return cls(time = (array[:,0]-array[0,0]), 
                   amplitude = array[:,1],
                   timeMultiplier = timeMultiplier, 
                   amplitudeMultiplier = amplitudeMultiplier)  # time is refered to "zero"   

    
    #Methods of class waveForm:
    
    def interpolateSignal(self, ninterpolated):  # Interpolates signal to ninterpolated samples
        # TODO: Add tests
        """ """
        time_interp = np.linspace(0, (self.time[-1]-self.time[0]), ninterpolated)
        amplitude_interp = np.interp(time_interp, self.time, self.amplitude)
        
        return time_interp, amplitude_interp         
        
        
        
 
# Creating class with multiple constructors in Python
# circle.py

#import math
#
#class Circle:
#    def __init__(self, radius):
#        self.radius = radius
#
#    @classmethod
#    def from_diameter(cls, diameter):
#        return cls(radius=diameter / 2)

    # ...
# And then:


#>> from circle import Circle
#
#>> Circle.from_diameter(84)
#Circle(radius=42.0)
#
#>> circle.area()
#5541.769440932395
#>> circle.perimeter()
#263.89378290154264
 
