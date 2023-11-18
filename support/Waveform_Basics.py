import numpy as np
import vector


def linPolEField(E0Inc, thetaInRad, phiInRad, alphaInRad):
    """Creates a 'vector' object representing an electric field with rms value E0Inc; 
    the field belongs to a linearly polarized plane wave impinging towards the 
    origin on a 3D space with incidence angles phi and theta.
    Moreover, the field, in a generalized way, is characterized by an angle
    alpha between the electric field itself and the incidence plane."""
    
    # TODO:
    # Add checks on the angle values: 0<=thetaInRad<=np.pi ; 0<=phiInRad<2*np.pi; 0<=alphaInRad<2*np.pi 
    # (alpha can also be defined negative, as long as within similar limits)
    
    if((thetaInRad==0)|(thetaInRad==np.pi)):
        EIncVector = vector.obj(x = E0Inc*np.cos(phiInRad), y = E0Inc*np.sin(phiInRad), z = 0.0)
    else:         
        # Projection of EIncVector on alpha plane (i.e. plane containing EIncVector)
        Eu_i = E0Inc*np.cos(alphaInRad)
        Ev_i = E0Inc*np.sin(alphaInRad)
        # Projection of Eu_i and Ev_i on xy-plane and z:
        Euxy_i = Eu_i*np.cos(thetaInRad)
        Euz_i = Eu_i*np.sin(thetaInRad)
        Evxy_i = Ev_i  #  Ev_i is invariant with thetaIn (Ev_i is perpendicular to pz-plane by construction)
        Evz_i = 0.0   #  Ev_i is perpendicular to pz-plane by construction (i.e. this is no magic-number)
        # Projection of Euxy_i, Euz_i, Evxy_i and Evz_i on xyz system:
        EzProj = Euz_i + Evz_i
        ExProj = (-1.0*Euxy_i*np.cos(phiInRad)  + (-1.0*Evxy_i*np.sin(phiInRad)))   
        EyProj = (-1.0*Euxy_i*np.sin(phiInRad)  + Evxy_i*np.cos(phiInRad))
        
        EIncVector = vector.obj(x = ExProj, y = EyProj, z = EzProj)

    print("Ex: ",EIncVector.x)
    print("Ey: ",EIncVector.y)
    print("Ez: ",EIncVector.z)

    return EIncVector

class planeWave:
    """The planeWave class contains the definition of a plane wave propagating in space.
        The plane wave is defined here in the frequency domain (i.e. omega), and a harmonic
        field is assumed.
        The object is defined in terms of a vector object 'E0Vector' ('Vector' library) and 
        another vector object 'waveVector'."""
    
    def __init__(self, E0Vector, waveVector):
        self.E0Vector = E0Vector
        self.waveVector = waveVector
    def __str__(self):
        return f"{self.E0Vector}*(exp(1j*dot({self.waveVector},position)))"   
    
    #Methods of class planeWave:
    
    def obsEField(self, obs):  # Compute E field at location obs
        """'self.obsEField(obs)' which returns the electric field (including the phase term, therefore 
        complex) of the plane wave at location 'obs'. """
        return self.E0Vector*np.exp(1j*self.waveVector.dot(obs))
        

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
    
    def __init__(self, E0Vector, waveVector):
        self.E0Vector = E0Vector
        self.waveVector = waveVector
    def __str__(self):
        return f"{self.E0Vector}*(exp(1j*dot({self.waveVector},position)))"   
    
    #Methods of class planeWave:
    
    def obsEField(self, obs):  # Compute E field at location obs
        """'self.obsEField(obs)' which returns the electric field (including the phase term, therefore 
        complex) of the plane wave at location 'obs'. """
        return self.E0Vector*np.exp(1j*self.waveVector.dot(obs))        
        
 
# Creating class with multiple constructors in Python
# circle.py

import math

class Circle:
    def __init__(self, radius):
        self.radius = radius

    @classmethod
    def from_diameter(cls, diameter):
        return cls(radius=diameter / 2)

    # ...
# And then:
>> from circle import Circle

>> Circle.from_diameter(84)
Circle(radius=42.0)

>> circle.area()
5541.769440932395
>> circle.perimeter()
263.89378290154264
 