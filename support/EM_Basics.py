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