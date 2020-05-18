# -*- coding: utf-8 -*-
"""
Demo for Fresnel propagation through multiple thin scattering layers. Generates a point
source at fixed position, then N thin layers with random phase distributions
at fixed distances from the source. Source emits an spherical wave that travels
through the whole system. Finally, an image of the speckle pattern is captured
with a CCD (intensity)

@author: Fer
"""
# %% Import stuff
# for plotting
from matplotlib import pyplot as plt 
# for manipulating images/matrices
from skimage import transform as tf
# for math/physics
import numpy as np
from numpy import matlib as mb
from numpy.fft import fft2 as fft2
from numpy.fft import ifft2 as ifft2
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift


#%% Function definitions

def fresnelProp(z0,dx,dy,propdist,wavelength,padding):
    """
    fresnelProp does the Fresnel propagation of a field.
    z0: initial field at the aperture plane
    dx: size of the pixel in x-dim (meters)
    dy: size of the pixel in y-dim (meters)
    propdist: propagation distance (meters)
    wavelength: wavelenght of illumination (meters)
    padding: padding size (factor, 1 means no padding)
    
    returns the propagated field (zProp) and the Fresnel number (fNum)
    """
    #%%Physical definitions
    k = 2*np.pi/wavelength
    #%%Sizes of the field and aperture
    (x,y) = z0.shape #Field size (in pixels)
    (Sx,Sy) = (x*dx,y*dy) #Sizes of the aperture (length units)
    #%%Generate frequency mesh (after padding)
    (dfx,dfy) = (1/Sx,1/Sy) #%Frequency element in kx,ky axes
    #Frequencies in both axes
    Fx = mb.repmat(np.array([dfx*np.linspace(-0.5*x,0.5*x,x*padding)],ndmin=2),x*padding,1)
    Fy = mb.repmat(np.array([dfy*np.linspace(0.5*y,-0.5*y,y*padding)],ndmin=2).T,1,y*padding)
    pxNumBig = Fx.size; #Total number of elements (pixels) after padding
    #%%Calculate pad size and apply pad to the initial field
    (padsizeX,padsizeY) = np.array((round((padding-1)*x/2),round((padding-1)*y/2))).astype(int)
    z0padded = np.pad(z0,((padsizeX, padsizeY),(padsizeX, padsizeY)),'constant')
    #%% Propagation in Fourier space
    #Build Fresnel propagator
    H = np.exp(1j*k*propdist)*np.exp(-1j*np.pi*wavelength*propdist*(Fx**2 + Fy**2));
    #Do the FT of initial field
    zFT = ifftshift(fft2(fftshift(z0padded)))/np.sqrt(pxNumBig)
    #Do the propagation in Fourier domain. Multiply in fourier domain (convolution in spatial domain), then IFT
    zProp = fftshift(ifft2(ifftshift(zFT*H)))
    #Calculate Fresnel Number
    fNum = (Sx/2)**2/(wavelength*propdist);
    return zProp,fNum

def circAp(totalSize,radius):
    """
    circAp generates a circular aperture, given the total size of the mask and
    the radius of the aperture (between 0 and 1). Also provides the coordinates
    of the mask in polar coordinate system
    """
    x = np.linspace(-1,1,totalSize)
    y = np.linspace(1,-1,totalSize)
    xx, yy = np.meshgrid(x,y)
    
    r = np.abs(xx + 1j*yy)
    theta = np.angle(xx + 1j*yy)
    
    mask = r<radius
    
    return mask,r,theta


# %% Define global parameters
class coordinates:
    size = 128 #number of spatial positions (in pixels)
    x = np.linspace(-1,1,size)  #generate axis (in pixels)
    apsize = 2e-3    #aperture_size (in meters)
    dx = apsize/size #pixel pitch in X
    dy = dx #pixel pitch in Y
    #generate coordinates (in length units)
    X,Y = np.meshgrid(x*apsize/2,x*apsize/2)
    # apradius = 1 #relative size of the aperture radius
    # #Calculate pupil, get polar coords
    # idx, theta, r = circAp(size,apradius) 
    pass

coords = coordinates

# %% Generate point source
class pointsource:
    zPos = 20e-2 #distance between source and first layer
    xPos, yPos = (0*1e-5,0*1e-5) #position in the XY plane
    wavelength = 532e-9 #define wavelength
    k = 2*np.pi/wavelength #define k
    amplitude = 1e0 #aplitude of the point source
    Rdist = coords.apsize**2/(2*wavelength) #Rayleigh distance
    pass

#initialize point source properties
source = pointsource

# %% Definition of scattering layer properties and preallocation
np.random.seed(0) #Set random seed so we always get the same scattering layers

#class for storing scattering layers properties
class scatprops:
    #number of superpixels over the FoV 
    #(to make random phases with detail size bigger than 1 pixel)
    pxfrac = 32
    numlayers = 6   #number of scattering layers
    #padding for each fresnel propagation between scat.layers. 
    #length must be: numlayers-1. Manual tuning for now
    paddings = np.array((2,1,1,1,1),dtype='i4') 
    interlayerdist = 5e-3   #interlayer distance
    propdist = 15e-2    #distance between last scattering layer and camera

#initialize scatterig layers properties
scprops = scatprops

#define class for storing each scattering layer
class scat:
    seed = np.array((),dtype='f4')  #random seed
    efield = np.array((),dtype='f4')    #field using random seed
    pass

#Initialize layers (empty)
sclayer = [scat() for i in range(scprops.numlayers)]

# %% Build spherical wavefront at the plane of the first scattering layer
class wave:
    dt = np.array((),dtype='f4')
    amplitude = np.array((),dtype='f4')
    phase = np.array((),dtype='f4')
    efield = np.array((),dtype='f4')
    pass

#Generate waves at all crucial planes (empty)
waves = [wave() for i in range(scprops.numlayers*2+1)]

#Generate first field (from point source to first scat.layer)
#Generate distance travelled (dt) between source and points in the first scattering plane
waves[0].dt = np.sqrt(source.zPos**2 + (coords.X-source.xPos)**2 + (coords.Y-source.yPos)**2)
waves[0].amplitude = source.amplitude/waves[0].dt #Calculate wave amplitude
waves[0].phase = source.k*waves[0].dt #Calculate wave phase
#Combine amplitude and phase into field expression
waves[0].efield = waves[0].amplitude*np.exp(1j*waves[0].phase) 

# %% Propagation of the field between scat.layers

k = 1   #dummy variable because number of layers != number of fields I store
for i in range(scprops.numlayers):
    #Generate random seed between -1 and 1
    #[resize] changes 'spatial resolution' of the seed
    sclayer[i].seed = np.random.randn(scprops.pxfrac,scprops.pxfrac)
    sclayer[i].seed = tf.resize(sclayer[i].seed,waves[k-1].efield.shape,order=0)
    #Generate phase distribution of the scattering layer
    sclayer[i].efield = np.ones(sclayer[i].seed.shape,dtype=float)*np.exp(1j*sclayer[i].seed*np.pi)
    #Multiply field before the scat.layer with the field of the scat.layer
    waves[k].efield = waves[k-1].efield*sclayer[i].efield
    k = k+1
    #Propagate field to the next layer (except if you are on the last layer)
    if i<scprops.numlayers:
        waves[k].efield, _ = fresnelProp(waves[k-1].efield,coords.dx,coords.dy,scprops.interlayerdist,source.wavelength,scprops.paddings[i-1])
        k = k+1

# Final propagation (from last scat.layer to camera)
waves[k-1].efield, _ = fresnelProp(waves[k-2].efield,coords.dx,coords.dy,scprops.propdist,source.wavelength,2)
waves[k-1].intensity = np.abs(waves[k-1].efield)    
    
#%%Show speckle image 
plt.figure(dpi=300)
plt.imshow(waves[k-1].intensity,cmap='hot')
plt.axis('on')
plt.title('Final intensity on camera')
plt.colorbar()
plt.show()
