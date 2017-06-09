import healpy as hpy
import numpy as np

#Calculating the Healpix discretization
Nside=16
Npix=hpy.pixelfunc.nside2npix(Nside)
ipixmin=int(Npix/2)
ipixmin=0
X=np.transpose(np.array(hpy.pixelfunc.pix2ang(Nside,np.arange(ipixmin,Npix,1))))
X[:,0] = np.pi/2-X[:,0]
#X[:,1] = X[:,1]-np.pi
#X*=180/np.pi
Num=np.size(X,0)
#---------------------------------

np.save('Discrete16.npy',X)

#Calculating the Healpix discretization
Nside=8
Npix=hpy.pixelfunc.nside2npix(Nside)
ipixmin=int(Npix/2)
ipixmin=0
X=np.transpose(np.array(hpy.pixelfunc.pix2ang(Nside,np.arange(ipixmin,Npix,1))))
X[:,0] = np.pi/2-X[:,0]
#X[:,1] = X[:,1]-np.pi
#X*=180/np.pi
Num=np.size(X,0)
#---------------------------------

np.save('Discrete8.npy',X)

#Calculating the Healpix discretization
Nside=32
Npix=hpy.pixelfunc.nside2npix(Nside)
ipixmin=int(Npix/2)
ipixmin=0
X=np.transpose(np.array(hpy.pixelfunc.pix2ang(Nside,np.arange(ipixmin,Npix,1))))
X[:,0] = np.pi/2-X[:,0]
#X[:,1] = X[:,1]-np.pi
#X*=180/np.pi
Num=np.size(X,0)
#---------------------------------

np.save('Discrete32.npy',X)