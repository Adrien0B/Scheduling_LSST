
# coding: utf-8

# In[1]:

import healpy as hpy
import datetime
import os

from ACO_glouton import *
from functions import conversion
from functions import deconversion
from matplotlib import animation
from animation_path import *

date_aux="2017/03/09 23:00:00"


# In[2]:

#Initializing the observer instance
obs=ephem.Observer()
obs.lat="-33:27:00"
obs.lon="-70:40:00"
obs.date=date_aux
#-------------------


# In[3]:

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


# In[4]:

#Times since last visit
#T = npr.randint(6, size=Num)+1
if Num == 768:
    T=np.load('Times8_uniform.npy')
elif Num == 1536:
	T=np.load('Times.npy')
elif Num == 6144:
	T=np.load('Times32.npy')
elif Num == 3072:
    T=np.load('Times_complete.npy')
elif Num == 12288:
    T=np.load('Times32_complete.npy')
else:
    print 'Times since last observations created randomly'
    T=npr.randint(0,6,Num)
#-------------------------


# In[ ]:

ACO=ACOScheduleGlouton(X,obs,15,T)
print ACO.m


# In[ ]:

ACO.RunACO_Pheromone(0,1)
print len(ACO.BPS)


# In[ ]:

if Nside==32:
    markersize=25
elif Nside==8:
    markersize=400
else:
    markersize=100


timenow = datetime.datetime.now()


# In[ ]:

def print_list(a):
    print '[',
    for i in a[:-1]:
        print i,',',
    print a[-1],']'


a = np.sort(np.transpose(np.array([int(i[1]*10)/10.0 for i in ACO.BPS])))
b = np.sort(np.transpose(np.array([i[2] for i in ACO.BPS])))[::-1]
print_list(a)
print_list(b)
title = "videos/%s-%s-%s_%s-%s-%s_%sIter_Pareto_Values" % (timenow.year,timenow.month,timenow.day,timenow.hour,timenow.minute,timenow.second,ACO.AntIterations)
np.savetxt("qhflsdjfn.txt",(a,b),fmt='%s',delimiter=',',newline=']\n[',header=' ')
lines = open("qhflsdjfn.txt").readlines()
open(title, 'w').writelines(lines[1:-1])
os.remove("qhflsdjfn.txt")


# In[ ]:

if(0):
    print timenow
    fig = ACO.PlotParetoFront(title="%s-%s-%s_%s-%s-%s_%sIter_Pareto_Front" % (timenow.year,timenow.month,timenow.day,timenow.hour,timenow.minute,timenow.second,ACO.AntIterations))
    fig.savefig("videos/%s-%s-%s_%s-%s-%s_%sIter_Pareto_Front" % (timenow.year,timenow.month,timenow.day,timenow.hour,timenow.minute,timenow.second,ACO.AntIterations))
    plt.clf()
    for i in range(len(ACO.BPS)):
        [schedAA,schedDR]=ACO.AZALT(ACO.BPS[i][0])
        print "Saving equatorial for solution "+str(i+1)+"/"+str(len(ACO.BPS))
        animation_path(schedDR,(-10,-145,0),"%s-%s-%s_%s-%s-%s_%sIter_Sol%s_%sObservations_T%s_O%s_equatorial" % (timenow.year,timenow.month,timenow.day,timenow.hour,timenow.minute,timenow.second,ACO.AntIterations,i,np.size(ACO.BPS[i][0],0),ACO.BPS[i][2],int(ACO.BPS[i][1]*10)/10.0),marker='.',s=markersize)
        print "Saving horizontal for solution "+str(i+1)+"/"+str(len(ACO.BPS))
        plt.clf()
        animation_path(schedAA,(0,90,0),"%s-%s-%s_%s-%s-%s_%sIter_Sol%s_%sObservations_T%s_O%s_horizontal" % (timenow.year,timenow.month,timenow.day,timenow.hour,timenow.minute,timenow.second,ACO.AntIterations,i,np.size(ACO.BPS[i][0],0),ACO.BPS[i][2],int(ACO.BPS[i][1]*10)/10.0),marker='.',s=markersize)
        plt.clf()
        plt.close('all')

