# -*- coding: cp1252 -*-


import matplotlib
from matplotlib.axes import Axes
from matplotlib.patches import Circle
from matplotlib.path import Path
from matplotlib.ticker import NullLocator, Formatter, FixedLocator
from matplotlib.transforms import Affine2D, BboxTransformTo, Transform
from matplotlib.projections import register_projection
import matplotlib.spines as mspines
import matplotlib.axis as maxis
from transformation import transform

import ephem
import numpy as np
import numpy.random as npr
#from ranking import ranking

def conversion2(lat,dec,ar,LST):
	"Computes [Azimuth,Altitude]. It can recieve  vectors."
	pi = np.pi
	deg2rad = pi / 180.
	hr2rad = pi / 12.
	ha = np.mod(np.array( LST - ar),2*pi)
	cosz = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(ha)
	sinz = np.sqrt(1 - cosz**2)
	z = np.arccos(cosz)
	ALT=pi/2-z

	cosAZ = (np.sin(dec) - cosz * np.sin(lat)) / (sinz * np.cos(lat))
	AZ =np.array( np.arccos(cosAZ))
	for i in range(AZ.size):
		if (AZ.size==1):
			if (ha<=pi):
				AZ=-AZ
			break
		if (ha[i] <= pi) :
			AZ[i] = -AZ[i]

	return np.array([AZ,ALT])

def conversion4(lat, dec, ra, LST):
    "Computes [Azimuth,Altitude]. It can recieve  vectors."
    H = LST-ra

    tanAzix = np.cos(lat)*np.sin(dec)-np.sin(lat)*np.cos(dec)*np.cos(H)
    tanAziy = np.cos(dec)*np.sin(H)

    AZ = -np.arctan2(tanAziy,tanAzix)

    sinh = np.sin(lat)*np.sin(dec)+np.cos(lat)*np.cos(dec)*np.cos(H)

    ALT = np.arcsin(sinh)

    return np.array([AZ, ALT])

def conversion5(obs, dec, ra):
    """
    AZALT=np.array(map(lambda x,y : transform((x,y),'equatorial','horizon',obs),np.array(ra),np.array(dec)))
    AZ = AZALT[:,0]
    ALT = AZALT[:,1]
    return np.array([AZ,ALT])
    """
    try:
        AZ = []
        ALT = []
        for i in range(len(dec)):
            a = transform((ra[i],dec[i]),'equatorial','horizon',obs)
            AZ.append(a[0])
            ALT.append(a[1])
        return np.array([AZ,ALT])
    except:
        a = transform((ra, dec), 'equatorial', 'horizon', obs)
        return np.array([a[0],a[1]])

def deconversion4(obs, az, h):
    try:
        RA = []
        DEC = []
        for i in range(len(az)):
            a = obs.radec_of(az[i], h[i])
            RA.append(a[0])
            DEC.append(a[1])

        return np.array([RA, DEC])
    except:
        a = obs.radec_of(az, h)
        return np.array([a[0], a[1]])

def deconversion5(obs, az, h):
    try:
        RA = []
        DEC = []
        for i in range(len(az)):
            a = transform((az[i],h[i]),'horizon','equatorial',obs)
            RA.append(a[0])
            DEC.append(a[1])
        return np.array([RA, DEC])
    except:
        a = transform((az, h), 'horizon', 'equatorial', obs)
        return np.array([a[0],a[1]])

def conversion(obs, dec, ra):
    #return conversion2(obs.lat,dec,ra,obs.sidereal_time())
    return conversion4(obs.lat,dec,ra,obs.sidereal_time())
    #return conversion5(obs,dec,ra)

def deconversion(obs,az,h):
    #return deconversion4(obs,az,h)
    return deconversion5(obs,az,h)

def factibles( X,obs,actual ):
	"Returns factible points from a list X=[DEC,RA] for the observer. It considers angle distance from the moon and the angle above the horizon."
	# Parameters
	dist_moon=40*np.pi/180
	above_horizon=40*np.pi/180

	#Initialize Moon instance
	moon=ephem.Moon()
	moon.compute(obs)

	#Compute [AZ,ALT]
	X_obs=np.transpose(np.array(conversion(obs,X[:,0],X[:,1]))) #[AZ,ALT] of X for the observer.

	#One if factible, zero if not
	I=np.ones(np.size(X,0)) #Idicates if factible or not.	
	for i in range(0,np.size(I,0)):
		if (ephem.separation((X_obs[i,0],X_obs[i,1]),(moon.az,moon.alt))<=dist_moon or X_obs[i,1]<=above_horizon):
		#if (X_obs[i,1]<=above_horizon):
			I[i]=0


	return I

def Time_dist(X,observer,actual):
    az_vel=10./(24*3600)
    alt_vel=10./(24*3600)
    X_obs=np.transpose(np.array(conversion(observer,X[:,0],X[:,1]))) #[AZ,ALT] of X for the observer.
    actual_AA=np.transpose(np.array(conversion(observer,actual[0],actual[1]))) #[AZ,ALT] of actual point.
    diference=np.abs(X_obs - np.ones((np.size(X,0),1))*actual_AA)
    diference[:,0]=np.min((diference[:,0],2*np.pi-diference[:,0]),0)
    return np.maximum(az_vel*diference[:,0],alt_vel*diference[:,1])

def Time_dist_Multi(X,obs,actual,i,j):
    observer = ephem.Observer()
    observer.date = obs[0]
    observer.epoch = obs[1]
    observer.lat = obs[2]
    observer.lon = obs[3]
    observer.elevation = obs[4]
    observer.temp = obs[5]
    observer.pressure = obs[6]
    az_vel=10./(24*3600)
    alt_vel=10./(24*3600)
    X_obs=np.transpose(np.array(conversion(observer,X[:,0],X[:,1]))) #[AZ,ALT] of X for the observer.
    actual_AA=np.transpose(np.array(conversion(observer,actual[0],actual[1]))) #[AZ,ALT] of actual point.
    diference=np.abs(X_obs - np.ones((np.size(X,0),1))*actual_AA)
    diference[:,0]=np.min((diference[:,0],2*np.pi-diference[:,0]),0)
    return [i,j,np.maximum(az_vel*diference[:,0],alt_vel*diference[:,1]),X,observer,actual]

def Time_dist_map(X,obs,actual,i,j):
    actual = actual[j,:]
    observer = ephem.Observer()
    observer.date = obs[0]
    observer.epoch = obs[1]
    observer.lat = obs[2]
    observer.lon = obs[3]
    observer.elevation = obs[4]
    observer.temp = obs[5]
    observer.pressure = obs[6]
    az_vel=10./(24*3600)
    alt_vel=10./(24*3600)
    X_obs=np.transpose(np.array(conversion(observer,X[:,0],X[:,1]))) #[AZ,ALT] of X for the observer.
    actual_AA=np.transpose(np.array(conversion(observer,actual[0],actual[1]))) #[AZ,ALT] of actual point.
    diference=np.abs(X_obs - np.ones((np.size(X,0),1))*actual_AA)
    diference[:,0]=np.min((diference[:,0],2*np.pi-diference[:,0]),0)
    return [i,j,np.maximum(az_vel*diference[:,0],alt_vel*diference[:,1]),X,observer,actual]

def Time_dist_pipe(X,obs,actual,i,j,conn):
    observer = ephem.Observer()
    observer.date = obs[0]
    observer.epoch = obs[1]
    observer.lat = obs[2]
    observer.lon = obs[3]
    observer.elevation = obs[4]
    observer.temp = obs[5]
    observer.pressure = obs[6]
    az_vel=10./(24*3600)
    alt_vel=10./(24*3600)
    X_obs=np.transpose(np.array(conversion(observer,X[:,0],X[:,1]))) #[AZ,ALT] of X for the observer.
    actual_AA=np.transpose(np.array(conversion(observer,actual[0],actual[1]))) #[AZ,ALT] of actual point.
    diference=np.abs(X_obs - np.ones((np.size(X,0),1))*actual_AA)
    diference[:,0]=np.min((diference[:,0],2*np.pi-diference[:,0]),0)
    conn.send([i,j,np.maximum(az_vel*diference[:,0],alt_vel*diference[:,1]),X,observer,actual])

def Factible_Schedule(X,observer,T):
	observation_time=30
	Init_time=observer.date
	sun=ephem.Sun()
	sun.compute(observer)
	EndNight=observer.next_rising(sun)
	LengthNight=EndNight-observer.date
	print LengthNight	
	print EndNight


	#n=np.int(np.floor(np.sqrt(N)))
	N=np.ceil(LengthNight*3600*24/30)
	Sched=np.zeros((N,2))
	SchedAZALT=np.zeros((N,2))
	
	
	Sched[1,:]=[0,np.pi/2]  #The first step is the zenit. It may be unfeasible. It will be removed.
	Rank1=np.zeros((np.size(X,0),10))

	for ii in range(0,10):
		[Rank1_aux,Rank2]=ranking(X,observer,np.size(X,0),T,1)
                Rank1[:,ii]=Rank1_aux
		observer.date+=LengthNight/10
	observer.date-=LengthNight

	i=0
	while observer.date<EndNight:
		i+=1
		Night_discrete=np.floor(10*(observer.date-Init_time)/LengthNight)
		Fact=factibles(X,observer,Sched[i-1,:])
		#calculate distance
		Dist=Time_dist(X[Fact,:],observer,Sched[i-1,:])
		#choose randomly
		Valor=(((T[Fact]>=3.65)+0.01)*T[Fact])*1000/(5*np.abs(Rank1[Fact,Night_discrete])+5*Dist+1)

		j=np.argmax(Valor)
		Sched[i,:]=X[Fact[j],:]
		[SchedAZALT[i,0],SchedAZALT[i,1]]=conversion(observer, Sched[i,0], Sched[i,1])
		Rank2[Fact[j]]=np.size(X,0)+i
		desplacement_time=Dist[j]
		T[Fact[j]]=0
		
		observer.date+= (desplacement_time + observation_time)/(3600*24)
		T+= (desplacement_time + observation_time)/(3600*24)
	return [Sched[1:i-1,:], SchedAZALT[1:i-1,:]]



# def ranking(X, obs, N, T ,vel):
# #Ranking 1: seg�n angulos zenitales.
# #Ranking 2: segun tiempo de ultima visita.
#
# #X: arreglo de pares (DEC, RA).
# #obs: observador (tiempo, latitud, longitud)
# #N: n�mero de puntos a obtener.
# #T: arreglo con los tiempos desde la �ltima visita.
# #vel: velocidad de movimiento del telescopio.
#
# #output: rk1, rk2.
#
#     lat = obs.lat
#     LST = obs.sidereal_time()
#
#     #calcular el �ngulo cenital de cada punto en X.
#     dec = X[:,0]
#     ra = X[:,1]
#     [az, alt] = conversion(obs, dec, ra)
#     zen1 = np.pi/2 - alt #�ngulo cenital
#
# #Using lists zen, find Ranking 1. We want the largest values, penalizing great slew time:
#
#     rk1 = np.zeros((1,N), dtype=np.int) #initialize rk1
#     zen = np.transpose(zen1)
#
#     #ind = np.argpartition(zen, -N)[-N:] #revisar comando.
#     #alternativa: usar modulo bottleneck. Debiera ser m�s r�pido.
#     temp = np.argpartition(-zen, N-1)
#     ind = temp[:N]
#     rk1 = np.transpose(ind)
#
# #Using list T, find Ranking 2.
#     Y = T
#     #idea: si tengo 2 puntos de mismo tiempo sin visita (T), prefiero el de menor
#     #tiempo de movimiento. Si tengo 2 puntos a mismo tiempo de movimiento, prioriza
#     #el que fue visitado hace m�s tiempo.
#
#     rk2 = np.zeros((1,N), dtype=np.int) #initialize rk2
#
#     #recorrer Y y determinar los �ndices de ranking.
#     #ind = np.argpartition(Y, -N)[-N:]
#     ## rk2 = ind[np.argsort(Y[ind])]
#
#     temp = np.argpartition(Y, N-1)
#     ind = temp[:N]
#     rk2 = np.transpose(ind)
#
#     return [zen1, rk2]




# class HammerAxes(Axes):
#     """
#     A custom class for the Aitoff-Hammer projection, an equal-area map
#     projection.
#
#     http://en.wikipedia.org/wiki/Hammer_projection
#     """
#     # The projection must specify a name.  This will be used be the
#     # user to select the projection, i.e. ``subplot(111,
#     # projection='custom_hammer')``.
#     name = 'custom_hammer'
#
#     def __init__(self, *args, **kwargs):
#         Axes.__init__(self, *args, **kwargs)
#         self.set_aspect(0.5, adjustable='box', anchor='C')
#         self.cla()
#
#     def _init_axis(self):
#         self.xaxis = maxis.XAxis(self)
#         self.yaxis = maxis.YAxis(self)
#         # Do not register xaxis or yaxis with spines -- as done in
#         # Axes._init_axis() -- until HammerAxes.xaxis.cla() works.
#         #self.spines['hammer'].register_axis(self.yaxis)
#         self._update_transScale()
#
#     def cla(self):
#         """
#         Override to set up some reasonable defaults.
#         """
#         # Don't forget to call the base class
#         Axes.cla(self)
#
#         # Set up a default grid spacing
#         self.set_longitude_grid(30)
#         self.set_latitude_grid(15)
#         self.set_longitude_grid_ends(75)
#
#         # Turn off minor ticking altogether
#         self.xaxis.set_minor_locator(NullLocator())
#         self.yaxis.set_minor_locator(NullLocator())
#
#         # Do not display ticks -- we only want gridlines and text
#         self.xaxis.set_ticks_position('none')
#         self.yaxis.set_ticks_position('none')
#
#         # The limits on this projection are fixed -- they are not to
#         # be changed by the user.  This makes the math in the
#         # transformation itself easier, and since this is a toy
#         # example, the easier, the better.
#         Axes.set_xlim(self, -np.pi, np.pi)
#         Axes.set_ylim(self, -np.pi / 2.0, np.pi / 2.0)
#
#     def _set_lim_and_transforms(self):
#         """
#         This is called once when the plot is created to set up all the
#         transforms for the data, text and grids.
#         """
#         # There are three important coordinate spaces going on here:
#         #
#         #    1. Data space: The space of the data itself
#         #
#         #    2. Axes space: The unit rectangle (0, 0) to (1, 1)
#         #       covering the entire plot area.
#         #
#         #    3. Display space: The coordinates of the resulting image,
#         #       often in pixels or dpi/inch.
#
#         # This function makes heavy use of the Transform classes in
#         # ``lib/matplotlib/transforms.py.`` For more information, see
#         # the inline documentation there.
#
#         # The goal of the first two transformations is to get from the
#         # data space (in this case longitude and latitude) to axes
#         # space.  It is separated into a non-affine and affine part so
#         # that the non-affine part does not have to be recomputed when
#         # a simple affine change to the figure has been made (such as
#         # resizing the window or changing the dpi).
#
#         # 1) The core transformation from data space into
#         # rectilinear space defined in the HammerTransform class.
#         self.transProjection = self.HammerTransform()
#
#         # 2) The above has an output range that is not in the unit
#         # rectangle, so scale and translate it so it fits correctly
#         # within the axes.  The peculiar calculations of xscale and
#         # yscale are specific to a Aitoff-Hammer projection, so don't
#         # worry about them too much.
#         xscale = 2.0 * np.sqrt(2.0) * np.sin(0.5 * np.pi)
#         yscale = np.sqrt(2.0) * np.sin(0.5 * np.pi)
#         self.transAffine = Affine2D() \
#             .scale(0.5 / xscale, 0.5 / yscale) \
#             .translate(0.5, 0.5)
#
#         # 3) This is the transformation from axes space to display
#         # space.
#         self.transAxes = BboxTransformTo(self.bbox)
#
#         # Now put these 3 transforms together -- from data all the way
#         # to display coordinates.  Using the '+' operator, these
#         # transforms will be applied "in order".  The transforms are
#         # automatically simplified, if possible, by the underlying
#         # transformation framework.
#         self.transData = \
#             self.transProjection + \
#             self.transAffine + \
#             self.transAxes
#
#         # The main data transformation is set up.  Now deal with
#         # gridlines and tick labels.
#
#         # Longitude gridlines and ticklabels.  The input to these
#         # transforms are in display space in x and axes space in y.
#         # Therefore, the input values will be in range (-xmin, 0),
#         # (xmax, 1).  The goal of these transforms is to go from that
#         # space to display space.  The tick labels will be offset 4
#         # pixels from the equator.
#         self._xaxis_pretransform = \
#             Affine2D() \
#             .scale(1.0, np.pi) \
#             .translate(0.0, -np.pi)
#         self._xaxis_transform = \
#             self._xaxis_pretransform + \
#             self.transData
#         self._xaxis_text1_transform = \
#             Affine2D().scale(1.0, 0.0) + \
#             self.transData + \
#             Affine2D().translate(0.0, 4.0)
#         self._xaxis_text2_transform = \
#             Affine2D().scale(1.0, 0.0) + \
#             self.transData + \
#             Affine2D().translate(0.0, -4.0)
#
#         # Now set up the transforms for the latitude ticks.  The input to
#         # these transforms are in axes space in x and display space in
#         # y.  Therefore, the input values will be in range (0, -ymin),
#         # (1, ymax).  The goal of these transforms is to go from that
#         # space to display space.  The tick labels will be offset 4
#         # pixels from the edge of the axes ellipse.
#         yaxis_stretch = Affine2D().scale(2*np.pi, 1.0).translate(-np.pi, 0.0)
#         yaxis_space = Affine2D().scale(1.0, 1.1)
#         self._yaxis_transform = \
#             yaxis_stretch + \
#             self.transData
#         yaxis_text_base = \
#             yaxis_stretch + \
#             self.transProjection + \
#             (yaxis_space +
#              self.transAffine +
#              self.transAxes)
#         self._yaxis_text1_transform = \
#             yaxis_text_base + \
#             Affine2D().translate(-8.0, 0.0)
#         self._yaxis_text2_transform = \
#             yaxis_text_base + \
#             Affine2D().translate(8.0, 0.0)
#
#     def get_xaxis_transform(self, which='grid'):
#         """
#         Override this method to provide a transformation for the
#         x-axis grid and ticks.
#         """
#         assert which in ['tick1', 'tick2', 'grid']
#         return self._xaxis_transform
#
#     def get_xaxis_text1_transform(self, pixelPad):
#         """
#         Override this method to provide a transformation for the
#         x-axis tick labels.
#
#         Returns a tuple of the form (transform, valign, halign)
#         """
#         return self._xaxis_text1_transform, 'bottom', 'center'
#
#     def get_xaxis_text2_transform(self, pixelPad):
#         """
#         Override this method to provide a transformation for the
#         secondary x-axis tick labels.
#
#         Returns a tuple of the form (transform, valign, halign)
#         """
#         return self._xaxis_text2_transform, 'top', 'center'
#
#     def get_yaxis_transform(self, which='grid'):
#         """
#         Override this method to provide a transformation for the
#         y-axis grid and ticks.
#         """
#         assert which in ['tick1', 'tick2', 'grid']
#         return self._yaxis_transform
#
#     def get_yaxis_text1_transform(self, pixelPad):
#         """
#         Override this method to provide a transformation for the
#         y-axis tick labels.
#
#         Returns a tuple of the form (transform, valign, halign)
#         """
#         return self._yaxis_text1_transform, 'center', 'right'
#
#     def get_yaxis_text2_transform(self, pixelPad):
#         """
#         Override this method to provide a transformation for the
#         secondary y-axis tick labels.
#
#         Returns a tuple of the form (transform, valign, halign)
#         """
#         return self._yaxis_text2_transform, 'center', 'left'
#
#     def _gen_axes_patch(self):
#         """
#         Override this method to define the shape that is used for the
#         background of the plot.  It should be a subclass of Patch.
#
#         In this case, it is a Circle (that may be warped by the axes
#         transform into an ellipse).  Any data and gridlines will be
#         clipped to this shape.
#         """
#         return Circle((0.5, 0.5), 0.5)
#
#     def _gen_axes_spines(self):
#         return {'custom_hammer': mspines.Spine.circular_spine(self,
#                                                               (0.5, 0.5), 0.5)}
#
#     # Prevent the user from applying scales to one or both of the
#     # axes.  In this particular case, scaling the axes wouldn't make
#     # sense, so we don't allow it.
#     def set_xscale(self, *args, **kwargs):
#         if args[0] != 'linear':
#             raise NotImplementedError
#         Axes.set_xscale(self, *args, **kwargs)
#
#     def set_yscale(self, *args, **kwargs):
#         if args[0] != 'linear':
#             raise NotImplementedError
#         Axes.set_yscale(self, *args, **kwargs)
#
#     # Prevent the user from changing the axes limits.  In our case, we
#     # want to display the whole sphere all the time, so we override
#     # set_xlim and set_ylim to ignore any input.  This also applies to
#     # interactive panning and zooming in the GUI interfaces.
#     def set_xlim(self, *args, **kwargs):
#         Axes.set_xlim(self, -np.pi, np.pi)
#         Axes.set_ylim(self, -np.pi / 2.0, np.pi / 2.0)
#     set_ylim = set_xlim
#
#     def format_coord(self, lon, lat):
#         """
#         Override this method to change how the values are displayed in
#         the status bar.
#
#         In this case, we want them to be displayed in degrees N/S/E/W.
#         """
#         lon = np.degrees(lon)
#         lat = np.degrees(lat)
#         if lat >= 0.0:
#             ns = 'N'
#         else:
#             ns = 'S'
#         if lon >= 0.0:
#             ew = 'E'
#         else:
#             ew = 'W'
#         # \u00b0 : degree symbol
#         return '%f\u00b0%s, %f\u00b0%s' % (abs(lat), ns, abs(lon), ew)
#
#     class DegreeFormatter(Formatter):
#         """
#         This is a custom formatter that converts the native unit of
#         radians into (truncated) degrees and adds a degree symbol.
#         """
#
#         def __init__(self, round_to=1.0):
#             self._round_to = round_to
#
#         def __call__(self, x, pos=None):
#             degrees = round(np.degrees(x) / self._round_to) * self._round_to
#             # \u00b0 : degree symbol
#             return "%d\u00b0" % degrees
#
#     def set_longitude_grid(self, degrees):
#         """
#         Set the number of degrees between each longitude grid.
#
#         This is an example method that is specific to this projection
#         class -- it provides a more convenient interface to set the
#         ticking than set_xticks would.
#         """
#         # Set up a FixedLocator at each of the points, evenly spaced
#         # by degrees.
#         number = (360.0 / degrees) + 1
#         self.xaxis.set_major_locator(
#             plt.FixedLocator(
#                 np.linspace(-np.pi, np.pi, number, True)[1:-1]))
#         # Set the formatter to display the tick labels in degrees,
#         # rather than radians.
#         self.xaxis.set_major_formatter(self.DegreeFormatter(degrees))
#
#     def set_latitude_grid(self, degrees):
#         """
#         Set the number of degrees between each longitude grid.
#
#         This is an example method that is specific to this projection
#         class -- it provides a more convenient interface than
#         set_yticks would.
#         """
#         # Set up a FixedLocator at each of the points, evenly spaced
#         # by degrees.
#         number = (180.0 / degrees) + 1
#         self.yaxis.set_major_locator(
#             FixedLocator(
#                 np.linspace(-np.pi / 2.0, np.pi / 2.0, number, True)[1:-1]))
#         # Set the formatter to display the tick labels in degrees,
#         # rather than radians.
#         self.yaxis.set_major_formatter(self.DegreeFormatter(degrees))
#
#     def set_longitude_grid_ends(self, degrees):
#         """
#         Set the latitude(s) at which to stop drawing the longitude grids.
#
#         Often, in geographic projections, you wouldn't want to draw
#         longitude gridlines near the poles.  This allows the user to
#         specify the degree at which to stop drawing longitude grids.
#
#         This is an example method that is specific to this projection
#         class -- it provides an interface to something that has no
#         analogy in the base Axes class.
#         """
#         longitude_cap = np.radians(degrees)
#         # Change the xaxis gridlines transform so that it draws from
#         # -degrees to degrees, rather than -pi to pi.
#         self._xaxis_pretransform \
#             .clear() \
#             .scale(1.0, longitude_cap * 2.0) \
#             .translate(0.0, -longitude_cap)
#
#     def get_data_ratio(self):
#         """
#         Return the aspect ratio of the data itself.
#
#         This method should be overridden by any Axes that have a
#         fixed data ratio.
#         """
#         return 1.0
#
#     # Interactive panning and zooming is not supported with this projection,
#     # so we override all of the following methods to disable it.
#     def can_zoom(self):
#         """
#         Return True if this axes support the zoom box
#         """
#         return False
#
#     def start_pan(self, x, y, button):
#         pass
#
#     def end_pan(self):
#         pass
#
#     def drag_pan(self, button, key, x, y):
#         pass

    # Now, the transforms themselves.

    # class HammerTransform(Transform):
    #     """
    #     The base Hammer transform.
    #     """
    #     input_dims = 2
    #     output_dims = 2
    #     is_separable = False
    #
    #     def transform_non_affine(self, ll):
    #         """
    #         Override the transform_non_affine method to implement the custom
    #         transform.
    #
    #         The input and output are Nx2 numpy arrays.
    #         """
    #         longitude = ll[:, 0:1]
    #         latitude = ll[:, 1:2]
    #
    #         # Pre-compute some values
    #         half_long = longitude / 2.0
    #         cos_latitude = np.cos(latitude)
    #         sqrt2 = np.sqrt(2.0)
    #
    #         alpha = 1.0 + cos_latitude * np.cos(half_long)
    #         x = (2.0 * sqrt2) * (cos_latitude * np.sin(half_long)) / alpha
    #         y = (sqrt2 * np.sin(latitude)) / alpha
    #         return np.concatenate((x, y), 1)
    #
    #     # This is where things get interesting.  With this projection,
    #     # straight lines in data space become curves in display space.
    #     # This is done by interpolating new values between the input
    #     # values of the data.  Since ``transform`` must not return a
    #     # differently-sized array, any transform that requires
    #     # changing the length of the data array must happen within
    #     # ``transform_path``.
    #     def transform_path_non_affine(self, path):
    #         ipath = path.interpolated(path._interpolation_steps)
    #         return Path(self.transform(ipath.vertices), ipath.codes)
    #     transform_path_non_affine.__doc__ = \
    #         Transform.transform_path_non_affine.__doc__
    #
    #     if matplotlib.__version__ < '1.2':
    #         # Note: For compatibility with matplotlib v1.1 and older, you'll
    #         # need to explicitly implement a ``transform`` method as well.
    #         # Otherwise a ``NotImplementedError`` will be raised. This isn't
    #         # necessary for v1.2 and newer, however.
    #         transform = transform_non_affine
    #
    #         # Similarly, we need to explicitly override ``transform_path`` if
    #         # compatibility with older matplotlib versions is needed. With v1.2
    #         # and newer, only overriding the ``transform_path_non_affine``
    #         # method is sufficient.
    #         transform_path = transform_path_non_affine
    #         transform_path.__doc__ = Transform.transform_path.__doc__
    #
    #     def inverted(self):
    #         return HammerAxes.InvertedHammerTransform()
    #     inverted.__doc__ = Transform.inverted.__doc__
    #
    # class InvertedHammerTransform(Transform):
    #     input_dims = 2
    #     output_dims = 2
    #     is_separable = False
    #
    #     def transform_non_affine(self, xy):
    #         x = xy[:, 0:1]
    #         y = xy[:, 1:2]
    #
    #         quarter_x = 0.25 * x
    #         half_y = 0.5 * y
    #         z = np.sqrt(1.0 - quarter_x*quarter_x - half_y*half_y)
    #         longitude = 2*np.arctan((z*x)/(2.0*(2.0*z*z - 1.0)))
    #         latitude = np.arcsin(y*z)
    #         return np.concatenate((longitude, latitude), 1)
    #     transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__
    #
    #     # As before, we need to implement the "transform" method for
    #     # compatibility with matplotlib v1.1 and older.
    #     if matplotlib.__version__ < '1.2':
    #         transform = transform_non_affine
    #
    #     def inverted(self):
    #         # The inverse of the inverse is the original transform... ;)
    #         return HammerAxes.HammerTransform()
    #     inverted.__doc__ = Transform.inverted.__doc__
    #
    #
    #
