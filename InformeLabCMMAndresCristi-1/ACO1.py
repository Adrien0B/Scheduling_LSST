import ephem
import numpy as np
import numpy.random as npr
import healpy as hpy
import datetime
import matplotlib.pyplot as plt
import pickle
from functions import conversion
from functions import factibles
from functions import Time_dist 


class ACOSchedule(object):
	def __init__(self,X,observer,ND,T):
		#This is the inicialization function. It recieves an array X of points in (DEC,RA) coordinates (a two column numpy array),
		# an observer instance of the ephem library, the number of intervals ND for the discretization of the night, and the times
		# since the last visit in the array T with same length as X.

		print datetime.datetime.now()
		self.NightDisc=ND 	#Number of intervals for discretization of the night.

		# ACS parameters---------
		self.q_0=0.9		#Parameter for step choosing in ACS
		self.chi=0.1		#Parameter for pheromone wasting in ACS
		self.beta=2		#Parameter of ACS
		self.rho=1		#Parameter for pheromone updating in ACS
		self.NormObsQ=1	
		self.NormTimeQ=1
		self.m=10		#Number of ants per iteration.
		#------------------------

		#Compute Night parameters
		sun=ephem.Sun()
		sun.compute(observer)
		self.NightEnd=observer.next_rising(sun)
		observer.date=self.NightEnd
		sun.compute(observer)
		self.NightBeg=observer.previous_setting(sun)
		observer.date=self.NightBeg
		self.NightLength=self.NightEnd-self.NightBeg
		self.Interval=self.NightLength/ND
		print self.NightLength
		#----------------------

		self.X=X					#DEC-RA of discretization of the sky.
		self.NX=np.size(X,0)
		self.obs=observer				#Observer instance of pyephem.
		self.T=T					#Time since last visit for X.
		self.Ph_ObsQ=self.init_Pheromone()		#Pheromone for quality of the observations (little zenith angle).
		self.Ph_TimeQ=self.init_Pheromone()		#Pheromone for correct distance in time between obs of the same point.
		self.Dist=self.DiscreteDistances()		#Distances between points
		self.ZenAn=self.DiscreteZenithAngle()
		self.AzAn=self.DiscreteAzimuth()
		self.Fact=self.DiscreteFactibles()
		self.Epsilon=0.5
		self.ObservationTime=30./(24*3600)
		print self.ObservationTime

		#Best (in Pareto's sense) solutions
		self.BPS= []	
			#This list stores the solutions. The solutions are also lists with the form [path,ObsQ,TimeQ].
		#---------------------------------

		self.IterationNumber=-1
		self.ParetoHistorial=[]

		print('Construccion Completa\n')
		print datetime.datetime.now()


	def RunACO_Pheromone(self,SuperAntIt,AntIterations):
		# Runs the algorithm. Recieves the number of SuperAnts, and the number of iterations of Ants.
		
		Lambda=0.5		#Parameter in (0,1) for choosing the pheromone matrix to be used (used by Ants).

		print '*****************Super Ants******************'
		for i in range(SuperAntIt):
			print i,len(self.BPS), datetime.datetime.now()
			[Path,ObsQ,TimeQ]=self.SuperAnts()
			#self.EvaporatePh(0.99,0.99)
			self.Update_BPS(Path,ObsQ,TimeQ)
		self.Update_Pheromone()

		print '*****************Ants******************'
		for j in range(AntIterations):
			print j,len(self.BPS), datetime.datetime.now()
			self.IterationNumber=j
			self.Colony(self.m,Lambda)
		return


	def Colony(self,m,Lambda):
		#Runs an iteration of the ACS algorithm. m is the number of Ants and Lambda is for choosing the matrix of pheromones.
		for j in range(m):
			[Path,ObsQ,TimeQ]=self.Ants(Lambda)
			self.Update_BPS(Path,ObsQ,TimeQ)
		self.Update_Pheromone()
		return


	def Ants(self,Lambda):
		#This function runs one ant for the ACS algorithm.

		TotalT=0
		Path=[]
		NotVisited=np.ones(np.size(self.Fact[0]))
		(Fact0no0,)=np.nonzero(self.Fact[0])
		Path.append([Fact0no0[npr.randint(0,Fact0no0.size)],0])
		ActualPoint=Path[0][0]
		ActualPeriod=0
		PhRowAux=np.zeros(np.size(self.X,0))


		while TotalT<=self.NightLength:
			U=npr.rand()
			#eta=np.power(((np.power(2,self.T[ActualPeriod])+128)/192.)/((self.Dist[ActualPeriod][ActualPoint,:]+0.0001)+(self.ZenAn[ActualPeriod]+np.pi)/(3*np.pi/2)),1)
			eta=np.power(((np.power(2,self.T[ActualPeriod])+128)/192.)/(self.Dist[ActualPeriod][ActualPoint,:]+0.0001),self.beta)
			
			if U<=Lambda:
				index=self.ChooseStep(self.Ph_ObsQ[ActualPeriod][ActualPoint,:]*self.Fact[ActualPeriod]*NotVisited*eta)

			else:
				index=self.ChooseStep(self.Ph_TimeQ[ActualPeriod][ActualPoint,:]*self.Fact[ActualPeriod]*NotVisited*eta)

			TotalT+=self.Dist[ActualPeriod][ActualPoint,index]+self.ObservationTime
			NewPeriod=min(int(np.floor(TotalT/self.Interval)),self.NightDisc-1)

			self.Ph_ObsQ[ActualPeriod][ActualPoint,index]*=(1-self.chi)
			self.Ph_ObsQ[ActualPeriod][ActualPoint,index]+=self.chi
			self.Ph_TimeQ[ActualPeriod][ActualPoint,index]*=(1-self.chi)
			self.Ph_TimeQ[ActualPeriod][ActualPoint,index]+=self.chi

			ActualPoint=index
			ActualPeriod=NewPeriod
			Path.append([ActualPoint,ActualPeriod])
			NotVisited[ActualPoint]=0
			if sum(NotVisited*self.Fact[ActualPeriod])==0: 
				#print 'no more factible points'
				NotVisited+=1
				break

		[ObsQ,TimeQ]=self.ObjectiveValues(np.array(Path,dtype=int))

		return [np.array(Path,dtype=int),ObsQ,TimeQ]


	def ChooseStep(self,values):
		#This function chooses the next step of an ant, using values as the probabilities p_ij.
		U=npr.rand()
		if U<= self.q_0:
			index=np.argmax(values)
		else:
			SumValues=np.cumsum(values)
			Sum=SumValues[-1]
			U2=npr.rand()*Sum
			index=np.searchsorted(SumValues,U2)
		return index

	def SuperAnts(self):
		#This function computes a (good) factible path.

		TotalT=0
		Path=[]
		NotVisited=np.ones(np.size(self.Fact[0]))
		(Fact0no0,)=np.nonzero(self.Fact[0])
		Path.append([Fact0no0[npr.randint(0,Fact0no0.size)],0])
		ActualPoint=Path[0][0]
		ActualPeriod=0
		while TotalT<=self.NightLength:
			#while part_sum<U2:
			#	part_sum+=PhRowAux[index]
			Valor=(((self.T*self.Fact[ActualPeriod]>=3.65)+0.01)*self.T*self.Fact[ActualPeriod]*NotVisited)*1000/(5*np.abs(self.ZenAn[ActualPeriod])+5*self.Dist[ActualPeriod][ActualPoint,:]+1)
			RR=npr.randint(6)
			index=np.argpartition(-Valor,RR)[RR]
			#	index+=1
			TotalT+=self.Dist[ActualPeriod][ActualPoint,index]+self.ObservationTime
			NewPeriod=min(int(np.floor(TotalT/self.Interval)),self.NightDisc-1)
			ActualPoint=index
			ActualPeriod=NewPeriod
			Path.append([ActualPoint,ActualPeriod])
			
			NotVisited[ActualPoint]=0
			if sum(NotVisited*self.Fact[ActualPeriod])==0: 
				#print 'no more factible points'
				NotVisited+=1
				break

		[ObsQ,TimeQ]=self.ObjectiveValues(np.array(Path,dtype=int))

		return [np.array(Path,dtype=int),ObsQ,TimeQ]



	def Update_BPS(self,path,ObsQ,TimeQ):
		#Update the list of Pareto eficient paths, with path and it's objective values ObsQ and TimeQ (only if it is non dominated).
		if len(self.BPS)==0:
			self.BPS.append([path,ObsQ,TimeQ])
			self.NormObsQ=ObsQ
			self.NormTimeQ=TimeQ
		else:
			BPS_Aux=[]
			flag_pareto=1
			while 0<len(self.BPS):
				A=self.BPS.pop()
				if (((ObsQ>A[1])*(TimeQ>=A[2]))+((ObsQ>=A[1])*(TimeQ>A[2])))==0:
					BPS_Aux.append(A)
					if (ObsQ<=A[1])*(TimeQ<=A[2]):
						flag_pareto=0
			if flag_pareto:
				BPS_Aux.append([path,ObsQ,TimeQ])
				self.ParetoHistorial.append([ObsQ,TimeQ,self.IterationNumber])
				print 'new non dominated solution',datetime.datetime.now()
			self.BPS=BPS_Aux
		return

	def Update_Pheromone(self):
		for j in range(len(self.BPS)):
			addObsQ=self.BPS[j][1]/self.NormObsQ
			addTimeQ=self.BPS[j][2]/self.NormTimeQ
			path=self.BPS[j][0]
			for i in range(np.size(path,0)-1):
				self.Ph_ObsQ[path[i,1]][path[i,0],path[i+1,0]]*=(1-self.rho)
				self.Ph_ObsQ[path[i,1]][path[i,0],path[i+1,0]]+=self.rho*addObsQ
				self.Ph_TimeQ[path[i,1]][path[i,0],path[i+1,0]]*=(1-self.rho)
				self.Ph_TimeQ[path[i,1]][path[i,0],path[i+1,0]]+=self.rho*addTimeQ
		return


	def ObjectiveValues(self,SCH):
		#Calculate the value of the objective functions for a factible schedule SCH.
		ObsQ=0
		TimeQ=0
		for i in range(np.size(SCH,0)):
			ObsQ+=np.pi/2- self.ZenAn[SCH[i,1]][SCH[i,0]]
			TimeQ+= self.TimeFunc(self.T[SCH[i,0]])
		return [ObsQ,TimeQ]

	def DiscreteDistances(self):
		#Calculate the distances for the discretization of time.
		D=[]
		NX=np.size(self.X,0)
		for i in range(self.NightDisc):
			self.obs.date=self.NightBeg+(i+0.5)*self.Interval
			D.append(np.zeros((NX,NX)))
			for j in range(NX):
				D[i][j,:]=Time_dist(self.X,self.obs,self.X[j,:])
		self.obs.date=self.NightBeg
		return D


	def DiscreteZenithAngle(self):
		#Calculate zenith angles for the discretization of time.
		ZA=[]
		NX=np.size(self.X,0)
		for i in range(self.NightDisc):
			self.obs.date=self.NightBeg+(i+0.5)*self.Interval
			[Azi,Alt]=conversion(self.obs.lat,self.X[:,0],self.X[:,1],self.obs.sidereal_time())
			ZA.append(np.pi/2-np.array(Alt))
		return ZA

	def DiscreteAzimuth(self):
		#Calculate Azimuth angles for the discretization of time.
		AZ=[]
		NX=np.size(self.X,0)
		for i in range(self.NightDisc):
			self.obs.date=self.NightBeg+(i+0.5)*self.Interval
			[Azi,Alt]=conversion(self.obs.lat,self.X[:,0],self.X[:,1],self.obs.sidereal_time())
			AZ.append(Azi)
		return AZ

	def DiscreteFactibles(self):
		#Calculate factible points for the discretization of time.
		Fact=[]
		for i in range(self.NightDisc):
			self.obs.date=self.NightBeg+(i+0.5)*self.Interval
			Fact.append(factibles(self.X,self.obs,0))
		return Fact

	def set_Ph_ObsQ(self,M):
		self.Ph_ObsQ=M
		return
	def set_Ph_TimeQ(self,M):
		self.Ph_TimeQ=M
		return
	def init_Pheromone(self):
		NX=np.size(self.X,0)
		P=[]
		for i in range(self.NightDisc):
			P.append(np.ones((NX,NX)))
		return P

	def TimeFunc(self,t):
		#This is the function that is evaluated when TimeQ is calculated. TimeQ is the sum of TimeFunc(time since last visit)
		# over each observation spot of the schedule. The unit is [days].
		return np.power(2,t) 


	def EvaporatePh(self,alphaObsQ,alphaTimeQ):
		#Recieves alpha in (0,1) and evaporates pheromones multiplying by alpha the matrices.
		for i in range(self.NightDisc):
			self.Ph_ObsQ[i]*=alphaObsQ
			self.Ph_TimeQ[i]*=alphaTimeQ
		return


	def AZALT(self,Path):
		#Calculates the (AZ,ALT) coordinates for the points of a Path at the times at which they are visited.
		AA=[]
		for i in range(np.size(Path,0)):
			aa=[self.AzAn[Path[i,1]][Path[i,0]],np.pi/2-self.ZenAn[Path[i,1]][Path[i,0]]]
			AA.append(aa)
		return AA

	def PlotParetoHistorial(self):
		#Plots the objective values of the non dominated solutions found during the whole process, even if they were eliminated.
		# It shows the evolution of the solutions encountered.
		PHist=np.array(self.ParetoHistorial)
		plt.scatter(PHist[:,0],PHist[:,1],c=np.arange(np.size(PHist,0)),alpha=0.7,s=100)
		plt.show()





#Initializing the observer instance
obs=ephem.Observer()
obs.lat="-33:27:00"
obs.lon="-70:40:00"
obs.date="2015/03/12 23:00:00"
#-------------------



#Calculating the Healpix discretization
Nside=16
Npix=hpy.pixelfunc.nside2npix(Nside)
ipixmin=int(Npix/2)
X=np.transpose(np.array(hpy.pixelfunc.pix2ang(Nside,np.arange(ipixmin,Npix,1))))
X[:,0]=np.pi/2-X[:,0]
X[:,1]-=np.pi
Num=np.size(X,0)
#---------------------------------



#Times since last visit
#T = npr.randint(6, size=Num)+1
T=np.load('Times.npy')		#This only works for Nside=16.
#-------------------------


ACO=ACOSchedule(X,obs,10,T)
ACO.RunACO_Pheromone(10,1000)
print len(ACO.BPS)

schedAA=ACO.AZALT(ACO.BPS[0][0])
SS=np.size(schedAA,0)

colores=np.linspace(1,SS-1,SS-1)/SS
np.save('schedAA',schedAA)
np.save('colores',colores)
np.save('ParetoHistorial',np.array(ACO.ParetoHistorial))
ACO.PlotParetoHistorial()

with open('BPS','wb') as f:
	pickle.dump(ACO.BPS,f)
