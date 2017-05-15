import datetime
import matplotlib.pyplot as plt

import ephem
import numpy as np
import numpy.random as npr

from functions import Time_dist
from functions import conversion
from functions import deconversion
from functions import factibles


class ACOSchedule(object):
    def __init__(self, X, observer, ND, T):
        # This is the inicialization function. It recieves an array X of points in (DEC,RA) coordinates (a two column numpy array),
        # an observer instance of the ephem library, the number of intervals ND for the discretization of the night, and the times
        # since the last visit in the array T with same length as X.

        print datetime.datetime.now()
        self.NightDisc = ND  # Number of intervals for discretization of the night.

        # ACS parameters---------
        self.q_0 = 0.2  # Parameter for step choosing in ACS
        self.chi = 0.1  # Parameter for pheromone wasting in ACS
        self.beta = 2  # Parameter of ACS
        self.rho = 0.1  # Parameter for pheromone updating in ACS
        self.NormObsQ = 1
        self.NormTimeQ = 1
        self.m = 60  # Number of ants per iteration.
        # ------------------------

        # Compute Night parameters
        sun = ephem.Sun()
        sun.compute(observer)
        self.NightEnd = observer.next_rising(sun)
        observer.date = self.NightEnd
        sun.compute(observer)
        self.NightBeg = observer.previous_setting(sun)
        observer.date = self.NightBeg
        self.NightLength = self.NightEnd - self.NightBeg
        self.Interval = self.NightLength / ND
        print self.NightLength
        # ----------------------

        self.X = X  # DEC-RA of discretization of the sky.
        self.NX = np.size(X, 0)
        self.obs = observer  # Observer instance of pyephem.
        self.T = T  # Time since last visit for X.
        self.Ph_ObsQ = self.init_Pheromone()  # Pheromone for quality of the observations (little zenith angle).
        self.Ph_TimeQ = self.init_Pheromone()  # Pheromone for correct distance in time between obs of the same point.
        self.Dist = self.DiscreteDistances()  # Distances between points
        self.ZenAn = self.DiscreteZenithAngle()
        self.AzAn = self.DiscreteAzimuth()
        self.Fact = self.DiscreteFactibles()
        [self.Times,self.MAA] = self.MoonPosition()
        self.Epsilon = 0.5
        if np.size(X,0) == 3072:
            self.ObservationTime = 4*30. / (24 * 3600)
        elif np.size(X,0) == 768:
            self.ObservationTime = 16*30. / (24 * 3600)
        else:
            self.ObservationTime = 30. / (24 * 3600)
        print self.ObservationTime

        # Best (in Pareto's sense) solutions
        self.BPS = []
        # This list stores the solutions. The solutions are also lists with the form [path,ObsQ,TimeQ].
        # ---------------------------------

        self.IterationNumber = -1
        self.ParetoHistorial = []

        print('Construccion Completa\n')
        print datetime.datetime.now()

    def RunACO_Pheromone(self, SuperAntIt, AntIterations):
        # Runs the algorithm. Recieves the number of SuperAnts, and the number of iterations of Ants.

        Lambda = 0.5  # Parameter in (0,1) for choosing the pheromone matrix to be used (used by Ants).

        print '*****************Super Ants******************'
        for i in range(SuperAntIt):
            print i, len(self.BPS), datetime.datetime.now()
            [Path, ObsQ, TimeQ] = self.SuperAnts()
            # self.EvaporatePh(0.99,0.99)
            self.Update_BPS(Path, ObsQ, TimeQ)
        self.Update_Pheromone()

        print '*****************Ants******************'
        self.AntIterations = AntIterations
        for j in range(AntIterations):
            print j, len(self.BPS), datetime.datetime.now()
            self.IterationNumber = j
            self.Colony(self.m, Lambda)
        return

    def Colony(self, m, Lambda):
        # Runs an iteration of the ACS algorithm. m is the number of Ants and Lambda is for choosing the matrix of pheromones.
        """Sum = np.cumsum(self.Ph_ObsQ[0][1400, :] * self.Fact[0])
        eta = 1./(self.Dist[0][1400,:]+np.min(filter(lambda x:x>0,self.Dist[0][1400,:]))/10.)
        #print Sum[-1]
        maxS = Sum[-1]
        maxeta = np.max(eta)
        mineta = np.min(eta)
        for i in range(1450,1470):
            eta = 1./(self.Dist[0][i,:]+np.min(filter(lambda x:x>0,self.Dist[0][i,:]))/10.)
            Ceta = np.cumsum(eta)
            Ceta2 = np.cumsum(eta*self.Fact[0])
            S = np.cumsum(self.Ph_ObsQ[0][i, :] * self.Fact[0])
            Seta = np.cumsum(self.Ph_ObsQ[0][i, :] * self.Fact[0] * np.power(eta,self.beta))
            maxS = max(S[-1],maxS)
            maxeta = max(maxeta,np.max(eta))
            mineta = min(mineta,np,min(eta))
            plt.figure(1)
            plt.plot(S/S[-1])
            plt.figure(2)
            plt.plot(Ceta2/Ceta2[-1])
            plt.figure(3)
            plt.plot(Seta/Seta[-1])
            plt.figure(4)
            plt.plot(Ceta2*S/(Ceta2[-1]*S[-1]))
        print maxS
        print mineta,maxeta
        plt.show()"""

        """etaZ1 = 1. / self.ZenAn[0]
        etaZ1 = etaZ1 * self.Fact[0]
        etaZ1 = etaZ1 / sum(etaZ1)  # print Sum[-1]
        etaZ2 = 3 * np.pi / (2 * (self.ZenAn[0] + np.pi))
        etaZ2 = etaZ2 * self.Fact[0]
        etaZ2 = etaZ2 / sum(etaZ2)
        plt.figure(1)
        plt.plot(np.cumsum(etaZ1))
        plt.figure(2)
        plt.plot(np.cumsum(etaZ2))
        plt.show()"""

        """etaT1 = self.T
        etaT1 = etaT1 * self.Fact[0]
        etaT1 = 1.0 * etaT1 / sum(etaT1)
        etaT2 = np.power(2,self.T)
        etaT2 = etaT2 * self.Fact[0]
        etaT2 = 1.0 * etaT2 / sum(etaT2)
        etaT3 = (np.power(2,self.T)+128.)/192.0
        etaT3 = etaT3 * self.Fact[0]
        etaT3 = 1.0 * etaT3 / sum(etaT3)
        plt.figure(1)
        plt.plot(np.cumsum(etaT1))
        plt.figure(2)
        plt.plot(np.cumsum(etaT2))
        plt.figure(3)
        plt.plot(np.cumsum(etaT3))
        plt.show()"""

        F = np.transpose(np.ones(np.size(self.X,0)).reshape(np.size(self.X,0),1)*self.Fact[0])

        A = np.zeros((np.size(self.X,0),np.size(self.X,0)))
        minA = []
        maxA = []
        for i in range(self.NightDisc):
            F = np.transpose(np.ones(np.size(self.X, 0)).reshape(np.size(self.X, 0), 1) * self.Fact[i])
            Ft = np.transpose(F)
            A += Ft*self.Ph_ObsQ[i]*F
            minA.append(np.min(self.Ph_ObsQ[i]))
            maxA.append(np.max(self.Ph_ObsQ[i]))
        plt.matshow(A)
        plt.show()
        #print minA
        #print maxA

        A = np.zeros((np.size(self.X,0),np.size(self.X,0)))
        minA = []
        maxA = []
        for i in range(self.NightDisc):
            F = np.transpose(np.ones(np.size(self.X, 0)).reshape(np.size(self.X, 0), 1) * self.Fact[i])
            Ft = np.transpose(F)
            A += Ft*self.Ph_TimeQ[i]*F
            minA.append(np.min(self.Ph_TimeQ[i]))
            maxA.append(np.max(self.Ph_TimeQ[i]))
        plt.matshow(A)
        plt.show()
        #print minA
        #print maxA


        self.iterationBPS = []
        for j in range(m):
            Lambda = j / (self.m - 1)
            [Path, ObsQ, TimeQ] = self.Ants(Lambda)
            self.Update_iterationBPS(Path, ObsQ, TimeQ)

        for j in range(len(self.iterationBPS)):
            B_aux = self.iterationBPS[j]
            self.Update_BPS(B_aux[0], B_aux[1], B_aux[2])

        self.Update_Pheromone_Iteration()
        #print np.min(self.Ph_ObsQ),np.max(self.Ph_ObsQ),np.min(self.Ph_TimeQ),np.max(self.Ph_TimeQ)
        return

    def Ants(self, Lambda):
        # This function runs one ant for the ACS algorithm.

        TotalT = 0
        Path = []
        NotVisited = np.ones(np.size(self.Fact[0]))
        (Fact0no0,) = np.nonzero(self.Fact[0])
        Path.append([Fact0no0[npr.randint(0, Fact0no0.size)], 0])
        ActualPoint = Path[0][0]
        ActualPeriod = 0
        # PhRowAux=np.zeros(np.size(self.X,0))


        while TotalT <= self.NightLength:
            U = npr.rand()
            #eta=np.power(((np.power(2,self.T[ActualPeriod])+128)/192.)/((self.Dist[ActualPeriod][ActualPoint,:]+0.0001)+(self.ZenAn[ActualPeriod]+np.pi)/(3*np.pi/2)),1)
            #eta = np.power(((np.power(2, self.T[ActualPeriod]) + 128) / 192.) / (self.Dist[ActualPeriod][ActualPoint, :] + 0.0001),self.beta)
            #eta = 1./(self.Dist[ActualPeriod][ActualPoint,:]+np.min(filter(lambda x:x>0,self.Dist[ActualPeriod][ActualPoint,:]))/10.)
            #eta = np.divide(self.ZenAn[ActualPeriod],self.Dist[ActualPeriod][ActualPoint,:])
            #eta = self.T/(self.Dist[ActualPeriod][ActualPoint,:]*self.ZenAn[ActualPeriod])
            #eta = np.ones(np.size(self.Dist[ActualPeriod][ActualPoint,:]))

            etaD = 1. / (self.Dist[ActualPeriod][ActualPoint, :] + np.min(
                filter(lambda x: x > 0, self.Dist[ActualPeriod][ActualPoint, :])) / 10.)
            etaD = etaD * self.Fact[ActualPeriod] * NotVisited
            etaD = etaD / np.sum(etaD)

            etaZ = 1. / self.ZenAn[ActualPeriod]
            #etaZ = 3*np.pi/(2*(self.ZenAn[ActualPeriod]+np.pi))
            etaZ = etaZ * self.Fact[ActualPeriod] * NotVisited
            etaZ = etaZ / sum(etaZ)

            #etaT = self.T
            etaT = np.power(2,self.T)
            #etaT = (np.power(2,self.T)+128.)/192.
            etaT = 1.0 * etaT * self.Fact[ActualPeriod] * NotVisited
            etaT = etaT / sum(etaT)

            eta = np.power(etaD,3) * np.power(etaZ,2) * np.power(etaT,1)

            if U <= Lambda:
                tau = self.Ph_ObsQ[ActualPeriod][ActualPoint, :] * self.Fact[ActualPeriod] * NotVisited
                tau = tau / np.sum(tau)
                index = self.ChooseStep(tau * eta)

            else:
                tau = self.Ph_TimeQ[ActualPeriod][ActualPoint, :] * self.Fact[ActualPeriod] * NotVisited
                tau = tau / np.sum(tau)
                index = self.ChooseStep(tau * eta)

            TotalT += self.Dist[ActualPeriod][ActualPoint, index] + self.ObservationTime
            NewPeriod = min(int(np.floor(TotalT / self.Interval)), self.NightDisc - 1)

            self.Ph_ObsQ[ActualPeriod][ActualPoint, index] *= (1 - self.chi)
            self.Ph_ObsQ[ActualPeriod][ActualPoint, index] += self.chi
            self.Ph_TimeQ[ActualPeriod][ActualPoint, index] *= (1 - self.chi)
            self.Ph_TimeQ[ActualPeriod][ActualPoint, index] += self.chi

            ActualPoint = index
            ActualPeriod = NewPeriod
            Path.append([ActualPoint, ActualPeriod])
            NotVisited[ActualPoint] = 0
            if sum(NotVisited * self.Fact[ActualPeriod]) == 0:
                print 'no more factible points'
                NotVisited += 1
                break

        [ObsQ, TimeQ] = self.ObjectiveValues(np.array(Path, dtype=int))

        return [np.array(Path, dtype=int), ObsQ, TimeQ]

    def ChooseStep(self, values):
        # This function chooses the next step of an ant, using values as the probabilities p_ij.
        U = npr.rand()
        if U <= self.q_0:
            index = np.argmax(values)
        else:
            SumValues = np.cumsum(values)
            Sum = SumValues[-1]
            U2 = npr.rand() * Sum
            index = np.searchsorted(SumValues, U2)
        return index

    def SuperAnts(self):
        # This function computes a (good) factible path.

        TotalT = 0
        Path = []
        NotVisited = np.ones(np.size(self.Fact[0]))
        (Fact0no0,) = np.nonzero(self.Fact[0])
        Path.append([Fact0no0[npr.randint(0, Fact0no0.size)], 0])
        ActualPoint = Path[0][0]
        ActualPeriod = 0
        while TotalT <= self.NightLength:
            # while part_sum<U2:
            #	part_sum+=PhRowAux[index]
            Valor = (((self.T * self.Fact[ActualPeriod] >= 3.65) + 0.01) * self.T * self.Fact[
                ActualPeriod] * NotVisited) * 1000 / (
                    5 * np.abs(self.ZenAn[ActualPeriod]) + 5 * self.Dist[ActualPeriod][ActualPoint, :] + 1)
            RR = npr.randint(6)
            index = np.argpartition(-Valor, RR)[RR]
            #	index+=1
            TotalT += self.Dist[ActualPeriod][ActualPoint, index] + self.ObservationTime
            NewPeriod = min(int(np.floor(TotalT / self.Interval)), self.NightDisc - 1)
            ActualPoint = index
            ActualPeriod = NewPeriod
            Path.append([ActualPoint, ActualPeriod])

            NotVisited[ActualPoint] = 0
            if sum(NotVisited * self.Fact[ActualPeriod]) == 0:
                # print 'no more factible points'
                NotVisited += 1
                break

        [ObsQ, TimeQ] = self.ObjectiveValues(np.array(Path, dtype=int))

        return [np.array(Path, dtype=int), ObsQ, TimeQ]

    def Update_BPS(self, path, ObsQ, TimeQ):
        # Update the list of Pareto eficient paths, with path and it's objective values ObsQ and TimeQ (only if it is non dominated).
        if len(self.BPS) == 0:
            self.BPS.append([path, ObsQ, TimeQ])
            self.NormObsQ = ObsQ
            self.NormTimeQ = TimeQ
            self.ParetoHistorial.append([ObsQ, TimeQ, self.IterationNumber])
            print 'new non dominated solution', datetime.datetime.now()
        else:
            BPS_Aux = []
            flag_pareto = 1
            while 0 < len(self.BPS):
                A = self.BPS.pop()
                if (((ObsQ > A[1]) * (TimeQ >= A[2])) + ((ObsQ >= A[1]) * (TimeQ > A[2]))) == 0:
                    BPS_Aux.append(A)
                    if (ObsQ <= A[1]) * (TimeQ <= A[2]):
                        flag_pareto = 0
            if flag_pareto:
                BPS_Aux.append([path, ObsQ, TimeQ])
                self.ParetoHistorial.append([ObsQ, TimeQ, self.IterationNumber])
                print 'new non dominated solution', datetime.datetime.now()
            self.BPS = BPS_Aux
        return

    def Update_iterationBPS(self, path, ObsQ, TimeQ):
        # Update the list of Pareto eficient paths, with path and it's objective values ObsQ and TimeQ (only if it is non dominated).
        if len(self.iterationBPS) == 0:
            self.iterationBPS.append([path, ObsQ, TimeQ])
            self.NormObsQ = ObsQ
            self.NormTimeQ = TimeQ
        else:
            BPS_Aux = []
            flag_pareto = 1
            while 0 < len(self.iterationBPS):
                A = self.iterationBPS.pop()
                if (((ObsQ > A[1]) * (TimeQ >= A[2])) + ((ObsQ >= A[1]) * (TimeQ > A[2]))) == 0:
                    BPS_Aux.append(A)
                    if (ObsQ <= A[1]) * (TimeQ <= A[2]):
                        flag_pareto = 0
            if flag_pareto:
                BPS_Aux.append([path, ObsQ, TimeQ])
                #self.ParetoHistorial.append([ObsQ, TimeQ, self.IterationNumber])
                #print 'new non dominated solution', datetime.datetime.now()
            self.iterationBPS = BPS_Aux
        return

    def Update_Pheromone(self):
        for j in range(len(self.BPS)):
            addObsQ = self.BPS[j][1] / self.NormObsQ
            addTimeQ = self.BPS[j][2] / self.NormTimeQ
            path = self.BPS[j][0]
            for i in range(np.size(path, 0) - 1):
                self.Ph_ObsQ[path[i, 1]][path[i, 0], path[i + 1, 0]] *= (1 - self.rho)
                # self.Ph_ObsQ[path[i,1]][path[i,0],path[i+1,0]]+=self.rho*addObsQ
                self.Ph_ObsQ[path[i, 1]][path[i, 0], path[i + 1, 0]] += addObsQ
                self.Ph_TimeQ[path[i, 1]][path[i, 0], path[i + 1, 0]] *= (1 - self.rho)
                # self.Ph_TimeQ[path[i,1]][path[i,0],path[i+1,0]]+=self.rho*addTimeQ
                self.Ph_TimeQ[path[i, 1]][path[i, 0], path[i + 1, 0]] += addTimeQ
        return

    def Update_Pheromone_Iteration(self):
        for j in range(len(self.iterationBPS)):
            addObsQ = 10*self.iterationBPS[j][1] / self.NormObsQ
            addTimeQ = 10*self.iterationBPS[j][2] / self.NormTimeQ
            path = self.iterationBPS[j][0]
            for i in range(np.size(path, 0) - 1):
                self.Ph_ObsQ[path[i, 1]][path[i, 0], path[i + 1, 0]] *= (1 - self.rho)
                self.Ph_ObsQ[path[i, 1]][path[i, 0], path[i + 1, 0]] += addObsQ
                self.Ph_TimeQ[path[i, 1]][path[i, 0], path[i + 1, 0]] *= (1 - self.rho)
                self.Ph_TimeQ[path[i, 1]][path[i, 0], path[i + 1, 0]] += addTimeQ
        return

    def ObjectiveValues(self, SCH):
        # Calculate the value of the objective functions for a factible schedule SCH.
        ObsQ = 0
        TimeQ = 0
        for i in range(np.size(SCH, 0)):
            ObsQ += np.pi / 2 - self.ZenAn[SCH[i, 1]][SCH[i, 0]]
            TimeQ += self.TimeFunc(self.T[SCH[i, 0]])
        return [ObsQ, TimeQ]

    def DiscreteDistances(self):
        # Calculate the distances for the discretization of time.
        D = []
        NX = np.size(self.X, 0)
        for i in range(self.NightDisc):
            self.obs.date = self.NightBeg + (i + 0.5) * self.Interval
            D.append(np.zeros((NX, NX)))
            for j in range(NX):
                D[i][j, :] = Time_dist(self.X, self.obs, self.X[j, :])
        self.obs.date = self.NightBeg
        return D

    def DiscreteZenithAngle(self):
        # Calculate zenith angles for the discretization of time.
        ZA = []
        NX = np.size(self.X, 0)
        for i in range(self.NightDisc):
            self.obs.date = self.NightBeg + (i + 0.5) * self.Interval
            [Azi, Alt] = conversion(self.obs, self.X[:, 0], self.X[:, 1])
            ZA.append(np.pi / 2 - np.array(Alt))
        return ZA

    def DiscreteAzimuth(self):
        # Calculate Azimuth angles for the discretization of time.
        AZ = []
        NX = np.size(self.X, 0)
        for i in range(self.NightDisc):
            self.obs.date = self.NightBeg + (i + 0.5) * self.Interval
            [Azi, Alt] = conversion(self.obs, self.X[:, 0], self.X[:, 1])
            AZ.append(Azi)
        return AZ

    def DiscreteFactibles(self):
        # Calculate factible points for the discretization of time.
        Fact = []
        for i in range(self.NightDisc):
            self.obs.date = self.NightBeg + (i + 0.5) * self.Interval
            Fact.append(factibles(self.X, self.obs, 0))
        return Fact

    def set_Ph_ObsQ(self, M):
        self.Ph_ObsQ = M
        return

    def set_Ph_TimeQ(self, M):
        self.Ph_TimeQ = M
        return

    def init_Pheromone(self):
        P = []
        for i in range(self.NightDisc):
            P.append(np.ones((self.NX, self.NX)))
        return P

    def TimeFunc(self, t):
        # This is the function that is evaluated when TimeQ is calculated. TimeQ is the sum of TimeFunc(time since last visit)
        # over each observation spot of the schedule. The unit is [days].
        return np.power(2, t)

    def EvaporatePh(self, alphaObsQ, alphaTimeQ):
        #Obsolete
        print 'evaporate'
        # Recieves alpha in (0,1) and evaporates pheromones multiplying by alpha the matrices.
        for i in range(self.NightDisc):
            self.Ph_ObsQ[i] *= alphaObsQ
            self.Ph_TimeQ[i] *= alphaTimeQ
        return

    def AZALT(self, Path):
        # Calculates the (AZ,ALT) coordinates for the points of a Path at the times at which they are visited.
        AA = []
        DR = []
        for i in range(np.size(Path, 0)):
            aa = [self.AzAn[Path[i, 1]][Path[i, 0]], np.pi / 2 - self.ZenAn[Path[i, 1]][Path[i, 0]]]
            AA.append(aa)

            self.obs.date = self.NightBeg + (Path[i, 1] + 0.5) * self.Interval
            dr = deconversion(self.obs, aa[0], aa[1])
            DR.append(dr)
        return [np.array(AA), np.array(DR)]

    def MoonPosition(self):
        MAA = []
        Times = []
        moon = ephem.Moon()
        for i in range(self.NightDisc):
            self.obs.date = self.NightBeg + (i + 0.5) * self.Interval
            print self.obs.date
            moon.compute(self.obs)
            Times.append(self.obs.date)
            MAA.append([moon.ra, moon.dec])
        return [np.array(Times), np.array(MAA)]

    def PlotParetoHistorial(self):
        # Plots the objective values of the non dominated solutions found during the whole process, even if they were eliminated.
        # It shows the evolution of the solutions encountered.
        PHist = np.array(self.ParetoHistorial)
        plt.scatter(PHist[:, 0], PHist[:, 1], c=np.arange(np.size(PHist, 0)), alpha=0.7, s=100)
        plt.show()

    def PlotParetoFront(self):
        # Plots the the non dominated solutions found in the end
        PFX = [self.BPS[i][1] for i in range(len(self.BPS))]
        PFY = [self.BPS[i][2] for i in range(len(self.BPS))]
        plt.scatter(PFX, PFY)
        plt.show()