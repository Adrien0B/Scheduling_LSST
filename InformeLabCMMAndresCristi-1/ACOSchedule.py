import datetime
import matplotlib.pyplot as plt
import multiprocessing as mp

import ephem
import numpy as np
import numpy.random as npr
import random as rd
import os
import time

from functions import *


class ACOSchedule(object):
    def __init__(self):
        return

    def initialise(self, X, observer,start_date,end_date, ND, T):
        # This is the inicialization function. It recieves an array X of points in (DEC,RA) coordinates (a two column numpy array),
        # an observer instance of the ephem library, the number of intervals ND for the discretization of the night, and the times
        # since the last visit in the array T with same length as X.

        print datetime.datetime.now()
        # self.NightDisc = ND  # Number of intervals for discretization of the night.

        # ACS parameters---------
        self.q_0 = 0.2  # Parameter for step choosing in ACS
        self.chi = 0.1  # Parameter for pheromone wasting in ACS
        self.rho = 0.1  # Parameter for pheromone updating in ACS
        self.NormObsQ = 1
        self.NormTimeQ = 1
        self.m = 10  # Number of ants per iteration.
        # self.Filters = 1    # Number of Filters used
        # ------------------------

        # Compute Night parameters
        self.ND = ND
        self.obs = observer  # Observer instance of pyephem.
        # sun = ephem.Sun()
        # sun.compute(observer)
        # self.NightEnd = observer.next_rising(sun)
        # observer.date = self.NightEnd
        # sun.compute(observer)
        # self.NightBeg = observer.previous_setting(sun)
        # observer.date = self.NightBeg
        # self.NightLength = self.NightEnd - self.NightBeg
        # self.Interval = self.NightLength / ND
        self.start_date = start_date
        self.end_date = end_date
        self.Times = self.CreateTimes(start_date,end_date)
        self.NightLength = sum(self.Times[:,1])
        print self.NightLength
        # ----------------------

        self.X = X  # DEC-RA of discretization of the sky.
        self.Fact = self.DiscreteFactibles()
        self.TotalFact = sum(self.Fact)
        self.TotalFact = self.TotalFact!=0
        self.X = self.X[self.TotalFact,:]
        self.Fact = self.DiscreteFactibles()
        self.X_obs = []
        for time in self.Times[:,0]:
            self.obs.date = time
            self.X_obs.append(np.transpose(np.array(conversion(self.obs, self.X[:, 0], self.X[:, 1]))))  # [AZ,ALT] of X for the observer.
        self.NX = np.size(self.X, 0)
        self.T = np.array(T[self.TotalFact],dtype = float)  # Time since last visit for X.
        #self.T = T # Time since last visit for X.
        self.Ph_ObsQ = self.init_Pheromone()  # Pheromone for quality of the observations (little zenith angle).
        self.Ph_TimeQ = self.init_Pheromone()  # Pheromone for correct distance in time between obs of the same point.
        self.Ph_ObsQBeg = self.Ph_ObsQ[1,:]
        self.Ph_TimeQBeg = self.Ph_TimeQ[1,:]
        self.Dist = self.DiscreteDistances()  # Distances between points
        self.ZenAn = [np.pi / 2 - i[:, 1] for i in self.X_obs]
        self.AzAn = [i[:, 0] for i in self.X_obs]
        #self.ZenAn = self.DiscreteZenithAngle()
        #self.AzAn = self.DiscreteAzimuth()
        self.MAA = self.MoonPosition()
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
        self.AntIterations = 0
        self.ChangedIterations = []
        self.ParetoHistorial = []


        print('Construction Complete\n')
        print datetime.datetime.now()

    def initialise_from_file(self,title):
        # Copy all data from a previous run
        ACO = np.load(title)[0]
        print datetime.datetime.now()

        # ACS parameters---------
        self.q_0 = ACO.q_0  # Parameter for step choosing in ACS
        self.chi = ACO.chi  # Parameter for pheromone wasting in ACS
        self.rho = ACO.rho  # Parameter for pheromone updating in ACS
        self.NormObsQ = ACO.NormObsQ
        self.NormTimeQ = ACO.NormTimeQ
        self.m = ACO.m  # Number of ants per iteration.
        # ------------------------

        #Night parameters
        self.ND = ACO.ND
        self.obs = ACO.obs  # Observer instance of pyephem.
        self.start_date = ACO.start_date
        self.end_date = ACO.end_date
        self.Times = ACO.Times
        self.NightLength = ACO.NightLength
        print self.NightLength
        # ----------------------

        self.X = ACO.X  # DEC-RA of discretization of the sky.
        self.Fact = ACO.Fact
        self.TotalFact = ACO.TotalFact
        self.X_obs = ACO.X_obs
        self.NX = ACO.NX
        self.T = ACO.T  # Time since last visit for X.
        self.Ph_ObsQ = ACO.Ph_ObsQ  # Pheromone for quality of the observations (little zenith angle).
        self.Ph_TimeQ = ACO.Ph_TimeQ  # Pheromone for correct distance in time between obs of the same point.
        self.Ph_ObsQBeg = ACO.Ph_ObsQBeg
        self.Ph_TimeQBeg = ACO.Ph_TimeQBeg
        self.Dist = ACO.Dist  # Distances between points
        self.ZenAn = ACO.ZenAn
        self.AzAn = ACO.AzAn
        self.MAA = ACO.MAA
        self.Epsilon = ACO.Epsilon
        self.ObservationTime = ACO.ObservationTime
        print self.ObservationTime

        # Best (in Pareto's sense) solutions
        self.BPS = ACO.BPS
        # This list stores the solutions. The solutions are also lists with the form [path,ObsQ,TimeQ].
        # ---------------------------------

        self.IterationNumber = ACO.IterationNumber
        self.AntIterations = ACO.AntIterations
        self.ChangedIterations = list(ACO.ChangedIterations)
        self.ParetoHistorial = ACO.ParetoHistorial

        print('Construction Complete\n')
        print datetime.datetime.now()

    def RunACO_Pheromone(self, SuperAntIt, AntIterations):
        # Runs the algorithm. Recieves the number of SuperAnts, and the number of iterations of Ants.

        # Lambda = 0.5  # Parameter in (0,1) for choosing the pheromone matrix to be used (used by Ants).

        print '*****************Super Ants******************'
        for i in range(SuperAntIt):
            print i, len(self.BPS), datetime.datetime.now()
            [Path, ObsQ, TimeQ] = self.SuperAnts()
            # self.EvaporatePh(0.99,0.99)
            self.Update_BPS(Path, ObsQ, TimeQ)
        self.Update_Pheromone()

        print '*****************Ants******************'
        for j in range(AntIterations):
            print j+self.AntIterations, len(self.BPS), datetime.datetime.now()
            self.IterationNumber = j+self.AntIterations
            self.Colony(self.m)
        self.AntIterations = self.AntIterations + AntIterations
        return

    def Colony(self, m):
        # Runs an iteration of the ACS algorithm. m is the number of Ants and Lambda is for choosing the matrix of pheromones.


        # """Show heuristic cumulative probabilities"""
        # fig = plt.figure()
        # for i in range(100):
        #     etaD = 1. / (self.Dist[0][i, :] + np.min(
        #         filter(lambda x: x > 0, self.Dist[0][i, :])) / 10.)
        #     etaD = etaD * self.Fact[0]
        #     etaD = etaD / np.sum(etaD)
        #
        #     etaZ = 1. / self.ZenAn[0]
        #     # etaZ = 3*np.pi/(2*(self.ZenAn[ActualPeriod]+np.pi))
        #     etaZ = etaZ * self.Fact[0]
        #     etaZ = etaZ / sum(etaZ)
        #
        #     # etaT = self.T
        #     etaT = np.power(2, self.T) - 1.0
        #     # etaT = (np.power(2,self.T)+128.)/192.
        #     etaT = 1.0 * etaT * self.Fact[0]
        #     etaT = etaT / sum(etaT)
        #
        #     eta = np.power(etaD, 2) * np.power(etaZ, 2) * np.power(etaT, 0)
        #     eta = eta / np.sum(eta)
        #     # eta_sum = np.cumsum(eta)/sum(eta)
        #
        #     tau = self.Ph_ObsQ[i, :] * self.Fact[0]
        #     tau = tau / np.sum(tau)
        #
        #     prob = eta*tau
        #     prob_sum = np.cumsum(prob) / sum(prob)
        #     plt.plot(prob_sum)
        # # plt.show()
        # fig.savefig("videos/Probability/Probability_{:0>5}".format(self.IterationNumber))
        # plt.close()
        # plt.clf()

        # """Show and Save Pheromone matrices"""
        # save_figure=0
        # show_figure=1
        # plt.clf()
        # plt.matshow(self.Ph_ObsQ)
        # fig = plt.gcf()
        # if(save_figure):
        #     fig.savefig("videos/Pheromones/Pheromone_Observation_"+str(self.IterationNumber))
        # if(show_figure):
        #     plt.show(block=False)
        #     print np.min(self.Ph_ObsQ), np.max(self.Ph_ObsQ)
        # plt.clf()
        # plt.matshow(self.Ph_TimeQ)
        # fig = plt.gcf()
        # if (save_figure):
        #     fig.savefig("videos/Pheromones/Pheromone_Time_"+str(self.IterationNumber))
        # if (show_figure):
        #     plt.show(block=False)
        #     print np.min(self.Ph_TimeQ), np.max(self.Ph_TimeQ)
        # plt.clf()


        self.iterationBPS = []

        # for j in range(m):
        #     # print "Ant " + str(j+1)
        #     Lambda = j / (self.m - 1)
        #     [Path, ObsQ, TimeQ] = self.Ants(Lambda)
        #     # print ObsQ,TimeQ
        #     self.Update_iterationBPS(Path, ObsQ, TimeQ)

        Pipes = [mp.Pipe() for j in range(self.m)]
        processes = [mp.Process(target=Ant_Pipe, args=(self, j / (self.m - 1.0), Pipes[j][1])) for j in range(self.m)]
        for p in processes:
            p.start()
        # print mp.active_children()
        results = []
        j = 0
        for p in processes:
            results.append(Pipes[j][0].recv())
            j += 1
            p.join()
        for [Path, ObsQ, TimeQ] in results:
            # print ObsQ,TimeQ
            self.Update_iterationBPS(Path, ObsQ, TimeQ)

        for j in range(len(self.iterationBPS)):
            B_aux = self.iterationBPS[j]
            Change = self.Update_BPS(B_aux[0], B_aux[1], B_aux[2])

        if(Change):
            self.ChangedIterations.append(self.IterationNumber)

        self.Update_Pheromone_Iteration()
        return

    def Ants(self, Lambda):
        # This function runs one ant for the ACS algorithm.

        TotalT = 0
        LocalT = 0
        T_ant = np.copy(self.T)
        Path = []

        U = npr.rand()
        if U <= Lambda:
            tau = self.Ph_ObsQBeg * self.Fact[0]
            tau = tau / np.sum(tau)
            index = self.ChooseStep(tau)
        else:
            tau = self.Ph_TimeQBeg * self.Fact[0]
            tau = tau / np.sum(tau)
            index = self.ChooseStep(tau)

        Path.append([index, 0])
        TotalT += self.ObservationTime
        ActualPoint = index
        ActualPeriod = 0
        T_ant[index] = 0.0
        TimeQ = 0.0
        ObsQ = 0.0

        flagOverNight = False

        while TotalT <= self.NightLength:
            U = npr.rand()

            etaD = 1. / (self.Dist[ActualPeriod][ActualPoint, :] + np.min(
                filter(lambda x: x > 0, self.Dist[ActualPeriod][ActualPoint, :])) / 10.)
            etaD = etaD * self.Fact[ActualPeriod]
            etaD = etaD / np.sum(etaD)

            etaZ = 1. / self.ZenAn[ActualPeriod]
            #etaZ = 3*np.pi/(2*(self.ZenAn[ActualPeriod]+np.pi))
            etaZ = etaZ * self.Fact[ActualPeriod]
            etaZ = etaZ / sum(etaZ)

            #etaT = T_ant
            etaT = np.power(2,T_ant)-1.0
            #etaT = (np.power(2,T_ant)+128.)/192.
            etaT = 1.0 * etaT * self.Fact[ActualPeriod]
            etaT = etaT / sum(etaT)

            eta = np.power(etaD, 2) * np.power(etaZ, 2) * np.power(etaT, 1)
            #eta = etaD * etaZ * etaT

            if U <= Lambda:
                tau = self.Ph_ObsQ[ActualPoint, :] * self.Fact[ActualPeriod]
                tau = tau / np.sum(tau)
                index = self.ChooseStep(tau * eta)

            else:
                tau = self.Ph_TimeQ[ActualPoint, :] * self.Fact[ActualPeriod]
                tau = tau / np.sum(tau)
                index = self.ChooseStep(tau * eta)

            TotalT += self.Dist[ActualPeriod][ActualPoint, index] + self.ObservationTime
            LocalT += self.Dist[ActualPeriod][ActualPoint, index] + self.ObservationTime
            if(LocalT >= self.Times[ActualPeriod,1]):
                LocalT -= self.Times[ActualPeriod,1]
                NewPeriod = ActualPeriod+1
                T_ant += self.Times[ActualPeriod,1]
                if(NewPeriod%self.ND == 0):
                    flagOverNight = True
                if(NewPeriod >= len(self.Times)):
                    break
            else :
                NewPeriod = ActualPeriod

            if(not(flagOverNight)):
                ObsQ += np.pi / 2 - self.ZenAn[NewPeriod][index]        #Observation Quality based on zenith angle
                # ObsQ += 1/np.cos(self.ZenAn[NewPeriod][index])        #Observation quality based on air mass
                TimeQ += self.TimeFunc(T_ant[index])
                T_ant[index] = 0.0
                # self.Ph_ObsQ[ActualPoint, index] *= (1 - self.chi)
                # self.Ph_ObsQ[ActualPoint, index] += self.chi
                # self.Ph_TimeQ[ActualPoint, index] *= (1 - self.chi)
                # self.Ph_TimeQ[ActualPoint, index] += self.chi

                ActualPoint = index
                ActualPeriod = NewPeriod
                Path.append([ActualPoint, ActualPeriod])
            else:
                flagOverNight = False
                TotalT = sum(self.Times[0:NewPeriod,1])
                LocalT = 0.0
                ActualPeriod = NewPeriod


        return [np.array(Path, dtype=int), ObsQ, TimeQ]

    def ChooseStep(self, values):
        # This function chooses the next step of an ant, using values as the probabilities p_ij.
        # plt.figure()
        # plt.plot(np.cumsum(values)/sum(values))
        # plt.show()
        # U = npr.rand()
        U = rd.random()
        if U <= self.q_0:
            index = np.argmax(values)
        else:
            SumValues = np.cumsum(values)
            Sum = SumValues[-1]
            # U2 = npr.rand() * Sum
            U2 = rd.random() * Sum
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
        Change = False
        if len(self.BPS) == 0:
            self.BPS.append([path, ObsQ, TimeQ])
            self.ParetoHistorial.append([ObsQ, TimeQ, self.IterationNumber])
            Change = True
            print 'new non dominated solution', datetime.datetime.now(), 'Obs = ', ObsQ, 'Time = ',TimeQ
        else:
            BPS_Aux = []
            flag_pareto = 1
            while 0 < len(self.BPS):
                A = self.BPS.pop()
                # if (((ObsQ > A[1]) * (TimeQ >= A[2])) + ((ObsQ >= A[1]) * (TimeQ > A[2]))) == 0:            #Test with zenith angle quality
                if (((ObsQ < A[1]) * (TimeQ >= A[2])) + ((ObsQ <= A[1]) * (TimeQ > A[2]))) == 0:        #Test with air mass quality
                    BPS_Aux.append(A)
                    # if (ObsQ <= A[1]) * (TimeQ <= A[2]):            #Test with zenith angle quality
                    if (ObsQ >= A[1]) * (TimeQ <= A[2]):        #Test with air mass quality
                        flag_pareto = 0
            if flag_pareto:
                BPS_Aux.append([path, ObsQ, TimeQ])
                self.ParetoHistorial.append([ObsQ, TimeQ, self.IterationNumber])
                Change = True
                print 'new non dominated solution', datetime.datetime.now(), 'Obs = ', ObsQ, 'Time = ',TimeQ
            self.BPS = BPS_Aux
        return Change

    def Update_iterationBPS(self, path, ObsQ, TimeQ):
        # Update the list of Pareto eficient paths, with path and it's objective values ObsQ and TimeQ (only if it is non dominated).
        if len(self.iterationBPS) == 0:
            self.iterationBPS.append([path, ObsQ, TimeQ])
            self.NormObsQ = ObsQ
            self.NormTimeQ = TimeQ
        else:
            # self.NormObsQ = max(self.NormObsQ,ObsQ)         #Norm for zenith angle quality
            self.NormObsQ = min(self.NormObsQ,ObsQ)          #Norm for air mass quality
            self.NormTimeQ = max(self.NormTimeQ,TimeQ)
            BPS_Aux = []
            flag_pareto = 1
            while 0 < len(self.iterationBPS):
                A = self.iterationBPS.pop()
                # if (((ObsQ > A[1]) * (TimeQ >= A[2])) + ((ObsQ >= A[1]) * (TimeQ > A[2]))) == 0:            #Test with zenith angle quality
                if (((ObsQ < A[1]) * (TimeQ >= A[2])) + ((ObsQ <= A[1]) * (TimeQ > A[2]))) == 0:        #Test with air mass quality
                    BPS_Aux.append(A)
                    # if (ObsQ <= A[1]) * (TimeQ <= A[2]):            #Test with zenith angle quality
                    if (ObsQ >= A[1]) * (TimeQ <= A[2]):        #Test with air mass quality
                        flag_pareto = 0
            if flag_pareto:
                BPS_Aux.append([path, ObsQ, TimeQ])
                #self.ParetoHistorial.append([ObsQ, TimeQ, self.IterationNumber])
                #print 'new non dominated solution', datetime.datetime.now()
            self.iterationBPS = BPS_Aux
        return

    def Update_Pheromone(self):
        for j in range(len(self.BPS)):
            addObsQ = 10*self.BPS[j][1] / self.NormObsQ
            addTimeQ = 10*self.BPS[j][2] / self.NormTimeQ
            path = self.BPS[j][0]
            for i in range(np.size(path, 0) - 1):
                self.Ph_ObsQ[path[i, 0], path[i + 1, 0]] *= (1 - self.rho)
                # self.Ph_ObsQ[path[i,0],path[i+1,0]]+=self.rho*addObsQ
                self.Ph_ObsQ[path[i, 0], path[i + 1, 0]] += addObsQ
                self.Ph_TimeQ[path[i, 0], path[i + 1, 0]] *= (1 - self.rho)
                # self.Ph_TimeQ[path[i,0],path[i+1,0]]+=self.rho*addTimeQ
                self.Ph_TimeQ[path[i, 0], path[i + 1, 0]] += addTimeQ
        return

    def Update_Pheromone_Iteration(self):
        for j in range(len(self.iterationBPS)):
            # addObsQ = 10*self.iterationBPS[j][1] / self.NormObsQ            #Pheromone with zenith angle quality
            addObsQ = 10*self.NormObsQ / self.iterationBPS[j][1]        #Pheromone with air mass quality
            addTimeQ = 10*self.iterationBPS[j][2] / self.NormTimeQ
            path = self.iterationBPS[j][0]
            self.Ph_ObsQBeg[path[0]] *= (1 - self.rho)
            self.Ph_ObsQBeg[path[0]] += addObsQ
            self.Ph_TimeQBeg[path[0]] *= (1 - self.rho)
            self.Ph_TimeQBeg[path[0]] += addTimeQ
            for i in range(np.size(path, 0) - 1):
                self.Ph_ObsQ[path[i, 0], path[i + 1, 0]] *= (1 - self.rho)
                self.Ph_ObsQ[path[i, 0], path[i + 1, 0]] += addObsQ
                self.Ph_TimeQ[path[i, 0], path[i + 1, 0]] *= (1 - self.rho)
                self.Ph_TimeQ[path[i, 0], path[i + 1, 0]] += addTimeQ
        return

    def DiscreteDistances(self):
        # Calculate the distances for the discretization of time.

        t_start = time.time()
        # D_p = []
        # NX = np.size(self.X, 0)
        # i = 0
        # results = []
        # pool = mp.Pool(processes=4)
        # for time_l in self.Times[:, 0]:
        #     self.obs.date = time_l
        #     obs = (self.obs.date, self.obs.epoch, float(self.obs.lat), float(self.obs.lon), self.obs.elevation, self.obs.temp,self.obs.pressure)
        #     D_p.append(np.zeros((NX, NX)))
        #
        #     # results.append(pool.imap(lambda j : Time_dist_map(self.X_obs[i], obs, self.X_obs[i],i,j), range(NX)))
        #     results.append([pool.apply_async(Time_dist_Multi, args=(self.X_obs[i], obs, self.X_obs[i][j, :],i, j)) for j in range(NX)])
        #     i += 1
        #
        # pool.close()
        # pool.join()
        # pool.terminate()
        # output = [[p.get() for p in r] for r in results]
        # for output_aux in output:
        #     for o in output_aux:
        #         D_p[o[0]][o[1], :] = o[2]
        # print "Time for parallel : " + str(time.time() - t_start)

        # t_start = time.time()
        # D_p = []
        # NX = np.size(self.X, 0)
        # i = 0
        # results = []
        # Pipes = []
        # processes = []
        # for time_l in self.Times[:, 0]:
        #     self.obs.date = time_l
        #     obs = (self.obs.date, self.obs.epoch, float(self.obs.lat), float(self.obs.lon), self.obs.elevation, self.obs.temp,self.obs.pressure)
        #     D_p.append(np.zeros((NX, NX)))
        #     for j in range(NX):
        #         Pipes.append(mp.Pipe())
        #         processes.append(mp.Process(target=Time_dist_pipe, args=(self.X_obs[i], obs, self.X_obs[i][j, :], i, j, Pipes[-1][0])))
        #         processes[-1].start()
        #     i += 1
        #
        # print "Starting Processes : ",len(processes)
        #
        # # j=0
        # # for p in processes:
        # #     if(j%100 == 0):
        # #         print j
        # #     j += 1
        # #     p.start()
        #
        # print "All Processes Started"
        #
        # j = 0
        # for p in processes:
        #     results.append(Pipes[j][1].recv())
        #     j += 1
        #     p.join()
        #
        #
        # for o in results:
        #     D_p[o[0]][o[1], :] = o[2]
        # print "Time for parallel : " + str(time.time() - t_start)


        t_start = time.time()

        D = []
        NX = np.size(self.X, 0)
        i = 0
        for time_l in self.Times[:, 0]:
            self.obs.date = time_l
            D.append(np.zeros((NX, NX)))
            for j in range(NX):
                D[i][j, :] = Time_dist(self.X_obs[i], self.obs, self.X_obs[i][j, :])
            i+=1
        print "Time for sequential : " + str(time.time()-t_start)
        # print np.array_equal(np.array(D),np.array(D_p))
        return D

    def DiscreteFactibles(self):
        # Calculate factible points for the discretization of time.
        Fact = []
        for time in self.Times[:,0]:
            self.obs.date = time
            Fact.append(factibles(self.X, self.obs, 0))
        return Fact

    def set_Ph_ObsQ(self, M):
        self.Ph_ObsQ = M
        return

    def set_Ph_TimeQ(self, M):
        self.Ph_TimeQ = M
        return

    def init_Pheromone(self):
        P = np.ones((self.NX, self.NX))
        return P

    def TimeFunc(self, t):
        # This is the function that is evaluated when TimeQ is calculated. TimeQ is the sum of TimeFunc(time since last visit)
        # over each observation spot of the schedule. The unit is [days].
        return np.power(2, t)

    def AZALT(self, Path):
        # Calculates the (AZ,ALT) coordinates for the points of a Path at the times at which they are visited.
        AA = []
        DR = []
        for i in range(np.size(Path, 0)):
            aa = [self.AzAn[Path[i, 1]][Path[i, 0]], np.pi / 2 - self.ZenAn[Path[i, 1]][Path[i, 0]]]
            AA.append(aa)

            dr = self.X[Path[i,0]]
            dr = [dr[1],dr[0]]
            DR.append(dr)
        return [np.array(AA), np.array(DR)]

    def MoonPosition(self):
        MAA = []
        moon = ephem.Moon()
        for time in self.Times[:,0]:
            self.obs.date = time
            moon.compute(self.obs)
            MAA.append([moon.ra, moon.dec])
        return np.array(MAA)

    def CreateTimes(self,start_date,end_date):
        self.obs.date = end_date
        end_date = self.obs.date
        sun = ephem.Sun()
        self.obs.date = start_date
        sun.compute(self.obs)
        NightEnd = self.obs.next_rising(sun)
        self.obs.date = NightEnd
        sun.compute(self.obs)
        NightBeg = self.obs.previous_setting(sun)
        Times = []

        while(NightBeg < end_date):
            self.obs.date = NightBeg
            ThisNightLength = NightEnd - NightBeg
            Interval = ThisNightLength / self.ND



            for i in range(self.ND):
                self.obs.date = NightBeg + (i + 0.5) * Interval
                print self.obs.date
                Times.append([self.obs.date,Interval])
            self.obs.date = NightEnd
            NightEnd = self.obs.next_rising(sun)
            self.obs.date = NightEnd
            sun.compute(self.obs)
            NightBeg = self.obs.previous_setting(sun)

        return np.array(Times)

    def PlotParetoHistorial(self,title="",show=False):
        # Plots the objective values of the non dominated solutions found during the whole process, even if they were eliminated.
        # It shows the evolution of the solutions encountered.
        PHist = np.array(self.ParetoHistorial)
        fig = plt.figure()
        plt.scatter(PHist[:, 0], PHist[:, 1], c=np.arange(np.size(PHist, 0)), alpha=0.7, s=100)
        plt.title(title)
        if(show):
            plt.show(block=False)
        return fig

    def PlotParetoFront(self,title="",show=False):
        # Plots the the non dominated solutions found in the end
        PFX = [self.BPS[i][1] for i in range(len(self.BPS))]
        PFY = [self.BPS[i][2] for i in range(len(self.BPS))]
        fig = plt.figure()
        plt.scatter(PFX, PFY)
        plt.title(title)
        if(show):
            plt.show(block=False)
        return fig

    def set_time(self):
        self.timenow = datetime.datetime.now()

    def saveACO(self,title=""):
        if(title == ""):
            title = "videos/%s-%s-%s_%s-%s-%s_ACO_Save" % (self.timenow.year,self.timenow.month,self.timenow.day,self.timenow.hour,self.timenow.minute,self.timenow.second)
        if(title[-4:]!=".npy"):
            title = title + ".npy"
        np.save(title,[self])
        return title


# def Ant_Multi(ACO,Lambda):
#     return Ants_ext(ACO,Lambda)
#
# def Ant_Pool(ACO,Lambda,output):
#     # print "Ant %d : Begining" % (Lambda*(ACO.m-1),)
#     # print os.getppid()
#     output.put(Ants_ext(ACO,Lambda))
#     # print "Ant %d : End" % (Lambda*(ACO.m-1),)

def Ants_ext(ACO, Lambda):
    # This function runs one ant for the ACS algorithm.
    rd.seed(os.urandom(624))

    TotalT = 0
    LocalT = 0
    T_ant = np.copy(ACO.T)
    Path = []

    U = rd.random()
    if U <= Lambda:
        tau = ACO.Ph_ObsQBeg * ACO.Fact[0]
        tau = tau / np.sum(tau)
        index = ACO.ChooseStep(tau)
    else:
        tau = ACO.Ph_TimeQBeg * ACO.Fact[0]
        tau = tau / np.sum(tau)
        index = ACO.ChooseStep(tau)

    Path.append([index, 0])
    ActualPoint = index
    ActualPeriod = 0
    T_ant[index] = 0.0
    TimeQ = 0.0
    ObsQ = 0.0

    flagOverNight = False

    while TotalT <= ACO.NightLength:
        # U = rd.random()

        etaD = 1. / (ACO.Dist[ActualPeriod][ActualPoint, :] + np.min(
            filter(lambda x: x > 0, ACO.Dist[ActualPeriod][ActualPoint, :])) / 10.)
        etaD = etaD * ACO.Fact[ActualPeriod]
        etaD = etaD / np.sum(etaD)

        etaZ = 1. / ACO.ZenAn[ActualPeriod]
        #etaZ = 3*np.pi/(2*(ACO.ZenAn[ActualPeriod]+np.pi))
        etaZ = etaZ * ACO.Fact[ActualPeriod]
        etaZ = etaZ / sum(etaZ)

        #etaT = T_ant
        etaT = np.power(2,T_ant)-1.0
        #etaT = (np.power(2,T_ant)+128.)/192.
        etaT = 1.0 * etaT * ACO.Fact[ActualPeriod]
        etaT = etaT / sum(etaT)

        eta = np.power(etaD, 2) * np.power(etaZ, 2) * np.power(etaT, 1)
        #eta = etaD * etaZ * etaT

        if U <= Lambda:
            tau = ACO.Ph_ObsQ[ActualPoint, :] * ACO.Fact[ActualPeriod]
            tau = tau / np.sum(tau)
            index = ACO.ChooseStep(tau * eta)

        else:
            tau = ACO.Ph_TimeQ[ActualPoint, :] * ACO.Fact[ActualPeriod]
            tau = tau / np.sum(tau)
            index = ACO.ChooseStep(tau * eta)

        TotalT += ACO.Dist[ActualPeriod][ActualPoint, index] + ACO.ObservationTime
        LocalT += ACO.Dist[ActualPeriod][ActualPoint, index] + ACO.ObservationTime
        if(LocalT >= ACO.Times[ActualPeriod,1]):
            LocalT -= ACO.Times[ActualPeriod,1]
            NewPeriod = ActualPeriod+1
            T_ant += ACO.Times[ActualPeriod,1]
            if(NewPeriod%ACO.ND == 0):
                flagOverNight = True
            if(NewPeriod >= len(ACO.Times)):
                break
        else :
            NewPeriod = ActualPeriod

        if(not(flagOverNight)):
            # ObsQ += np.pi / 2 - ACO.ZenAn[NewPeriod][index]         #Observation quality based on zenith angle
            ObsQ += 1 / np.cos(ACO.ZenAn[NewPeriod][index])        #Observation quality based on air mass
            TimeQ += ACO.TimeFunc(T_ant[index])
            T_ant[index] = 0.0
            # ACO.Ph_ObsQ[ActualPoint, index] *= (1 - ACO.chi)
            # ACO.Ph_ObsQ[ActualPoint, index] += ACO.chi
            # ACO.Ph_TimeQ[ActualPoint, index] *= (1 - ACO.chi)
            # ACO.Ph_TimeQ[ActualPoint, index] += ACO.chi

            ActualPoint = index
            ActualPeriod = NewPeriod
            Path.append([ActualPoint, ActualPeriod])
        else:
            flagOverNight = False
            TotalT = sum(ACO.Times[0:NewPeriod,1])
            LocalT = 0.0
            ActualPeriod = NewPeriod

    # [ObsQ, TimeQ] = ACO.ObjectiveValues(np.array(Path, dtype=int))

    return [np.array(Path, dtype=int), ObsQ, TimeQ]

def Ant_Pipe(ACO,Lambda,conn):
    A = Ants_ext(ACO,Lambda)
    conn.send(A)
