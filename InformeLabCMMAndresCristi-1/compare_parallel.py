import numpy as np
import matplotlib.pyplot as plt
from data_files import get_data_from_file

Colors = ['b','m','orange','black','r','pink','yellow','cyan','g']

P_x = []
P_y = []

a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-44-45_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-45-12_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-45-12_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-46-10_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-46-41_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-47-10_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-47-39_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-48-12_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-48-46_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/parallel/2017-6-19_12-49-14_100Iter_Pareto_Values")
P_x.append(a[0,:])
P_y.append(a[1,:])



S_x = []
S_y = []

a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-26-57_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-27-49_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-28-43_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-29-37_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-30-25_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-31-15_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-32-9_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-33-3_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-33-55_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/sequential/2017-6-16_9-34-46_100Iter_Pareto_Values")
S_x.append(a[0,:])
S_y.append(a[1,:])



Pp_x = []
Pp_y = []

a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_14-52-19_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_14-54-24_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_14-56-29_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_14-58-30_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_15-0-31_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_15-2-20_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_15-4-13_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_15-5-59_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_15-7-46_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pool/2017-6-19_15-9-32_100Iter_Pareto_Values")
Pp_x.append(a[0,:])
Pp_y.append(a[1,:])


PPipe_x = []
PPipe_y = []

a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-7-59_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-8-26_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-8-54_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-9-21_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-9-47_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-10-15_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-10-42_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-11-7_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-11-34_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])
a = get_data_from_file("videos/compare_parallel/para_pipe/2017-6-20_17-12-3_100Iter_Pareto_Values")
PPipe_x.append(a[0,:])
PPipe_y.append(a[1,:])


fig=plt.figure()
for i in range(np.size(P_x,0)):
    plt.plot(P_x[i],P_y[i],marker='+',label=str(i+1),c=Colors[0])
for i in range(np.size(S_x,0)):
    plt.plot(S_x[i],S_y[i],marker='+',label=str(i+1),c=Colors[1])
for i in range(np.size(Pp_x,0)):
    plt.plot(Pp_x[i],Pp_y[i],marker='+',label=str(i+1),c=Colors[2])
for i in range(np.size(PPipe_x,0)):
    plt.plot(PPipe_x[i],PPipe_y[i],marker='+',label=str(i+1),c=Colors[3])
plt.show()

def oneD(A):
    a = []
    for i in A:
        for j in i:
            a.append(j)
    return np.array(a)

def medianCurve(Eta_x,Eta_y):
    A = oneD(Eta_x)
    B = oneD(Eta_y)
    a = np.array(A,dtype='int')
    U = np.unique(a)

    X = []
    Y = []
    X_err = []
    Y_err = []

    for u in U:
        I = a==u
        N = np.sum(I)
        Ai = A[I]
        Bi = B[I]
        X.append(np.mean(Ai))
        Y.append(np.mean(Bi))
        X_err.append(np.std(Ai)/np.sqrt(N))
        Y_err.append(np.std(Bi)/np.sqrt(N))

    return [X,Y,X_err,Y_err]

[x1,y1,xerr1,yerr1] = medianCurve(P_x,P_y)
[x2,y2,xerr2,yerr2] = medianCurve(S_x,S_y)
[x3,y3,xerr3,yerr3] = medianCurve(Pp_x,Pp_y)
[x4,y4,xerr4,yerr4] = medianCurve(PPipe_x,PPipe_y)

plt.figure()
plt.errorbar(x1,y1,xerr=xerr1,yerr=yerr1,label='Parallel (100Iters)',c=Colors[0])
plt.errorbar(x2,y2,xerr=xerr2,yerr=yerr2,label='Sequential (100Iters)',c=Colors[1])
plt.errorbar(x3,y3,xerr=xerr3,yerr=yerr3,label='Parallel Pools(4) (100Iters)',c=Colors[2])
plt.errorbar(x4,y4,xerr=xerr4,yerr=yerr4,label='Parallel Pipes (100Iters)',c=Colors[3])
plt.legend()
plt.show()


