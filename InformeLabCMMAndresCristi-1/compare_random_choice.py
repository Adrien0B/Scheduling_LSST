import numpy as np
import matplotlib.pyplot as plt
from data_files import get_data_from_file

Colors = ['b','m','orange','black','r','pink','yellow','cyan','g']

R_x = []
R_y = []

a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_12-38-44_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_12-46-35_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_12-54-12_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_13-2-1_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_13-9-34_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_13-16-55_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_13-24-14_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_13-31-32_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_13-38-49_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_13-46-7_1000Iter_Pareto_Values")
R_x.append(a[0,:])
R_y.append(a[1,:])


NR_x = []
NR_y = []

a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_13-53-22_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_14-1-3_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_14-8-43_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_14-16-4_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_14-23-47_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_14-31-58_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_14-40-15_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_14-48-49_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_14-57-12_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-5-18_1000Iter_Pareto_Values")
NR_x.append(a[0,:])
NR_y.append(a[1,:])


R2_x = []
R2_y = []

a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-8-40_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-10-25_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-12-11_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-14-16_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-16-16_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-18-0_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-19-44_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-21-24_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-23-19_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-24-56_200Iter_Pareto_Values")
R2_x.append(a[0,:])
R2_y.append(a[1,:])


NR2_x = []
NR2_y = []

a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-26-33_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-28-13_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-29-52_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-31-31_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-33-8_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-34-49_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-36-28_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-38-11_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-39-53_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-41-34_200Iter_Pareto_Values")
NR2_x.append(a[0,:])
NR2_y.append(a[1,:])


R5_x = []
R5_y = []

a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-44-52_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-45-18_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-45-45_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-46-12_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-46-39_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-47-5_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-47-33_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-48-0_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-48-28_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/random/2017-6-15_15-48-57_50Iter_Pareto_Values")
R5_x.append(a[0,:])
R5_y.append(a[1,:])


NR5_x = []
NR5_y = []

a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-49-26_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-49-55_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-50-21_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-50-49_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-51-15_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-51-42_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-52-8_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-52-38_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-53-7_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])
a = get_data_from_file("videos/compare_choice_random/non_random/2017-6-15_15-53-35_50Iter_Pareto_Values")
NR5_x.append(a[0,:])
NR5_y.append(a[1,:])



fig=plt.figure()
for i in range(np.size(R_x,0)):
    plt.plot(R_x[i],R_y[i],marker='+',label=str(i+1),c=Colors[0])
for i in range(np.size(NR_x,0)):
    plt.plot(NR_x[i],NR_y[i],marker='+',label=str(i+1),c=Colors[1])
for i in range(np.size(R2_x,0)):
    plt.plot(R2_x[i],R2_y[i],marker='+',label=str(i+1),c=Colors[2])
for i in range(np.size(NR2_x,0)):
    plt.plot(NR2_x[i],NR2_y[i],marker='+',label=str(i+1),c=Colors[3])
for i in range(np.size(R5_x,0)):
    plt.plot(R5_x[i],R5_y[i],marker='+',label=str(i+1),c=Colors[4])
for i in range(np.size(NR5_x,0)):
    plt.plot(NR5_x[i],NR5_y[i],marker='+',label=str(i+1),c=Colors[5])
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

[x1,y1,xerr1,yerr1] = medianCurve(R_x,R_y)
[x2,y2,xerr2,yerr2] = medianCurve(NR_x,NR_y)
[x3,y3,xerr3,yerr3] = medianCurve(R2_x,R2_y)
[x4,y4,xerr4,yerr4] = medianCurve(NR2_x,NR2_y)
[x5,y5,xerr5,yerr5] = medianCurve(R5_x,R5_y)
[x6,y6,xerr6,yerr6] = medianCurve(NR5_x,NR5_y)

plt.figure()
plt.errorbar(x1,y1,xerr=xerr1,yerr=yerr1,label='Random choice (1000Iters)',c=Colors[0])
plt.errorbar(x2,y2,xerr=xerr2,yerr=yerr2,label='Deterministic choice (1000Iters)',c=Colors[1])
plt.errorbar(x3,y3,xerr=xerr3,yerr=yerr3,label='Random choice (200Iters)',c=Colors[2])
plt.errorbar(x4,y4,xerr=xerr4,yerr=yerr4,label='Deterministic choice (200Iters)',c=Colors[3])
plt.errorbar(x5,y5,xerr=xerr5,yerr=yerr5,label='Random choice (50Iters)',c=Colors[4])
plt.errorbar(x6,y6,xerr=xerr6,yerr=yerr6,label='Deterministic choice (50Iters)',c=Colors[5])
plt.legend()
plt.show()


