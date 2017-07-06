import os
import numpy as np
import matplotlib.pyplot as plt

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


def get_data_from_file(title):
    f = open(title,'r')
    x = f.readline()
    y = f.readline()
    f.close()
    x = x.replace('[','')
    x = x.replace(']\n','')
    y = y.replace('[','')
    y = y.replace(']\n','')
    x = x.split(',')
    y = y.split(',')
    x = np.array(map(float,x))
    y = np.array(map(float,y))
    return np.array([x,y])

def create_list(directory):
    if(directory[-1]!='/'):
        directory = directory + '/'
    F = os.listdir(directory)
    List_x = []
    List_y = []
    for f in F:
        if(f[-14:]=="_Pareto_Values"):
            a = get_data_from_file(directory+f)
            List_x.append(a[0,:])
            List_y.append(a[1,:])
    return [List_x,List_y]


def explore_directory(directory):
    Titles = []
    Data = []
    Median = []
    D = os.walk(directory)
    (_,dirnames,_) = D.next()
    for dir in dirnames:
        (dirpath, _, _) = D.next()
        L = create_list(dirpath)
        Titles.append(dir)
        Data.append(L)
        L_Err = medianCurve(L[0],L[1])
        Median.append(L_Err)
    return [Titles,Data,Median]

def plot_directory(directory):
    [Titles,Data,Median] = explore_directory(directory)
    Colors = ['b', 'm', 'orange', 'black', 'r', 'pink', 'yellow', 'cyan', 'g', 'salmon', 'slategrey', 'brown']
    fig = plt.figure()
    for k in range(len(Titles)):
        for i in range(np.size(Data[k][0], 0)):
            plt.plot(Data[k][0][i], Data[k][1][i], marker='+', label=Titles[k], c=Colors[k])
    plt.show()

    plt.figure()
    for k in range(len(Titles)):
        plt.errorbar(Median[k][0],Median[k][1],xerr=Median[k][2],yerr=Median[k][3],label=Titles[k],c=Colors[k])
    try:
        plt.legend()
    except:
        pass
    plt.show()
