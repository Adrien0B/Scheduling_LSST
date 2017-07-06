import numpy as np
import matplotlib.pyplot as plt

E0 = np.load("Exec_Times_0Iters.npy")
E20 = np.load("Exec_Times_20Iters.npy")
E50 = np.load("Exec_Times_50Iters.npy")
E100 = np.load("Exec_Times_100Iters.npy")
E200 = np.load("Exec_Times_200Iters.npy")
E500 = np.load("Exec_Times_500Iters.npy")
E1000 = np.load("Exec_Times_1000Iters.npy")

E1000 = np.append(E1000,[[92],[16437.633327007294]],axis=1)

A20 = E20 - E0
A50 = E50 - E0
A100 = E100 - E0
A200 = E200 - E0
A500 = E500 - E0
A1000 = E1000 - E0

N20 = A20/20.0
N20[0,:] = E20[0,:]
N50 = A50/50.0
N50[0,:] = E50[0,:]
N100 = A100/100.0
N100[0,:] = E100[0,:]
N200 = A200/200.0
N200[0,:] = E200[0,:]
N500 = A500/500.0
N500[0,:] = E500[0,:]
N1000 = A1000/1000.0
N1000[0,:] = E1000[0,:]


if(0):
    # Raw Data : Number of Nights vs Total Execution Time
    plt.plot(E0[0,:],E0[1,:],marker='+',label='0 Iterations')
    plt.plot(E20[0,:],E20[1,:],marker='+',label='20 Iterations')
    plt.plot(E50[0,:],E50[1,:],marker='+',label='50 Iterations')
    plt.plot(E100[0,:],E100[1,:],marker='+',label='100 Iterations')
    plt.plot(E200[0,:],E200[1,:],marker='+',label='200 Iterations')
    plt.plot(E500[0,:],E500[1,:],marker='+',label='500 Iterations')
    plt.plot(E1000[0,:],E1000[1,:],marker='+',label='1000 Iterations')
    plt.legend()
    plt.show()

if(0):
    # Number of Nights vs Execution Time for Iterations
    plt.plot(E20[0,:],A20[1,:],marker='+',label='20 Iterations')
    plt.plot(E50[0,:],A50[1,:],marker='+',label='50 Iterations')
    plt.plot(E100[0,:],A100[1,:],marker='+',label='100 Iterations')
    plt.plot(E200[0,:],A200[1,:],marker='+',label='200 Iterations')
    plt.plot(E500[0,:],A500[1,:],marker='+',label='500 Iterations')
    plt.plot(E1000[0,:],A1000[1,:],marker='+',label='1000 Iterations')
    plt.legend()
    plt.show()

if(0):
    # Number of nights vs Execution Time for One Iteration
    plt.plot(N20[0,:],N20[1,:],marker='+',label='20 Iterations')
    plt.plot(N50[0,:],N50[1,:],marker='+',label='50 Iterations')
    plt.plot(N100[0,:],N100[1,:],marker='+',label='100 Iterations')
    plt.plot(N200[0,:],N200[1,:],marker='+',label='200 Iterations')
    plt.plot(N500[0,:],N500[1,:],marker='+',label='500 Iterations')
    plt.plot(N1000[0,:],N1000[1,:],marker='+',label='1000 Iterations')
    plt.legend()
    plt.show()


from data_files import medianCurve


X = np.array([N20[0,:],N50[0,:],N100[0,:],N200[0,:],N500[0,:],N1000[0,:]])
Y = np.array([N20[1,:],N50[1,:],N100[1,:],N200[1,:],N500[1,:],N1000[1,:]])
# X = np.array([N20[0,:-1],N50[0,:-1],N100[0,:-1],N200[0,:-1],N500[0,:-1],N1000[0,:-1]])
# Y = np.array([N20[1,:-1],N50[1,:-1],N100[1,:-1],N200[1,:-1],N500[1,:-1],N1000[1,:-1]])
[x,y,xerr,yerr] = medianCurve(X,Y)

if(0):
    plt.errorbar(x,y,xerr=xerr,yerr=yerr)
    plt.show()

P1 = np.polyfit(x,y,1)
P2 = np.polyfit(x,y,2)
X_ = np.vstack((np.power(x,2), x, np.ones_like(x))).T
X_ = X_[:,:-1]
p1 = np.linalg.lstsq(X_[:,1:], y)[0]
p1 = np.append(p1,0.0)
p2 = np.linalg.lstsq(X_, y)[0]
p2 = np.append(p2,0.0)


if(1):
    # Execution Times and Fit to data (1st and 2nd degree polynomial)
    plt.scatter(N20[0, :], N20[1, :], marker='+', label='20 Iterations')
    plt.scatter(N50[0, :], N50[1, :], marker='+', label='50 Iterations')
    plt.scatter(N100[0, :], N100[1, :], marker='+', label='100 Iterations')
    plt.scatter(N200[0, :], N200[1, :], marker='+', label='200 Iterations')
    plt.scatter(N500[0, :], N500[1, :], marker='+', label='500 Iterations')
    plt.scatter(N1000[0, :], N1000[1, :], marker='+', label='1000 Iterations')
    plt.scatter(x,y,label='Average values')
    x = np.linspace(0,100,100)
    plt.plot(x,np.polyval(P1,x),label="Fit to data (1st degree)")
    plt.plot(x,np.polyval(p1,x),label="Fit to data (1st degree through origin)")
    plt.plot(x,np.polyval(P2,x),label="Fit to data (2nd degree)")
    plt.plot(x,np.polyval(p2,x),label="Fit to data (2nd degree through origin)")
    plt.legend()
    plt.show()



if(1):
    plt.plot([0,20,50,100,200,500,1000],np.array([E0[1,:],E20[1,:],E50[1,:],E100[1,:],E200[1,:],E500[1,:],E1000[1,:]]),marker='+')
    plt.show()

