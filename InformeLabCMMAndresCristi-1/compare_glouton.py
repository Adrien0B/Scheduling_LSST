import numpy as np
import matplotlib.pyplot as plt


Colors = ['b','m','orange','black','r','pink','yellow','cyan','g','salmon','slategrey','brown']

Greedy_x = np.array([99.9,99.9,102.3,102.8,103.0,103.4,104.0,104.2])

Greedy_y = np.array([3090.0,3078.0,3026.0,3014.0,2954.0,2912.0,2908.0,2738.0])


Iter_x = np.array([[99.1,99.9,100.7,100.8,101.9,102.1,102.3,103.8,105.6,106.9,107.9,108.3,109.0,109.8,109.9,110.1],
                   [96.6, 100.1, 104.2, 107.8, 107.8, 109.9, 110.4, 110.7],
                   [100.7, 102.7, 103.8, 107.5, 107.6, 107.8, 108.3, 110.0, 110.1, 110.3, 110.8, 111.8],
                   [102.8, 103.3, 105.0, 105.9, 106.1, 107.3, 107.3, 107.5, 107.8, 108.9, 109.9, 110.5, 110.6],
                   [101.2, 102.9, 104.0, 104.9, 106.8, 107.6, 107.7, 107.8, 108.0, 109.7, 110.5],
                   [97.9, 100.3, 101.4, 102.7, 102.7, 103.0, 103.1, 103.9, 104.9, 106.9, 108.0, 108.2, 109.0, 109.3,
                    109.4, 109.9, 110.8],
                   [99.1, 99.9, 102.8, 103.0, 103.8, 106.4, 106.5, 107.1, 107.3, 108.7, 108.8, 109.1, 110.4, 110.5],
                   [95.8, 101.6, 101.8, 102.5, 104.0, 104.6, 105.3, 106.8, 108.7, 110.3, 111.0],
                   [102.5, 104.4, 105.8, 107.6, 108.9, 109.1, 110.8]])

Iter_y = np.array([[3580.0,3522.0,3520.0,3482.0,3422.0,3404.0,3400.0,3396.0,3370.0,3330.0,3316.0,3214.0,3122.0,3038.0,2994.0,2768.0],
                   [3544.0, 3512.0, 3494.0, 3346.0, 3132.0, 3100.0, 2978.0, 2942.0],
                   [3534.0, 3484.0, 3372.0, 3330.0, 3276.0, 3182.0, 3150.0, 3128.0, 2926.0, 2876.0, 2842.0, 2824.0],
                   [3526.0, 3482.0, 3390.0, 3340.0, 3314.0, 3286.0, 3186.0, 3148.0, 3136.0, 3118.0, 3062.0, 2898.0,
                    2634.0],
                   [3548.0, 3512.0, 3440.0, 3358.0, 3300.0, 3178.0, 3140.0, 3136.0, 3128.0, 3106.0, 3092.0],
                   [3512.0, 3506.0, 3498.0, 3468.0, 3462.0, 3454.0, 3404.0, 3392.0, 3388.0, 3342.0, 3314.0, 3230.0,
                    3122.0, 3084.0, 3076.0, 2876.0, 2822.0],
                   [3594.0, 3530.0, 3494.0, 3434.0, 3408.0, 3404.0, 3280.0, 3270.0, 3226.0, 3164.0, 3128.0, 3074.0,
                    3036.0, 2842.0],
                   [3534.0, 3518.0, 3504.0, 3436.0, 3400.0, 3354.0, 3316.0, 3310.0, 3306.0, 3120.0, 2848.0],
                   [3452.0, 3432.0, 3368.0, 3346.0, 3196.0, 3194.0, 3092.0]])


fig=plt.figure()
plt.plot(Greedy_x,Greedy_y,marker='+',label='Greedy',c=Colors[0])
for i in range(np.size(Iter_x,0)):
    plt.plot(Iter_x[i], Iter_y[i], marker='+', label='ACO', c=Colors[1])
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

[x1,y1,xerr1,yerr1] = medianCurve(Iter_x,Iter_y)

plt.figure()
plt.errorbar(Greedy_x,Greedy_y,label='Greedy',c=Colors[0])
plt.errorbar(x1,y1,xerr=xerr1,yerr=yerr1,label='ACO',c=Colors[1])
plt.legend()
plt.show()



