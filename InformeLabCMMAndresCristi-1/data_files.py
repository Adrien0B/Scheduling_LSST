import numpy as np
import matplotlib.pyplot as plt


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
