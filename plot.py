from config import Config
import matplotlib.pyplot as plt
import numpy as np

def plot(c: Config, v):
    times = np.asarray(list(map(lambda t: t.time, v.times)))
    W = np.asarray(list(map(lambda t: t.W, v.times)))
    A = np.asarray(list(map(lambda t: t.A, v.times)))
    S = np.asarray(list(map(lambda t: t.S, v.times)))
    L = np.asarray(list(map(lambda t: t.L, v.times)))
    ubound = c.C * times - W
    print(S)
    print(times)

    fig, ax1 = plt.subplots(1, 1)

    ax1.plot(times, ubound, marker='o', c='black', label='')
    ax1.plot(times + c.D, ubound, marker='o', c='black', label='')
    ax1.plot(times, A, marker='o', c='blue', label='A')
    ax1.plot(times, A-L, marker='o', c='lightblue', label='A-L')
    ax1.plot(times, S, marker='o', c='red', label='S')
    ax1.plot(times, L, marker='o', c='yellow', label='L')

    plt.show()
