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

    fig, (ax1, ax2_rtt) = plt.subplots(2, 1)
    ax2_loss = ax2_rtt.twinx()

    ax1.set_xlabel("Time")
    ax2_rtt.set_xlabel("Time")
    ax2_rtt.set_ylabel("RTT")
    ax2_loss.set_ylabel("Loss")

    ax1.plot(times, ubound, marker='o', c='black', label='')
    ax1.plot(times + c.D, ubound, marker='o', c='black', label='')
    if not c.inf_buf:
        ax1.plot(times, ubound + v.buf, marker='o', c='green', label='')
    ax1.plot(times, A, marker='o', c='blue', label='A')
    ax1.plot(times, A-L, marker='o', c='lightblue', label='A-L')
    ax1.plot(times, S, marker='o', c='red', label='S')

    for f in range(c.F):
        # For some reason, z3 doesn't always populate RTT even if it is in the
        # constraints (must be a bug in z3). This is a temporary hack to get
        # around this
        for t in range(c.T):
            if 'rtt' not in v.times[t].flows[f].__dict__:
                v.times[t].flows[f].__dict__['rtt'] = -1

        rtts = np.asarray(list(map(lambda t: t.flows[f].rtt, v.times)))
        ax2_rtt.plot(times, rtts, marker='o', label=f'rtt{f}')

        f_L = np.asarray(list(map(lambda t: t.flows[f].L, v.times)))
        f_Ld = np.asarray(list(map(lambda t: t.flows[f].Ld, v.times)))

        ax2_loss.plot(times, f_L, marker='o', c='olive', label=f'L{f}')
        ax2_loss.plot(times, f_Ld, marker='o', c='yellow', label=f'Ld{f}')

    plt.legend()
    plt.show()
