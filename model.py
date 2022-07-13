from config import OCnfig
from pyz3_utils import MySolver
from z3 import Implies, Or, Sum


class Flow:
    def __init__(self, name: str, c: Config, s: MySolver):
        self.c = c
        self.s = s

        self.A = s.Real(f"{name}_A")
        self.S = s.Real(f"{name}_S")
        self.L = s.Real(f"{name}_L")
        self.Ld = s.Real(f"{name}_Ld")


class Timestep:
    def __init__(self, name: str, c: Config, s: MySolver):
        self.c = c
        self.s = s

        self.time = s.Real(f"{name}_time")
        self.flows = [Flow(f"{name}_flow{f}", s, c)
                      for f in range(c.F)]
        self.W = s.Real(f"{name}_W")

        # Totals of all the quantities across flows
        self.A = s.Real(f"{name}_A")
        self.S = s.Real(f"{name}_S")
        self.L = s.Real(f"{name}_L")

        self.total()

    def total(self):
        ''' Values of individual flows sum to the total '''
        c = self.c
        s = self.s

        s.add(self.A == Sum(self.flows[f].A for f in range(c.F)))
        s.add(self.S == Sum(self.flows[f].S for f in range(c.F)))
        s.add(self.L == Sum(self.flows[f].L for f in range(c.F)))


class Variables:
    def __init__(self, c: Config, s: MySolver):
        self.times = [Timestep(f"time{t}", c, s) for t in range(c.T)]

        if c.compose = False:
            # Upper bound on amount of data needed in queue for loss to occur
            self.epsilon = s.Real("epsilon")
            s.add(self.epsilon >= 0)

        if not c.inf_buf:
            self.buf = s.Real("buf")
            if c.buf_size is not None:
                s.add(self.buf == c.buf_size)
            s.add(self.buf > 0)

        self.monotone()
        self.initial()

    def monotone(self):
        c = self.c
        s = self.s
        times = self.times

        for t in range(1, c.T):
            pre = times[t-1]
            nex = times[t]

            s.add(pre.W <= times[t].W)

            s.add(pre.time < nex.time)
            s.add(pre.A <= nex.A)
            s.add(pre.S <= nex.S)
            s.add(pre.L <= nex.L)

            for f in range(c.F):
                pre = times[t-1].flows[f]
                nex = times[t].flows[f]

                s.add(pre.time < nex.time)
                s.add(pre.A <= nex.A)
                s.add(pre.S <= nex.S)
                s.add(pre.L <= nex.L)
                s.add(pre.Ld <= nex.Ld)

    def initial(self):
        ''' Set initial conditions '''
        c = self.c
        s = self.s
        times = self.times

        s.add(times[0].time == 0)

        # Since values are invariant to y-shift, we don't need to enfore
        # initial conditions. However we can always enforce one of them,
        # for pretty plots.
        s.add(times[0].S == 0)


        for f in range(c.F):
            init = times[0].flows[f]
            # What does negative loss even mean?
            s.add(init.L >= 0)
            s.add(init.Ld >= 0)

    def network(self):
        # Invariants that apply independently on each timestep
        c = self.c
        s = self.s
        times = self.times

        for t in range(c.T):
            ts = self.times[t]
            for f in range(c.F):
                fl = ts.flow[f]
                s.add(fl.S <= fl.A - fl.L)
            s.add(ts.S <= c.C * ts.time - ts.W)

            # Do things at time ts.time - c.D

            # To begin with, if ts.time - c.D > 0, it should exist in the past
            s.add(Or(ts.time < c.D,
                *[times[pt].time == ts.time for pt in range(t)]))

            # If ts.time - c.D < 0, then give maximum slack. This corresponds
            # to no wastage when t < 0
            s.add(c.C * (ts.time - c.D) - times[0].W <= ts.W)

            for pt in range(t):
                pts = self.times[t]
                s.add(Implies(ts.time - c.D == pts.time,
                              c.C * pts.time - pts.W <= ts.S))

            if c.compose:
                if t > 0:
                    s.add(Implies(ts.W > times[t-1].W,
                                  ts.A - ts.L <= c.C * ts.time - ts.W))
            else:
                if t > 0:
                    s.add(Implies(ts.W > times[t-1].W,
                                  ts.A - ts.L <= ts.S + self.epsilon))

            if not c.inf_buf:
                if t == 0:
                    continue
