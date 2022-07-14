from config import Config
from pyz3_utils import MySolver, Variables
from z3 import And, Implies, Or, Sum


class Flow(Variables):
    def __init__(self, name: str, c: Config, s: MySolver):
        self.c = c
        self.s = s

        self.A = s.Real(f"{name}_A")
        self.S = s.Real(f"{name}_S")
        self.L = s.Real(f"{name}_L")
        self.Ld = s.Real(f"{name}_Ld")
        self.rtt = s.Real(f"{name}_rtt")


class Timestep(Variables):
    def __init__(self, name: str, c: Config, s: MySolver):
        self.c = c
        self.s = s

        self.time = s.Real(f"{name}_time")
        self.flows = [Flow(f"{name}_flow{f}", c, s)
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

        s.add(self.A == Sum([self.flows[f].A for f in range(c.F)]))
        s.add(self.S == Sum([self.flows[f].S for f in range(c.F)]))
        s.add(self.L == Sum([self.flows[f].L for f in range(c.F)]))


class ModelVariables(Variables):
    def __init__(self, c: Config, s: MySolver):
        self.c = c
        self.s = s

        self.times = [Timestep(f"t{t}", c, s) for t in range(c.T)]

        if c.compose == False:
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
        self.network()
        self.measure_losses_and_rtt()

    def monotone(self):
        c = self.c
        s = self.s
        times = self.times

        for t in range(1, c.T):
            pre = times[t-1]
            nex = times[t]

            s.add(pre.W <= nex.W)
            s.add(c.C * pre.time - pre.W <= c.C * nex.time - nex.W)

            s.add(pre.time < nex.time)
            s.add(pre.A <= nex.A)
            s.add(pre.S <= nex.S)
            s.add(pre.L <= nex.L)

            for f in range(c.F):
                pre = times[t-1].flows[f]
                nex = times[t].flows[f]

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
        ''' The heart of the CCAC model '''
        c = self.c
        s = self.s
        times = self.times

        for t in range(c.T):
            ts = self.times[t]
            for f in range(c.F):
                fl = ts.flows[f]
                s.add(fl.S <= fl.A - fl.L)
            s.add(ts.S <= c.C * ts.time - ts.W)

            # Do things at time ts.time - c.D

            # To begin with, if ts.time - c.D > 0, it should exist in the past
            s.add(Or(ts.time < c.D,
                *[times[pt].time == ts.time - c.D for pt in range(t)]))

            # If ts.time - c.D < 0, then give maximum slack. This corresponds
            # to no wastage when t < 0
            s.add(c.C * (ts.time - c.D) - times[0].W <= ts.S)

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
                s.add(ts.A - ts.L <= c.C * ts.time - ts.W + self.buf)

                if t == 0:
                    continue
                # We can make loss deterministic since we assume all curves are
                # joined by straight lines. This does not lose generality since
                # `times[t].time` is a variable. Thus Z3 can approximate any
                # curve it likes with a piecewise linear curve (well...as long
                # as it fits within c.T points)

                s.add(Implies(ts.L > times[t-1].L,
                              ts.A - ts.L ==
                              c.C * ts.time - ts.W + self.buf))
            else:
                s.add(ts.L == times[0].L)

    def measure_losses_and_rtt(self):
        '''Constrain Ld based on A, S and L. If these constraints are added, it will
        automatically calculate delay since for every t, it will force t' such
        that A[t'] - L[t'] = S[t] if S[t] >= A[0]

        '''
        c = self.c
        s = self.s
        times = self.times

        for f in range(c.F):
            s.add(times[0].flows[f].rtt > 0)

        for t in range(0, c.T):
            for f in range(c.F):
                cur = times[t].flows[f]
                zero = times[0].flows[f]

                # For every S, ensure a corresponding A exists
                if t > 0:
                    s.add(Or(cur.S < zero.A - zero.L,
                             *[times[tp].flows[f].A - times[tp].flows[f].L
                               == cur.S
                               for tp in range(t)]))
                # Duh
                s.add(cur.Ld <= cur.L)

                # If S is too small, delay can stretch to -inf
                s.add(Implies(cur.S < zero.A - zero.L,
                              cur.rtt >= c.R + times[t].time - times[0].time))

                # Set Ld and delay
                for tp in range(0, t):
                    pre = times[tp].flows[f]
                    # if A-L does not change, we will pick the time corresponding
                    # to the first one
                    if tp > 0:
                        not_subseq = (pre.A - pre.L != times[tp-1].flows[f].A
                                      - times[tp-1].flows[f].L)
                    else:
                        not_subseq = True

                    s.add(
                        Implies(
                            And(
                                pre.A - pre.L == cur.S,
                                not_subseq
                            ),
                            And(
                                cur.Ld == pre.L,
                                cur.rtt == c.R + times[t].time - times[tp].time
                            )))

if __name__ == "__main__":
    from plot import plot
    from pyz3_utils import run_query

    c = Config()
    c.unsat_core = True
    c.T = 10
    c.inf_buf = False
    c.check()
    s = MySolver()
    v = ModelVariables(c, s)

    s.add(v.times[-1].time >= 5)
    s.add(v.times[3].flows[0].rtt == 1.5)
    for t in range(1, c.T):
        s.add(v.times[t].A >= v.times[t-1].A + c.C / 2)

    res = run_query(c, s, v)
    print(res.satisfiable)
    if res.satisfiable == "sat":
        plot(c, res.v)
