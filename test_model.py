from config import Config
from model import ModelVariables
from pyz3_utils import MySolver, Piecewise, run_query
import unittest
from z3 import Or

class TestModel(unittest.TestCase):
    def test_max_delta_t(self):
        ''' Gap between time cannot be lrager than min(c.R, c.D) '''
        for (R, D) in [(1, 1), (1, 2), (2, 0.5)]:
            c = Config()
            c.R = R
            c.D = D
            c.check()
            s = MySolver()
            v = ModelVariables(c, s)

            s.add(Or(*[v.times[t].time - v.times[t-1].time > min(c.D, c.R)
                       for t in range(1, c.T)]))
            res = run_query(c, s, v)
            self.assertEqual(res.satisfiable, "unsat")

    def test_delta_t(self):
        ''' Run `verify` on `delta_t`s which are `Piecewise` objects '''
        c = Config()
        c.check()
        s = MySolver()
        v = ModelVariables(c, s)
        for delta_t in v.delta_t[1:]:
            delta_t.verify(s)

    def test_A_L_monotone(self):
        ''' Test that A-L is monotonic '''
        c = Config()
        c.check()
        s = MySolver()
        v = ModelVariables(c, s)

        s.add(Or(*[v.times[t].A - v.times[t].L <
                   v.times[t-1].A - v.times[t-1].L
                   for t in range(1, c.T)]))
        res = run_query(c, s, v)
        self.assertEqual(res.satisfiable, "unsat")


if __name__ == '__main__':
    # Run using 'python3 -m unittest test_model.py'
    unittest.main()
