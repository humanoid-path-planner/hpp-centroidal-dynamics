import unittest

from hpp_centroidal_dynamics import (
    LP_STATUS_OPTIMAL,
    Equilibrium,
    EquilibriumAlgorithm,
    SolverLP,
)
from numpy import array, asmatrix, cross


class TestCentroidalDynamics(unittest.TestCase):
    def test_all(self):
        # testing constructors
        eq = Equilibrium("test", 54.0, 4)
        eq = Equilibrium("test", 54.0, 4, SolverLP.SOLVER_LP_QPOASES)
        eq = Equilibrium("test", 54.0, 4, SolverLP.SOLVER_LP_QPOASES)
        eq = Equilibrium("test", 54.0, 4, SolverLP.SOLVER_LP_QPOASES, False)
        eq = Equilibrium("test", 54.0, 4, SolverLP.SOLVER_LP_QPOASES, False, 1)
        eq = Equilibrium("test", 54.0, 4, SolverLP.SOLVER_LP_QPOASES, True, 1, True)

        # whether useWarmStart is enable (True by default)
        previous = eq.useWarmStart()
        # enable warm start in solver (only for QPOases)
        eq.setUseWarmStart(False)
        self.assertNotEqual(previous, eq.useWarmStart())

        # access solver name
        self.assertEqual(eq.getName(), "test")

        z = array([0.0, 0.0, 1.0])
        P = asmatrix(
            array([array([x, y, 0]) for x in [-0.05, 0.05] for y in [-0.1, 0.1]])
        )
        N = asmatrix(array([z for _ in range(4)]))

        # setting contact positions and normals, as well as friction coefficients
        eq.setNewContacts(
            asmatrix(P), asmatrix(N), 0.3, EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_LP
        )

        c = asmatrix(array([0.0, 0.0, 1.0])).T

        # computing robustness of a given configuration,
        # first with no argument (0 acceleration, static equilibrium)
        status, robustness = eq.computeEquilibriumRobustness(c)
        self.assertEqual(status, LP_STATUS_OPTIMAL, "LP should not fail")
        self.assertGreater(robustness, 0, "first test should be in equilibrirum")

        # computing robustness of a given configuration with non zero acceleration
        ddc = asmatrix(array([1000.0, 0.0, 0.0])).T
        status, robustness = eq.computeEquilibriumRobustness(c, ddc)
        self.assertEqual(status, LP_STATUS_OPTIMAL, "LP should not fail")
        self.assertLess(robustness, 0, "first test should NOT be in equilibrirum")

        # now, use polytope projection algorithm
        eq.setNewContacts(
            asmatrix(P), asmatrix(N), 0.3, EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP
        )
        H, h = eq.getPolytopeInequalities()

        # c= asmatrix(array([0.,0.,1.])).T
        status, stable = eq.checkRobustEquilibrium(c)
        # TODO: This works well unless we add a --coverage in CXXFLAGS
        # self.assertEqual(status, LP_STATUS_OPTIMAL,
        # "checkRobustEquilibrium should not fail")
        self.assertTrue(stable, "lat test should be in equilibrirum")

        def compute_w(
            c, ddc, dL=array([0.0, 0.0, 0.0]), m=54.0, g_vec=array([0.0, 0.0, -9.81])
        ):
            w1 = m * (ddc - g_vec)
            return array(w1.tolist() + (cross(c, w1) + dL).tolist())

        def is_stable(
            H,
            c=array([0.0, 0.0, 1.0]),
            ddc=array([0.0, 0.0, 0.0]),
            dL=array([0.0, 0.0, 0.0]),
            m=54.0,
            g_vec=array([0.0, 0.0, -9.81]),
            robustness=0.0,
        ):
            w = compute_w(c, ddc, dL, m, g_vec)
            return (H.dot(-w) <= -robustness).all()

        self.assertTrue(is_stable(H))


if __name__ == "__main__":
    unittest.main()
