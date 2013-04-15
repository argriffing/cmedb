from numpy import testing
import numpy as np
import scipy.linalg

import mmpp


class TestMe(testing.TestCase):

    def test_expm(self):

        # define some values
        # @param a: rate from off to on
        # @param w: rate from on to off
        # @param r: poisson event rate
        # @param t: elapsed time
        a = 1.2
        w = 2.3
        r = 0.9
        t = 0.5

        # get the P block using mmpp
        P_block_cython = mmpp.get_mmpp_block(a, w, r, t)

        # get the P block using expm
        pre_Q = np.array([
            [0, a, 0],
            [w, 0, r],
            [0, 0, 0],
            ], dtype=float)
        Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
        P = scipy.linalg.expm(Q*t)
        P_block_expm = P[:2, :2]

        # assert that the two blocks are close
        testing.assert_allclose(P_block_cython, P_block_expm)


if __name__ == '__main__':
    testing.run_module_suite()

