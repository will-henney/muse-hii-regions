import numpy as np
import pyneb as pn

hi = pn.RecAtom("H", 1)
tem, den = 12500, 100
R0 = (hi.getEmissivity(tem, den, wave=6563)
      / hi.getEmissivity(tem, den, wave=4861))
print(f"Intrinsic Balmer decrement: R0 = {R0:.3f}")

# Set up reddening law
rc = pn.RedCorr()
rc.R_V = 2.74
rc.FitzParams = [-4.96, 2.26, 0.39, 0.6, 4.6, 1.0]
rc.law = "F99"

decrements = np.array(
    [3.44, 3.23, 3.63, 3.62, 3.01, 3.43]
)

rc.setCorr(obs_over_theo=decrements / R0, wave1=6563.0, wave2=4861.0)
print(f"E(B-V) = {np.round(rc.E_BV, 2)}")
print(f"c(H beta) = {np.round(rc.cHbeta, 2)}")
for wav in 4000, 5000, 6000, 7000, 8000, 9000:
    factor = 10**(0.4 * rc.E_BV * (rc.X(wav) - rc.X(4861)))
    factor2 = rc.getCorr(wav) / rc.getCorr(4861)
    print(wav, np.round(factor, 2))
    print(wav, np.round(factor2, 2))
