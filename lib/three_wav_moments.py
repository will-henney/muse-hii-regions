import numpy as np
from scipy.special import erf
import scipy.stats as ss
def _E_erf(x, x0, sig):
    "Special case of erf for gaussian profile"
    return 0.5 * erf((x - x0) / (sig * np.sqrt(2)))

_PROFILE = ss.norm       # Gaussian
#_PROFILE = ss.cauchy            # Lorentzian
def _E_cdf(x, x0, sig):
    "General case of any profile via the CDF"
    return _PROFILE.cdf(x, loc=x0, scale=sig)

# Use the general CDF form so that functional form of profile can be
# changed (see the _PROFILE variable above)
E = _E_cdf

def M0(x0, sig):
    return E(1.5, x0, sig) - E(-1.5, x0, sig)

def M1(x0, sig):
    return E(1.5, x0, sig) + E(-1.5, x0, sig) - E(0.5, x0, sig) - E(-0.5, x0, sig)

def M2(x0, sig):
    return E(1.5, x0, sig) + E(-0.5, x0, sig) - E(0.5, x0, sig) - E(-1.5, x0, sig)
