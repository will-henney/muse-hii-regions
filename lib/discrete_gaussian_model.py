import numpy as np
import scipy.stats
import astropy.modeling

def _E_cdf(x, x0, sig):
    "General case of any profile via the CDF"
    return scipy.stats.norm.cdf(x, loc=x0, scale=sig)

@astropy.modeling.custom_model
def DiscreteGaussianModel(x, amplitude=1.0, mean=0.0, stddev=1.0, bin_width=1.0):
    """
    A Gaussian profile, but integrated over bins of width bin_width (in units of x)

    """
    return  amplitude * (
        _E_cdf(x + bin_width/2, mean, stddev) - _E_cdf(x - bin_width/2, mean, stddev)
    )
