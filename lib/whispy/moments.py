"""Functions to work with velocity moments of emission lines from MUSE cubes

The lines need to have been continuum-subtracted first - see extract.py

Author: Will Henney, IRyA-UNAM, 2021
"""
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from mpdaf.obj import Cube
import astropy.units as u
import astropy.constants as const
import pandas as pd

LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value
FIGPATH = Path(".")
SAVEPATH = Path(".")


def find_moments(cube):
    """
    Returns the normalized wavelength moments: mom0, mom1, mom2

    mom0 is sum over wavelength
    mom1 is mean wavelength
    mom2 is rms wavelength width
    """
    # TODO: calculate variance arrays
    wavcube = cube.clone(np.ones, np.zeros)
    wavcube.data *= cube.wave.coord()[:, None, None]
    wavcube.unit = u.angstrom
    # zeroth moment: sum
    mom0 = cube.sum(axis=0)
    # first moment: mean
    mom1 = mom0.copy()
    mom1.data = np.sum(cube.data * wavcube.data, axis=0) / mom0
    mom1.unit = u.angstrom
    # second moment: sigma
    mom2 = mom0.copy()
    mom2.data = np.sum(cube.data * (wavcube.data - mom1.data) ** 2, axis=0) / mom0
    mom2.data = np.sqrt(mom2.data)
    mom2.unit = u.angstrom
    return mom0, mom1, mom2


def save_moments_to_fits(
    moments,
    *,
    rebin=1,
    flabel="ion",
    label="0000",
    restwav,
    irange=None,
    vrange=None,
    srange=None,
):
    """Write FITS files of velocity moment maps"""

    # Optionally do rebinning at `rebin` x `rebin`
    mom0 = moments[0].rebin(rebin)
    mom1 = moments[1].rebin(rebin)
    mom2 = moments[2].rebin(rebin)
    # Convert from wavelength to velocity units for 1st and 2nd moment
    vmean = LIGHT_SPEED_KMS * (mom1 - restwav) / restwav
    sigma = LIGHT_SPEED_KMS * mom2 / restwav

    # Optionally mask out pixels that are outside of desired ranges
    if irange is not None:
        mom0.mask = mom0.mask | (mom0.data < irange[0]) | (mom0.data > irange[1])
    if vrange is not None:
        vmean.mask = vmean.mask | (vmean.data < vrange[0]) | (vmean.data > vrange[1])
    if srange is not None:
        sigma.mask = sigma.mask | (sigma.data < srange[0]) | (sigma.data > srange[1])

    prefix = f"{flabel}-{label}-bin{rebin:02d}"
    mom0.write(SAVEPATH / f"{prefix}-sum.fits", savemask="nan", checksum=True)
    vmean.write(SAVEPATH / f"{prefix}-vmean.fits", savemask="nan", checksum=True)
    sigma.write(SAVEPATH / f"{prefix}-sigma.fits", savemask="nan", checksum=True)


def moments_corner_plot(
    moments,
    *,
    rebin=1,
    ilabel="ION",
    flabel="ion",
    label="0000",
    restwav,
    irange,
    vrange,
    srange,
    hist_bins=100,
    image_bins=50,
    hist_color="r",
    image_cmap="rocket_r",
):
    """Make corner plot of velocity moments of emission line"""

    mom0 = moments[0].rebin(rebin)
    mom1 = moments[1].rebin(rebin)
    mom2 = moments[2].rebin(rebin)

    vmean = LIGHT_SPEED_KMS * (mom1.data - restwav) / restwav
    sigma = LIGHT_SPEED_KMS * mom2.data / restwav

    m = (
        mom0.mask
        | (mom0.data < irange[0])
        | (mom0.data > irange[1])
        | (vmean < vrange[0])
        | (vmean > vrange[1])
        | (sigma < srange[0])
        | (sigma > srange[1])
    )
    df = pd.DataFrame(
        {
            f"log10 I({label})": np.log10(mom0.data[~m]),
            f"V({label})": vmean[~m],
            f"sig({label})": sigma[~m],
        }
    )
    # Add in the min/max values so plot limits are reproducible
    df0 = pd.DataFrame(
        {
            f"log10 I({label})": np.log10(irange),
            f"V({label})": np.array(vrange),
            f"sig({label})": np.array(srange),
        }
    )
    df = df.append(df0)
    g = sns.pairplot(
        df,
        kind="hist",
        height=4,
        corner=True,
        plot_kws=dict(cmap=image_cmap, bins=image_bins),
        diag_kws=dict(
            color=hist_color,
            bins=hist_bins,
            stat="count",
            common_norm=False,
        ),
    )
    #    g.axes[0, 0].set_ylim(0., None)
    g.fig.suptitle(f"{ilabel} {label} velocity moments ({rebin} x {rebin} binning)")
    g.tight_layout(pad=0)
    g.fig.savefig(FIGPATH / f"{flabel}-{label}-moments-corner-bin{rebin:02d}.pdf")
    return g
