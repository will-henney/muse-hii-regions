# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# + [markdown] pycharm={"name": "#%% md\n"}
# # Compare Peter Zeidler's NGC 346 cube with ESO data archive cube of same observations

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Imports
#
# Standard library imports

# + pycharm={"name": "#%%\n"}
from pathlib import Path
from typing import Union

# + [markdown] pycharm={"name": "#%% md\n"}
# Third-party  imports

# + pycharm={"name": "#%%\n"}
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
import seaborn as sns
import cmasher as cmr
from mpdaf.obj import Cube, Image, Spectrum
import regions as rg

sns.set_context("talk")

# + [markdown] pycharm={"name": "#%% md\n"}
# Will's libraries

# + pycharm={"name": "#%%\n"}
import wcsfile

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Paths for data files
#
# Files in `small_data_path` should not be too large to be checked into git

# + pycharm={"name": "#%%\n"}
data_path = Path.home() / "Work/Muse-Hii-Data/SMC-NGC-346"
small_data_path = Path.cwd().parent.parent / "data"

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Spectrum of the whole cube

# + pycharm={"name": "#%%\n"}
cubeP = Cube(str(data_path / "PeterZeidler" / "DATACUBE_FINAL_fwhm_cor.fits"))

# + [markdown] pycharm={"name": "#%% md\n"}
# Also load the standard ESO cube for comparison.
# The ESO cube needs the astrometry adjusting slightly to align with gaia frame.
# I originally did this by reading in a file I had made in DS9, but it is easier to just copy the Peter WCS.

# + pycharm={"name": "#%%\n"}
cubeE = Cube(str(data_path / "ADP.2017-10-16T11_04_19.247.fits"))

# + pycharm={"name": "#%%\n"}
cubeE.set_wcs(wcs=cubeP.wcs)

# + [markdown] pycharm={"name": "#%% md\n"}
# The only slight problem is that the data arrays are slightly different shapes. This does not matter so long as the lower left pixel is aligned.

# + pycharm={"name": "#%%\n"}
cubeE.wcs, cubeP.wcs

# + [markdown] pycharm={"name": "#%% md\n"}
# Read my hand-tweaked WCS that I exported from DS9

# + [markdown] pycharm={"name": "#%% md\n"}
# Extract full-cube spectra. This next cell takes several seconds to execute.

# + pycharm={"name": "#%%\n"}
specP = cubeP.sum(axis=(1, 2))
specE = cubeE.sum(axis=(1, 2))

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Global comparison of full spectrum
#
# First look at the full wav range on a log intensity scale

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 4))
specP.plot(ax=ax, label="Peter")
specE.plot(ax=ax, label="ESO")
ax.legend()
ax.set(
    yscale="log",
)
...;

# + [markdown] pycharm={"name": "#%% md\n"}
# So Peter cube is consistently brighter, presumably due to better photometry. Also has more sky features visible in the spectrum.
#
# Next, look ast the ratio between the two on a linear scale:

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 4))
(specP / specE).plot(ax=ax, label="Ratio")
ax.legend()
avratio = np.nanmedian(specP.data / specE.data)
ax.axhline(avratio, color='r', zorder=-100)
ax.set(
    ylim=[0.0, None],
)
...;

# + pycharm={"name": "#%%\n"}
avratio


# + [markdown] pycharm={"name": "#%% md\n"}
# ### Split into multiple wavelength sections
#
# If we choose 12 sections, then each section spans about 400 Angstrom.

# + pycharm={"name": "#%%\n"}
def split_spectrum(spec, n_sections):
    wave_min, wave_max = spec.wave.get_range()
    d_wave = (wave_max - wave_min) / n_sections
    start_waves = wave_min + np.arange(n_sections) * d_wave
    stop_waves = start_waves + d_wave
    return start_waves, stop_waves, d_wave


# + pycharm={"name": "#%%\n"}
wave_min, wave_max = specP.wave.get_range()
n_sections = 12
start_waves, stop_waves, d_wave = split_spectrum(cubeP, n_sections)
start_waves

# + [markdown] pycharm={"name": "#%% md\n"}
# Scale the ESO cube by the median ratio, so we can compare in more detail.

# + pycharm={"name": "#%%\n"}
spec_sections_P = [specP.subspec(w0, w0 + d_wave) for w0 in start_waves]
spec_sections_E = [avratio * specE.subspec(w0, w0 + d_wave) for w0 in start_waves]

# + [markdown] pycharm={"name": "#%% md\n"}
# First compare on a logarithmic flux scale.

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for section_P, section_E, ax in zip(spec_sections_P, spec_sections_E, axes):
    section_P.plot(ax=ax)
    section_E.plot(ax=ax)
    ax.set(yscale="log", ylabel="", xlabel="")
fig.tight_layout()

# + [markdown] pycharm={"name": "#%% md\n"}
# So, the strongest discrepancy is the [O I] 5577 line, which will be extremely weak  from  the nebula, but is strong from the night sky. This is almost completely removed in the ESO cube but is still strong in the Peter cube.

# + [markdown] pycharm={"name": "#%% md\n"}
# The Peter cube is also higher at the position of real nebula lines, such as H beta and [O III].  In these cases, I think it is not true night sky emission, but is rather the outer parts of the nebula that are being picked up. These are over-subtracted in the ESO cube, resulting in the lines being apparently seen in absorption in fainter spaxels of the cube,  The Peter cube does not have this problem.

# + [markdown] pycharm={"name": "#%% md\n"}
# Next, we look on a linear scale at the two cubes after normalizing by the median in each spectral interval.  We zoom in on the range +/-20% so we can see the weak night sky lines and atmospheric absorption.

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for section_P, section_E, ax in zip(spec_sections_P, spec_sections_E, axes):
    (section_P / np.median(section_P.data)).plot(ax=ax, label='Peter')
    (section_E / np.median(section_E.data)).plot(ax=ax, label='ESO')
    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(ylim=[0.8, 1.2], ylabel="", xlabel="")
axes[0].set_title("Peter cube versus ESO cube: summed spectrum")
axes[0].legend()
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig('peter-eso-comparison-spectra-full-cube.pdf')

# + [markdown] pycharm={"name": "#%% md\n"}
# We see a lot of weak absorption features with EW of less than an Angstrom. Also, lots of emission lines. In most sections, it is possible to distinguish the two, but at longer wavelengths this is more difficult.
#
# We could look at the model sky spectra to try and help with this.  We probably need to be able to separate the sky emission from the absorption, since the first is additive, while the second is multiplicative. **Actually, this might be a way of doing it:** in principle the absorption lines should be spatially constant in a map of `cube/avspec`, whereas the emission lines should be constant in a map of `(cube - avspec)`

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Fitting the continuum
#
# I had an idea that we can use sigma clipping in the fit so that we do not have to bother with masking out the emission lines first.  I will try this on the entire cube first, since I  am not sure how fast it is going to be to apply it pixel by pixel.
#
# Another idea is to use median filtering to get the continuum. But this might break down in regions with lots of closely spaced lines.
#
# - [ ] TODO do the fits

# + pycharm={"name": "#%%\n"}
import scipy.ndimage as ndi


# + pycharm={"name": "#%%\n"}
def get_median_continuum(spec: Spectrum, window_size=11):
    spec_c = spec.copy()
    spec_c.data = ndi.median_filter(spec.data, size=(window_size,))
    return spec_c

def get_median_continuum_cube(cube: Cube, window_size=11):
    cube_c = cube.copy()
    cube_c.data[cube_c.mask] = np.nan
    cube_c.data = ndi.median_filter(cube_c.data, size=(window_size, 1, 1))
    return cube_c



# + pycharm={"name": "#%%\n"}
specP_cont = get_median_continuum(specP, 31)
specE_cont = get_median_continuum(specE, 31)

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 4))
specP_cont.plot(ax=ax, label="Peter")
specE_cont.plot(ax=ax, label="ESO")
ax.legend()
ax.set(
    yscale="log",
)
...;

# + pycharm={"name": "#%%\n"}


# + [markdown] pycharm={"name": "#%% md\n"}
# ## Make some maps

# + pycharm={"name": "#%%\n"}
white_map_P = cubeP.sum(axis=0)
white_map_E = cubeE.sum(axis=0)
# ha_map_P = cubeP.select_lambda(6555.0, 6573.0).sum(axis=0)
# ha_map_E = cubeE.select_lambda(6555.0, 6573.0).sum(axis=0)
ha_map_P = cubeP.get_image((6555.0, 6573.0), method="sum", subtract_off=True)
ha_map_E = cubeE.get_image((6555.0, 6573.0), method="sum", subtract_off=True)

# + pycharm={"name": "#%%\n"}
ha_map_P.shape, ha_map_E.shape

# + pycharm={"name": "#%%\n"}
ha_map_E -= ha_map_E.data.min()
white_map_E -= white_map_E.data.min()

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(2, 3, sharex="all", sharey="all", figsize=(12, 8))
ha_map_P.plot(scale="sqrt", ax=axes[0, 0])
ha_map_E.plot(scale="sqrt", ax=axes[1, 0])
white_map_P.plot(scale="log", ax=axes[0, 1])
white_map_E.plot(scale="log", ax=axes[1, 1])
(ha_map_P / white_map_P).plot(scale="sqrt", cmap="gray_r", ax=axes[0, 2], vmin=0, vmax=0.15)
(ha_map_E / white_map_E).plot(scale="sqrt", cmap="gray_r", ax=axes[1, 2], vmin=0, vmax=0.4)

[ax.set(xticks=[], yticks=[]) for ax in axes.flat]
fig.tight_layout(pad=0)

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Look at particular regions
# To start with I will use the bowshock ones that I made for Jesus to use for the Spitzer spectra. Later, I should do some that concentrate more on the filaments.

# + pycharm={"name": "#%%\n"}
region_file = "ngc346-jesus-icrs.reg"
sky_regions = rg.Regions.read(small_data_path / region_file)

# + [markdown] pycharm={"name": "#%% md\n"}
# Now make a dict keyed by the name of each region. We will keep them as sky coordinates and manage the conversion to pixels at the moment of extracting the spectra.

# + pycharm={"name": "#%%\n"}
region_dict = {reg.meta["label"]: reg for reg in sky_regions}
Table(rows=region_dict.items())

# + [markdown] pycharm={"name": "#%% md\n"}
# Also read in extra regions for use later:

# + pycharm={"name": "#%%\n"}
extra_region_dict = {reg.meta["label"]: reg for reg in rg.Regions.read(small_data_path / "ngc346-extra-icrs.reg")}
Table(rows=extra_region_dict.items())

# + pycharm={"name": "#%%\n"}


# + pycharm={"name": "#%%\n"}


# + [markdown] pycharm={"name": "#%% md\n"}
# Calculate the fraction of the full spectrum flux that is in the Ha line.  We will use this as a criterion for whether a given pixel is dominated by stars or nebular emission. For example, using a threshold of 0.04.

# + pycharm={"name": "#%%\n"}
ha_frac_P = ha_map_P / white_map_P
ha_frac_E = ha_map_E / white_map_E

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(8, 8))
ha_map_P.plot(ax=ax, zscale=True, cmap="gray_r", title="Positions of sample regions")
ax.contourf(ha_frac_P.data, levels=[0.00, 0.01, 0.02, 0.03, 0.04], cmap="Greens")
for label, reg in region_dict.items():
    is_bg = label.endswith("bg")
    linestyle = "dotted" if is_bg else "solid"
    fontweight = "normal" if is_bg else "bold"
    pixreg = reg.to_pixel(ha_map_P.wcs.wcs)
    pixreg.plot(ax=ax, facecolor='none', edgecolor='red',
                lw=2,

                linestyle=linestyle,
                )
    ax.text(pixreg.center.x, pixreg.center.y, label,
            ha="center", va="center", fontsize="small", fontweight=fontweight, color="w")
...;

# + [markdown] pycharm={"name": "#%% md\n"}
# The green blobs are have a low Ha fraction and are dominated by the stellar continuum. We could maybe mask them out.

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(8, 8))
ha_map_P.plot(ax=ax, zscale=True, cmap="gray_r", title="Positions of extra sample regions")
ax.contourf(ha_frac_P.data, levels=[0.00, 0.01, 0.02, 0.03, 0.04], cmap="Greens")
for label, reg in extra_region_dict.items():
    is_bg = label.endswith("bg")
    linestyle = "dotted" if is_bg else "solid"
    fontweight = "normal" if is_bg else "bold"
    pixreg = reg.to_pixel(ha_map_P.wcs.wcs)
    pixreg.plot(ax=ax, facecolor='none', edgecolor='c',
                lw=2,
                linestyle=linestyle,
                )
    ax.text(pixreg.center.x, pixreg.center.y, label,
            ha="center", va="center", fontsize="small", fontweight=fontweight, color="w")
...;


# + pycharm={"name": "#%%\n"}
def get_spectrum_from_region(
        cube: Cube,
        region: Union[rg.PixelRegion, rg.SkyRegion],
        reduction_method: callable = np.sum,
        debug: bool = False,
) -> Spectrum:
    try:
        region_mask = region.to_mask()
    except AttributeError:
        region_mask = region.to_pixel(cube.wcs.wcs).to_mask()
    nv, ny, nx = cube.shape
    # Slices into 2D arrays
    slices_large, slices_small = region_mask.get_overlap_slices((ny, nx))
    if debug:
        print('2D slice:', slices_large)
    slices_cube = (slice(None, None),) + slices_large
    image_mask_large = region_mask.to_image((ny, nx))
    image_mask_small = image_mask_large[slices_large]
    cube_cutout = cube.data[slices_cube]
    cube_cutout[cube.mask[slices_cube]] = 0.0
    spec = reduction_method(cube_cutout * image_mask_small[None, :, :], axis=(1, 2))
    return Spectrum(wave=cube.wave, data=spec, unit=cube.unit)



# + [markdown] pycharm={"name": "#%% md\n"}
# ### Separate the nebula from the stars

# + pycharm={"name": "#%%\n"}
mask_stars_P = ha_frac_P.data < 0.04
mask_stars_E = ha_frac_E.data < 0.04

# + [markdown] pycharm={"name": "#%% md\n"}
# This defines "star" pixels as being those where the net continuum-subtracted Ha line flux is less than 0.04 times the bolometric flux of the cube.
#
# We can probably convert this into an equivalent width. Something like:

# + pycharm={"name": "#%%\n"}
(wave_max - wave_min) * 0.04

# + [markdown] pycharm={"name": "#%% md\n"}
# So that is 190 Angstroms, but this is the equivalent width with respect to the average continuum over the full MUSE range.

# + [markdown] pycharm={"name": "#%% md\n"}
# #### What to do about the bright YSO?
#
# This has a spectrum that is very different from the more diffuse gas emission, and it is a significant fraction of the total flux in some lines. So, it is probably best to leave it out of the nebular average, so that YSO lines do not show up as absorption in the ratio spectra of the other regions.

# + pycharm={"name": "#%%\n"}
mask_yso_P = region_dict['YSO'].to_pixel(ha_map_P.wcs.wcs).to_mask().to_image(ha_map_P.data.shape).astype(bool)
mask_yso_E = region_dict['YSO'].to_pixel(ha_map_E.wcs.wcs).to_mask().to_image(ha_map_E.data.shape).astype(bool)
mask_yso_P

# + [markdown] pycharm={"name": "#%% md\n"}
# #### Save the masks to fits files
#

# + pycharm={"name": "#%%\n"}
Image(data=mask_stars_E.astype(int), wcs=ha_frac_E.wcs).write(
    small_data_path / "n346-mask-stars.fits",
    savemask="none",
)
Image(data=mask_yso_E.astype(int), wcs=ha_frac_E.wcs).write(
    small_data_path / "n346-mask-yso.fits",
    savemask="none",
)

# + [markdown] pycharm={"name": "#%% md\n"}
# #### Calculate mean spectra for stars and for diffuse nebula
#
# Remember that the cube mask is `True` for invalid pixels

# + pycharm={"name": "#%%\n"}
cube_P_neb = cubeP.copy()
cube_P_neb.mask = cube_P_neb.mask | mask_stars_P | mask_yso_P
spec_mean_P_neb = cube_P_neb.mean(axis=(1, 2))

# + pycharm={"name": "#%%\n"}
cube_P_yso = cubeP.copy()
cube_P_yso.mask = cube_P_yso.mask | mask_stars_P | (~mask_yso_P)
spec_mean_P_yso = cube_P_yso.mean(axis=(1, 2))
del cube_P_yso

# + pycharm={"name": "#%%\n"}
cube_P_stars = cubeP.copy()
cube_P_stars.mask = cube_P_stars.mask | (~mask_stars_P)
spec_mean_P_stars = cube_P_stars.mean(axis=(1, 2))
del cube_P_stars

# + pycharm={"name": "#%%\n"}
cube_E_neb = cubeE.copy()
cube_E_neb.mask = cube_E_neb.mask | mask_stars_E | mask_yso_E
spec_mean_E_neb = cube_E_neb.mean(axis=(1, 2))

# + pycharm={"name": "#%%\n"}
cube_E_yso = cubeE.copy()
cube_E_yso.mask = cube_E_yso.mask | mask_stars_E | (~mask_yso_E)
spec_mean_E_yso = cube_E_yso.mean(axis=(1, 2))
del cube_E_yso

# + pycharm={"name": "#%%\n"}
cube_E_stars = cubeE.copy()
cube_E_stars.mask = cube_E_stars.mask | (~mask_stars_E)
spec_mean_E_stars = cube_E_stars.mean(axis=(1, 2))
del cube_E_stars

# + [markdown] pycharm={"name": "#%% md\n"}
# I have deleted all the cubes except the nebula ones after use in order not to use too much memory. Previously I was getting past 64GB,  which is not too healthy. Now it is at only (!!) 38GB, which fits comfortably in RAM on my laptop

# + [markdown] pycharm={"name": "#%% md\n"}
# Now repeat the graphs that compare Peter with ESO, but this time for the nebula and stars separately

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for section_P, section_E, ax in zip(
        [spec_mean_P_neb.subspec(w0, w0 + d_wave) for w0 in start_waves],
        [avratio * spec_mean_E_neb.subspec(w0, w0 + d_wave) for w0 in start_waves],
        axes
):
    (section_P / np.median(section_P.data)).plot(ax=ax, label='Peter')
    (section_E / np.median(section_E.data)).plot(ax=ax, label='ESO')
    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(ylim=[0.8, 1.2], ylabel="", xlabel="")
axes[0].set_title("Peter cube versus ESO cube: Nebular spectrum")
axes[0].legend()
fig.tight_layout()
fig.savefig('peter-eso-comparison-spectra-nebula-mean.pdf')

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for section_P, section_E, ax in zip(
        [spec_mean_P_stars.subspec(w0, w0 + d_wave) for w0 in start_waves],
        [avratio * spec_mean_E_stars.subspec(w0, w0 + d_wave) for w0 in start_waves],
        axes
):
    (section_P / np.median(section_P.data)).plot(ax=ax, label='Peter')
    (section_E / np.median(section_E.data)).plot(ax=ax, label='ESO')
    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(ylim=[0.8, 1.2], ylabel="", xlabel="")
axes[0].set_title("Peter cube versus ESO cube: Star spectrum")
axes[0].legend()
fig.tight_layout()
fig.savefig('peter-eso-comparison-spectra-star-mean.pdf')

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for section_P, section_E, ax in zip(
        [spec_mean_P_yso.subspec(w0, w0 + d_wave) for w0 in start_waves],
        [avratio * spec_mean_E_yso.subspec(w0, w0 + d_wave) for w0 in start_waves],
        axes
):
    (section_P / np.median(section_P.data)).plot(ax=ax, label='Peter')
    (section_E / np.median(section_E.data)).plot(ax=ax, label='ESO')
    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', linewidth=0.5, alpha=0.5)
    ax.set(ylim=[0.8, 1.2], ylabel="", xlabel="")
axes[0].set_title("Peter cube versus ESO cube: YSO spectrum")
axes[0].legend()
fig.tight_layout()
fig.savefig('peter-eso-comparison-spectra-yso-mean.pdf')

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Test the function for obtaining 1D spectrum from region
#
# Use the bow shock region as an example

# + pycharm={"name": "#%%\n"}
specBS = get_spectrum_from_region(cube_P_neb, region_dict['BS'], reduction_method=np.nanmean, debug=True)
specBSbg = get_spectrum_from_region(cube_P_neb, region_dict['BS bg'], reduction_method=np.nanmean, debug=True)


# + [markdown] pycharm={"name": "#%% md\n"}
# Compare with the slow way of doing it

# + pycharm={"name": "#%%\n"}

def get_spectrum_from_region_slow(
        cube: Cube,
        region: Union[rg.PixelRegion, rg.SkyRegion],
) -> Spectrum:
    """Simpler but slower algorithm that simply multiplies cube by region mask and then does the reduction"""
    try:
        region_mask = region.to_mask()
    except AttributeError:
        region_mask = region.to_pixel(cube.wcs.wcs).to_mask()

    nv, ny, nx = cube.shape
    image_mask_large = region_mask.to_image((ny, nx))
    return (cube * image_mask_large).mean(axis=(1, 2)) * (ny * nx) / np.sum(image_mask_large)



# + pycharm={"name": "#%%\n"}
specBS_2 = get_spectrum_from_region_slow(cube_P_neb, region_dict['BS'])
# specBSbg_2 = get_spectrum_from_region_slow(cube_P_neb, region_dict['BS bg'])

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots()
(specBS_2 / specBS).plot()
ax.set(ylim=[0.85, 1.15])

# + [markdown] pycharm={"name": "#%% md\n"}
# Strangely, the ratio is a bit larger than unity. Also, it shows a very slight reduction in the blue

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
specBS.plot()
specBSbg.plot()
(specBS - specBSbg).plot()
ax.set(yscale='log', ylim=[1e0, 1e5])

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
(specBS - specBSbg).plot()
ax.set(yscale='linear', ylim=[-10, 100], xlim=[4600, 5200])

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
(specBS / specBSbg).plot()
ax.set(yscale='linear', ylim=[0.95, 1.1], xlim=[4600, 5200])

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
(specBS / spec_mean_P_neb).plot()
(specBSbg / spec_mean_P_neb).plot()

ax.set(yscale='linear', ylim=[0.7, 1.1], xlim=[7100, 7400])

# + [markdown] pycharm={"name": "#%% md\n"}
# So dividing by the mean nebula does seem  to work after all.  We will try it again below

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Try doing the median continuum

# + pycharm={"name": "#%%\n"}
specBS_cont = get_median_continuum(specBS, 31)
specBSbg_cont = get_median_continuum(specBSbg, 31)

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
specBS.plot()
specBSbg.plot()
specBS_cont.plot()
specBSbg_cont.plot()

ax.set(yscale='linear', ylim=[0, 500])

# + pycharm={"name": "#%%\n"}


# + pycharm={"name": "#%%\n"}


# + [markdown] pycharm={"name": "#%% md\n"}
# ### Repeat for the ESO spectrum

# + pycharm={"name": "#%%\n"}
specBS_E = get_spectrum_from_region(cube_E_neb, region_dict['BS'], reduction_method=np.nanmean, debug=True)
specBSbg_E = get_spectrum_from_region(cube_E_neb, region_dict['BS bg'], reduction_method=np.nanmean, debug=True)

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
specBS_E.plot(label='ESO BS')
specBSbg_E.plot(label='ESO BS bg')
specBS.plot(label='Peter BS')
specBSbg.plot(label='Peter BS bg')
# (specBS_E - specBSbg_E).plot()
ax.legend()
ax.set(yscale='log', ylim=[1e0, 1e5])

# + pycharm={"name": "#%%\n"}


# + pycharm={"name": "#%%\n"}
specBS_E_cont = get_median_continuum(specBS_E, 31)
specBSbg_E_cont = get_median_continuum(specBSbg_E, 31)

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
specBS_E.plot()
specBSbg_E.plot()
specBS_E_cont.plot()
specBSbg_E_cont.plot()

ax.set(yscale='linear', ylim=[0, 100])

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
((specBS_E - specBS_E_cont) - (specBSbg_E - specBSbg_E_cont)).plot()
((specBS - specBS_cont) - (specBSbg - specBSbg_cont) + 60).plot()
((specBSbg_E - specBSbg_E_cont) - 10).plot()
((specBSbg - specBSbg_cont) + 40).plot()
ax.set(yscale='linear', ylim=[-20, 100], xlim=[4600, 5200])

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
((specBS_E - specBS_E_cont) - (specBSbg_E - specBSbg_E_cont)).plot()
((specBS - specBS_cont) - (specBSbg - specBSbg_cont) + 60).plot()
((specBSbg_E - specBSbg_E_cont) - 10).plot()
((specBSbg - specBSbg_cont) + 40).plot()
ax.set(yscale='linear', ylim=[-20, 100], xlim=[6100, 6900])

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 6))
((specBS_E - specBS_E_cont) - (specBSbg_E - specBSbg_E_cont)).plot()
((specBS - specBS_cont) - (specBSbg - specBSbg_cont) + 60).plot()
((specBSbg_E - specBSbg_E_cont) - 10).plot()
((specBSbg - specBSbg_cont) + 40).plot()
ax.minorticks_on()
ax.set(yscale='linear', ylim=[-20, 100], xlim=[7100, 7400])

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Extract the spectra from all the regions

# + pycharm={"name": "#%%\n"}
spec_dict_P_neb = {
    label: get_spectrum_from_region(cube_P_neb, reg, reduction_method=np.nanmean)
    for label, reg in region_dict.items()
}
spec_dict_extra_P_neb = {
    label: get_spectrum_from_region(cube_P_neb, reg, reduction_method=np.nanmean)
    for label, reg in extra_region_dict.items()
}

# + pycharm={"name": "#%%\n"}
spec_dict_E_neb = {
    label: get_spectrum_from_region(cube_E_neb, reg, reduction_method=np.nanmean)
    for label, reg in region_dict.items()
}
spec_dict_extra_E_neb = {
    label: get_spectrum_from_region(cube_E_neb, reg, reduction_method=np.nanmean)
    for label, reg in extra_region_dict.items()
}

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 8))
for label, spec in spec_dict_P_neb.items():
    if label.endswith(" bg"):
        lw = 1
    else:
        lw = 2
    ratio = spec / spec_mean_P_neb
    norm = np.median(ratio.data)
    (ratio / norm).plot(label=label, linewidth=lw)
ax.legend(ncol=5)
ax.set(ylim=[0.5, 1.5])
...;

# + pycharm={"name": "#%%\n"}
colors = cmr.take_cmap_colors(
    'cmr.ember', len(spec_dict_extra_P_neb) // 2, cmap_range=(0.25, 0.95), return_fmt='hex')

fig, ax = plt.subplots(figsize=(12, 8))
for i, (label, spec) in enumerate(spec_dict_extra_P_neb.items()):
    if label.endswith(" bg"):
        lw = 1
    else:
        lw = 4
    ratio = spec / spec_mean_P_neb
    norm = np.median(ratio.data)
    (ratio).plot(label=label, linewidth=lw, color=colors[i // 2], alpha=0.6)
    # (ratio / norm).plot(label=label, linewidth=lw)
ax.legend(ncol=4)
ax.set(ylim=[0., 2.5])
...;

# + pycharm={"name": "#%%\n"}
source_regions = set(s.split()[0] for s in region_dict.keys())
source_regions

# + pycharm={"name": "#%%\n"}
source_regions_extra = set(s.split()[0] for s in extra_region_dict.keys())
source_regions_extra

# + [markdown] pycharm={"name": "#%% md\n"}
# Combining the two sets uses the `|` operator, rather than `+`

# + pycharm={"name": "#%%\n"}
source_regions | source_regions_extra

# + [markdown] pycharm={"name": "#%% md\n"}
# Nominal systemic velocity in km/s

# + pycharm={"name": "#%%\n"}
v_sys = 170.0

# + [markdown] pycharm={"name": "#%% md\n"}
# #### Plot per-region raw spectra for Peter cube

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    section_mean_stars = spec_mean_P_stars.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_P_neb[label]  #- spec_dict_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg  #/ section_mean
        norm = np.median(ratio.data)
        ((ratio / norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean / norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
    # Mean spectrum of the stars
    norm = np.median(section_mean_stars.data)
    ((section_mean_stars / norm) - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.3)
    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    ax.plot([w1, w2], [0.7, 0.7], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, 0.7, f"  {v_sys:.0f} km/s", ha='left', va='center')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[0.6, 1.7],
        # yscale="log",
        ylabel="", xlabel="")
axes[0].set_title("Peter cube, mean raw spectrum by nebular region, plus mean spectrum of gas and of stars")
axes[0].legend(ncol=7, fontsize="x-small", loc="lower left")
axes[-1].set(
    xlabel='Observed wavelength, Angstrom',
)
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-region-strip-spectra.pdf")

# + [markdown] pycharm={"name": "#%% md\n"}
# This shows lots of absorption features, which I do not think are real. Later, we try dividing by  the mean spectrum, which works very well for eliminating the majority of these,  but it would still be good to understand them and see if they can  be removed cleanly.
#
# Here are some observations and ideas:
#
# 1. In the blue region, the strongest observed absorptions seem to coincide with telluric atmospheric absorption lines that I see in sky models. For  instance, the trio at 5166, 5172, 5183. This suggests an under-correction for the telluric absorption.
# 2. In the redder regions, the observed absorption features correspond closely to *emission* lines in the night sky spectrum. For example at 7240, 7246. Or 7401, or 7523. This suggests an *over*-correction for the night sky emission.
# 3. There are also some apparent observed emission features, such as the hedgehog pattern at 7600–7670, which corresponds to the O_2 absorption A-band. This suggests an over-correction for telluric absorption.
#

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    section_mean_stars = spec_mean_P_stars.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions_extra:
        spec = spec_dict_extra_P_neb[label]  #- spec_dict_extra_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg  #/ section_mean
        norm = np.median(ratio.data)
        ((ratio / norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.05
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean / norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
    # Mean spectrum of the stars
    # norm = np.median(section_mean_stars.data)
    # ((section_mean_stars/norm) - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.3)
    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    ax.plot([w1, w2], [0.7, 0.7], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, 0.7, f"  {v_sys:.0f} km/s", ha='left', va='center')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[0.6, 1.7],
        # yscale="log",
        ylabel="", xlabel="")
axes[0].set_title("Peter cube, mean raw spectrum by nebular region, plus mean spectrum of gas and of stars")
axes[0].legend(ncol=5, fontsize="x-small", loc="lower left")
axes[-1].set(
    xlabel='Observed wavelength, Angstrom',
)
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-region-extras-strip-spectra.pdf")

# + [markdown] pycharm={"name": "#%% md\n"}
#

# + [markdown] pycharm={"name": "#%% md\n"}
# And do it again but for all regions except the YSO and with the lines on top of each other

# + pycharm={"name": "#%%\n"}
all_region_dict = {**region_dict, **extra_region_dict}
all_spec_dict_P_neb = {**spec_dict_P_neb, **spec_dict_extra_P_neb}

# + pycharm={"name": "#%%\n"}
all_source_regions = (source_regions | source_regions_extra) - {'YSO'}
all_source_regions

# + pycharm={"name": "#%%\n"}
region_types = set(_.split('-')[0] for _ in all_source_regions)
region_types

# + pycharm={"name": "#%%\n"}
color_dict = {'BS': "r", 'FIL': 'g', 'GLOB': 'b', 'MIP': 'm', 'NEUT': 'c'}

# + pycharm={"name": "#%%\n"}
n_col, n_row = 5, 12
start_waves, stop_waves, d_wave = split_spectrum(cubeP, n_col * n_row)

fig, axes = plt.subplots(n_row, n_col, figsize=(12, 2 * n_row), sharey="all")
for w0, w1, ax in zip(start_waves, stop_waves, axes.flat):
    section_mean = spec_mean_P_neb.subspec(w0, w1)
    section_mean_stars = spec_mean_P_stars.subspec(w0, w1)
    spec_list = []
    for label in all_source_regions:
        spec = all_spec_dict_P_neb[label]  #- all_spec_dict_P_neb[label + " bg"]
        sub_spec = spec.subspec(w0, w1)
        norm = np.median(sub_spec.data)
        color = color_dict[label.split('-')[0]]
        spec_list.append(sub_spec.data.data / norm)
        (sub_spec / norm).plot(ax=ax, label=label, linewidth=0.5, alpha=0.6, color=color)
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean / norm)).plot(ax=ax, label="mean nebula", linewidth=0.5, color="k")
    # Find std deviation between regions
    spec_stack = np.stack(spec_list)
    spec_std = np.nanstd(spec_stack, axis=0)
    spec_mean = np.nanmean(spec_stack, axis=0)
    ax.fill_between(section_mean.wave.coord(), 0.9 + spec_std / spec_mean, 0.9, color='y')

    # # Indicator for rest-to-observed wavelength transformation
    ww1 = w0 + 0.9 * d_wave
    ww2 = ww1 * (1.0 + v_sys / 3e5)
    ax.plot([ww1, ww2], [0.95, 0.95], color='m', lw=5, solid_capstyle='butt')
    # ax.text(w2, 0.7, f"  {v_sys:.0f} km/s", ha='left', va='center')

    ax.xaxis.set_major_locator(MaxNLocator(min_n_ticks=1, nbins=2))
    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[0.85, 1.15],
        # yscale="log",
        ylabel="", xlabel="")
axes[0, n_col // 2].set_title("Peter cube, mean raw spectrum by nebular region, plus mean spectrum of gas")
# axes[0].legend(ncol=5, fontsize="x-small", loc="lower left")
axes[-1, 0].set(
    xlabel='Observed wavelength, Angstrom',
)
fig.tight_layout(w_pad=0.05)

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-all-region-strip-spectra.pdf")

# + pycharm={"name": "#%%\n"}
np.info(section_mean.wave.coord)

# + pycharm={"name": "#%%\n"}
locator = MaxNLocator(min_n_ticks=1, nbins=2)

# + pycharm={"name": "#%%\n"}
for w0, w1 in zip(start_waves, stop_waves):
    print(locator.tick_values(w0, w1))

# + [markdown] pycharm={"name": "#%% md\n"}
# #### Plot per-region raw spectra for ESO cube

# + pycharm={"name": "#%%\n"}
start_waves, stop_waves, d_wave = split_spectrum(cubeE, n_sections)

fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_E_neb.subspec(w0, w0 + d_wave)
    section_mean_stars = spec_mean_E_stars.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_E_neb[label]  #- spec_dict_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg  #/ section_mean
        norm = np.median(ratio.data)
        ((ratio / norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean / norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
    # Mean spectrum of the stars
    norm = np.median(section_mean_stars.data)
    ((section_mean_stars / norm) - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.3)
    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    ax.plot([w1, w2], [0.7, 0.7], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, 0.7, f"  {v_sys:.0f} km/s", ha='left', va='center')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[0.6, 1.7],
        # yscale="log",
        ylabel="", xlabel="")
axes[0].set_title("ESO cube, mean raw spectrum by nebular region, plus mean spectrum of gas and of stars")
axes[0].legend(ncol=7, fontsize="large", loc="lower left")
axes[-1].set(
    xlabel='Observed wavelength, Angstrom',
)
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("eso-region-strip-spectra.pdf")

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Continuum-subtracted and BG-subtracted ESO spectra by region

# + [markdown] pycharm={"name": "#%% md\n"}
# Repeat the strip plots that we did previously, except subtracting the median continuum and subtracting the BG region from each source region.  This seems to work much better for the ESO cube than for the Peter cube, at least for the weaker lines

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Net BG-subtracted spectrum for each source region

# + pycharm={"name": "#%%\n"}
net_spec_dict_E = {
    label: spec_dict_E_neb[label] - spec_dict_E_neb[label + ' bg']
    for label in source_regions
}
net_spec_dict_E['YSO'] = spec_mean_E_yso
net_spec_dict_E |= {
    label: spec_dict_extra_E_neb[label] - spec_dict_extra_E_neb[label + ' bg']
    for label in source_regions_extra
}
net_spec_dict_E |= {'ALL': spec_mean_E_neb, 'STARS': spec_mean_E_stars}

# + [markdown] pycharm={"name": "#%% md\n"}
# We have combined all the spatial regions from the original Jesús group and the extras. TODO Something went wrong with the BG subtraction of the YSO spectrum, so I restored it to the original without subtracting BG.

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Median-filtered continuum
#
# Next, calculate continuum with median filter along wavelength axis. Set the width to 101 pixels   so that we smooth over even the bright lines.

# + pycharm={"name": "#%%\n"}
n_wav_filter_width = 101
n_wav_filter_width * cubeE.wave.get_step()

# + pycharm={"name": "#%%\n"}
net_cont_dict_E = {
    label: get_median_continuum(spec, n_wav_filter_width)
    for label, spec in net_spec_dict_E.items()
}

# + pycharm={"name": "#%%\n"}
net_cont_dict_E

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Continuum-subtracted spectra

# + pycharm={"name": "#%%\n"}
net_csub_dict_E = {
    label: net_spec_dict_E[label] - net_cont_dict_E[label]
    for label in net_spec_dict_E
}


# + [markdown] pycharm={"name": "#%% md\n"}
# We now have dictionaries of the net BG-subtracted spectrum of each region and of the median-filtered continuum of each of those, and also the continuum-subtracted spectra (csub).

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Types of source region
#
# Make function to get the first part of the region label, which specifies its type. Note that the distinction between `FIL` and `GLOB` is a bit debatable.

# + pycharm={"name": "#%%\n"}
def get_region_type(label):
    return label.split('-')[0]


# + [markdown] pycharm={"name": "#%% md\n"}
# Test it by getting the set of all unique region types.

# + pycharm={"name": "#%%\n"}
set(map(get_region_type, net_cont_dict_E))

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Division of full wavelength range into sections
#
# Re-calculate wave range division for 12 sections

# + pycharm={"name": "#%%\n"}
start_waves, stop_waves, d_wave = split_spectrum(cubeE, n_sections)


# + [markdown] pycharm={"name": "#%% md\n"}
# ### Normalization of the spectra.
#
# ~~Normalization of spectrum for each type of region.  The idea was to  get the Paschen lines at similar levels.~~ I do it differently now, with a separate normalization for each region rather than having all of the same type being the same.

# + pycharm={"name": "#%%\n"}
# norm_by_region = {
#     'BS': 100.0,
#     'FIL': 10.0,
#     'GLOB': 50.0,
#     'MIP': 60.0,
#     'NEUT': 30.0,
#     'YSO': 1000.0,
#     'ALL': 80.0,
#     'STARS': 200.0
# }

# + pycharm={"name": "#%%\n"}
def norm_by_lines(
        spectrum: Spectrum,
        rest_waves: list,
        wave_width: float = 5.0,
        v_shift: float = 170.0,
        stat_function: callable=np.nanmax,
) -> float:
    """
    Normalize a spectrum by the strongest of a list of emission lines

    :param spectrum: One-dimensional continuum-subtracted spectrum
    :param rest_waves: List of rest wavelengths of the emission lines to normalize by
    :param wave_width: Full width of wave window for line extraction
    :param v_shift: Systemic velocity for shifting the rest wave
    :param stat_function: Function to calculate statistics on wave window, e.g., np.nanmax, np,nansum, etc
    :return: Normalization factor `norm` such that `spectrum` / `norm` has a peak value of unity for the strongest emission line in the `rest_waves` list
    """
    peaks = []
    for wave in rest_waves:
        _wave = wave * (1 + v_shift / 300000)
        _spectrum = spectrum.subspec(_wave - wave_width / 2, _wave + wave_width / 2)
        peaks.append(stat_function(_spectrum.data))
    return max(peaks)



# + [markdown] pycharm={"name": "#%% md\n"}
# Normalize to the strongest of the [Cl III] and [Cl IV] lines. I chose these, since one or both is clearly visible in all of the regions, plus they are not affected by underlying stellar absorption like the H and He lines are. I take 10 times the normalization so that these lines have a peak of about 0.1, which is the vertical offset between adjacent spectra in the graph.
#
# Also, I make manual adjustments to the STARS normalization, and also to FIL-C since it has a sawtooth pattern due to the velocity mismatch with its background.
#
# Also, now do the same to NEUT-D , which has a similar behavior to FIL-C, although I am not sure why

# + pycharm={"name": "#%%\n"}
norm_by_region = {
    label: 10 * norm_by_lines(spec, [5518, 8046])
    for label, spec in net_csub_dict_E.items()
}
norm_by_region['STARS'] *= 5
norm_by_region['FIL-C'] *= 2
norm_by_region['NEUT-D'] *= 4
norm_by_region

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Strip graphs of the spectra for different groups of regions
#
# Make graphs by wave section of the original regions. Each spectrum is BG-subtracted and with the median-filtered continuum subtracted.

# + [markdown] pycharm={"name": "#%% md\n"}
# #### The original Jesús group

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    shift = 0.0
    for label in sorted(source_regions):
        spec = net_spec_dict_E[label]
        cont = net_cont_dict_E[label]
        section_reg = (spec - cont).subspec(w0, w0 + d_wave)
        norm = norm_by_region[label]
        (section_reg / norm + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    section_mean = (net_spec_dict_E['ALL'] - net_cont_dict_E['ALL']).subspec(w0, w0 + d_wave)
    norm = norm_by_region['ALL']
    (section_mean / norm + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")

    # Mean spectrum of the stars
    section_mean_stars = (net_spec_dict_E['STARS'] - net_cont_dict_E['STARS']).subspec(w0, w0 + d_wave)
    norm = norm_by_region['STARS']
    (section_mean_stars / norm - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.2)
    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    y0 = -0.1
    ax.plot([w1, w2], [y0, y0], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, y0, f"  {v_sys:.0f} km/s", ha='left', va='center')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[-0.15, 0.65],
        # yscale="log",
        ylabel="", xlabel="")
fig.suptitle(
    "ESO cube, continuum-subtracted spectrum by nebular region, plus mean spectrum of gas and of stars",
    fontsize='x-large',
    y=0.99,
)
axes[0].legend(ncol=7, fontsize="large", bbox_to_anchor=(1.0, 1.0), loc='lower right')
axes[-1].set(
    xlabel='Observed wavelength, Angstrom',
)
fig.tight_layout(rect=(0, 0, 1, 0.98))

# + pycharm={"name": "#%%\n"}
fig.savefig("eso-region-strip-csub-spectra.pdf")

# + [markdown] pycharm={"name": "#%% md\n"}
# #### The extra regions group
#
# Repeat for the extra regions. There are more of these, so the lines are more squashed together

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    shift = 0.0
    for label in sorted(source_regions_extra):
        spec = net_spec_dict_E[label]
        cont = net_cont_dict_E[label]
        section_reg = (spec - cont).subspec(w0, w0 + d_wave)
        norm = norm_by_region[label]
        (section_reg / norm + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    section_mean = (net_spec_dict_E['ALL'] - net_cont_dict_E['ALL']).subspec(w0, w0 + d_wave)
    norm = norm_by_region['ALL']
    (section_mean / norm + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")

    # Mean spectrum of the stars
    section_mean_stars = (net_spec_dict_E['STARS'] - net_cont_dict_E['STARS']).subspec(w0, w0 + d_wave)
    norm = norm_by_region['STARS']
    (section_mean_stars / norm - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.2)
    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    y0 = -0.1
    ax.plot([w1, w2], [y0, y0], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, y0, f"  {v_sys:.0f} km/s", ha='left', va='center')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[-0.15, shift + 0.05],
        ylabel="", xlabel="")
fig.suptitle(
    "ESO cube, continuum-subtracted spectrum by nebular region, plus mean spectrum of gas and of stars",
    fontsize='x-large',
    y=0.99,
)
axes[0].legend(ncol=7, fontsize="large", bbox_to_anchor=(1.0, 1.0), loc='lower right')
axes[-1].set(
    xlabel='Observed wavelength, Angstrom',
)
fig.tight_layout(rect=(0, 0, 1, 0.98))

# + pycharm={"name": "#%%\n"}
fig.savefig("eso-region-extras-strip-csub-spectra.pdf")


# + [markdown] pycharm={"name": "#%% md\n"}
# ### Graph combining all the regions, but in a logical order
#
# We can try and order them according to different line ratios

# + [markdown] pycharm={"name": "#%% md\n"}
# #### Ionization ratios

# + pycharm={"name": "#%%\n"}
def get_line_ratio_dict(
        spectra: dict[str, Spectrum],
        waves_1: list,
        waves_2: list,
        stat_function: callable=np.nanmax,
) -> dict[str, float]:
    return {
        _label: norm_by_lines(_spec, waves_1, stat_function=stat_function) / norm_by_lines(_spec, waves_2, stat_function=stat_function)
        for _label, _spec in spectra.items()
    }


# + [markdown] pycharm={"name": "#%% md\n"}
# Make some tables of various line ratios.  It turns out that pandas works better than astropy.Table for displaying the table in Dataspell, but ymmv in other editors.

# + pycharm={"name": "#%%\n"}
import pandas as pd
df = pd.DataFrame(
    {
        'ariv/ariii': get_line_ratio_dict(net_csub_dict_E, [4740], [7136]),
        'oiii/oii': get_line_ratio_dict(net_csub_dict_E, [4959], [7318]),
        'siii/sii': get_line_ratio_dict(net_csub_dict_E, [9069], [6731]),
        'oii/oi': get_line_ratio_dict(net_csub_dict_E, [7318], [8446]),
        'oi/fei': get_line_ratio_dict(net_csub_dict_E, [8446], [8151]),
        'ni/ci': get_line_ratio_dict(net_csub_dict_E, [8703, 8712], [8727]),
        'fei/ha': get_line_ratio_dict(net_csub_dict_E, [8151], [6563]),
    }
)
df

# + pycharm={"name": "#%%\n"}
df = pd.DataFrame(
    {
       'ha/hb': get_line_ratio_dict(net_csub_dict_E, [6563], [4861]),
       'nii': get_line_ratio_dict(net_csub_dict_E, [6583], [6548]),
       'oiii': get_line_ratio_dict(net_csub_dict_E, [5007], [4959]),
       'ariii': get_line_ratio_dict(net_csub_dict_E, [7751], [7136]),
       'oii': get_line_ratio_dict(net_csub_dict_E, [7318], [7330]),
    }
)
df

# + [markdown] pycharm={"name": "#%% md\n"}
# These are ratios that should be constant, except for extinction in the ha/hb case.
#
# Now I calculate ratios with respect to H beta. Also, add a column with the region type. This uses techniques I learned from Chapter 11 of Matt Harrison's _Effective Pandas_ book

# + pycharm={"name": "#%%\n"}
line_list = [
    '[O III] 5007',
    '[O II] 7318',
    '[O I] 6300',
    'O I 8446',
    'Fe I 8151',
    'C I 8727',
    'N I 8223',
    'He II 4686',
    '[Ar IV] 4740',
    '[Ar III] 7136',
    '[S III] 9069',
    '[S II] 6731',
]
df = pd.DataFrame(
    {
        label: get_line_ratio_dict(net_csub_dict_E, [float(label.split()[-1])], [4861])
        for label in line_list
    }
).assign(
    reg_type=lambda x: x.index.to_series()
                       .astype('string').str
                       .split('-', expand=True)
                       .iloc[:, 0].astype('category')
)
df

# + pycharm={"name": "#%%\n"}
g = sns.pairplot(df, hue="reg_type", diag_kind='hist')

# + pycharm={"name": "#%%\n"}
g.savefig("eso-region-line-ratio-pairplot.pdf")

# + pycharm={"name": "#%%\n"}
line_list_low = [
    '[S II] 6731',
    '[O I] 6300',
    'O I 8446',
    'N I 8223',
    'Fe I 8151',
    'C I 8727',
]

# + pycharm={"name": "#%%\n"}


# + pycharm={"name": "#%%\n"}
_min, _max = 1e-5, 2.0
N = 10
bins = np.logspace(np.log10(_min), np.log10(_max), N)
def update_offdiag(xdata, ydata, **kwds):
    ax = plt.gca()
    ax.set(
        xscale='log',
        yscale='log',
        xlim=[_min, _max],
        ylim=[_min, _max],
        xticks=[1e-4, 1e-2, 1e0],
        yticks=[1e-4, 1e-2, 1e0],
    )

g = sns.pairplot(
    df,
    vars=line_list_low,
    hue="reg_type",
    diag_kind='hist',
    diag_kws=dict(bins=bins),
)
g.diag_sharey = False
g.map_offdiag(update_offdiag)
...;

# + pycharm={"name": "#%%\n"}
g.savefig("eso-region-line-ratio-low-pairplot.pdf")

# + [markdown] pycharm={"name": "#%% md\n"}
# So from this it looks like forbidden oi and sii form a nice consistent sequence, although FIL-C is a bit of an outlier, presumably due to oversubtraction of the background, which is caused by the filament being foreground.

# + pycharm={"name": "#%%\n"}
df.corr()

# + pycharm={"name": "#%%\n"}
df[line_list_low].corr()

# + [markdown] pycharm={"name": "#%% md\n"}
# So, from the correlation table, the low ionization lines fall naturally into three pairs:
# - [S II] and [O I], with  r = 0.97
# - O I and N I with r = 0.93
# - Fe I and C I with r = 0.76
#
# With other correlations being significantly worse
#

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Better version of the plot with all the regions
#

# + pycharm={"name": "#%%\n"}


# + pycharm={"name": "#%%\n"}
all_source_regions = sorted(source_regions | source_regions_extra)
all_source_regions

# + pycharm={"name": "#%%\n"}
region_types = set(_.split('-')[0] for _ in all_source_regions)
region_types

# + pycharm={"name": "#%%\n"}
region_colors = {
    'BS': 'light:#f00',
    'MIP': 'light:#c60',
    'NEUT': 'light:#088',
    'GLOB': 'light:#660',
    'FIL': 'light:#080',
    'YSO': 'light:#808',
}
skip_colors = 3

# + pycharm={"name": "#%%\n"}
region_order = {
    'BS': [''],
    'MIP': ['', '-B'],
    'NEUT': ['-C', '-D', '-B', ''],
    'GLOB': ['-B', '-F', '-D', '-G', '-C', '-E', ''],
    'FIL': ['', '-D', '-C', '-B'],
    'YSO': [''],
}

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 5 * n_sections))
for w0, ax in zip(start_waves, axes):
    shift = 0.0
    for region_type, suffixes in region_order.items():
        colors = sns.color_palette(region_colors[region_type], n_colors=skip_colors + len(suffixes), desat=1.0)
        for suffix, color in zip(suffixes, colors[skip_colors:]):
            label = region_type + suffix
            spec = net_spec_dict_E[label]
            cont = net_cont_dict_E[label]
            section_reg = (spec - cont).subspec(w0, w0 + d_wave)
            norm = norm_by_region[label]
            (section_reg / norm + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6, color=color)
            ax.text(w0 + d_wave, shift, f" {label}", color=color, ha='left', va='center')
            shift += 0.1

    # Mean spectrum of the nebula
    section_mean = (net_spec_dict_E['ALL'] - net_cont_dict_E['ALL']).subspec(w0, w0 + d_wave)
    norm = norm_by_region['ALL']
    color = 'k'
    (section_mean / norm + shift).plot(ax=ax, label="mean nebula", linewidth=3, color=color)
    ax.text(w0 + d_wave, shift, " ALL", color=color, ha='left', va='center', fontweight='bold')

    # Mean spectrum of the stars
    section_mean_stars = (net_spec_dict_E['STARS'] - net_cont_dict_E['STARS']).subspec(w0, w0 + d_wave)
    norm = norm_by_region['STARS']
    (section_mean_stars / norm - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color=color, alpha=0.2)
    ax.text(w0 + d_wave, -0.1, " STARS", color=color, alpha=0.2, ha='left', va='center', fontweight='bold')

    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    y0 = -0.1
    ax.plot([w1, w2], [y0, y0], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, y0, f"  {v_sys:.0f} km/s", ha='left', va='center')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[-0.3, shift + 0.35],
        ylabel="", xlabel="")
fig.suptitle(
    "ESO cube, continuum-subtracted spectrum by nebular region, plus mean spectrum of gas and of stars",
    fontsize='x-large',
    y=1.001,
)
# axes[0].legend(ncol=7, fontsize="large", bbox_to_anchor=(1.0, 1.0), loc='lower right')
axes[-1].set(
    xlabel='Observed wavelength, Angstrom',
)
fig.tight_layout(rect=(0, 0, 1, 1))

# + pycharm={"name": "#%%\n"}
fig.savefig("eso-region-all-strip-csub-spectra.pdf")

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Identifying lines
#
# Two parts to this:
# 1. Labeling known lines on the graphs
# 2. Automatically finding lines form the spectra, so we can try to identify the unknown lines

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Apply median filter to the entire cube
#
# I am not sure how fast this is going to be, but I will try. Start with a smaller cube
#
#

# + pycharm={"name": "#%%\n"}
_c = cubeE.select_lambda(6350, 6750)

# + [markdown] pycharm={"name": "#%% md\n"}
# Check that

# + pycharm={"name": "#%%\n"}
nwin = 51
_c_cont = get_median_continuum_cube(_c, nwin)

# + pycharm={"name": "#%%\n"}
_spec1 = _c_cont.sum(axis=(1, 2))
_spec2 = get_median_continuum(_c.sum(axis=(1, 2)), nwin)

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 5))
_spec1.plot(label='Median filter, then sum pixels')
_spec2.plot(label='Sum pixels, then median filter', linewidth=3)
_c.sum(axis=(1, 2)).plot()
ax.legend()
# ax.set(yscale='log')
ax.set(ylim=[1.e7, 2.5e7])

# + pycharm={"name": "#%%\n"}
cube_data_path = Path.cwd().parent.parent / "big-data" / "ngc346new"

# + [markdown] pycharm={"name": "#%% md\n"}
# Now try it on the whole cube

# + pycharm={"name": "#%%\n"}
nwin = 11
cube_cont  = get_median_continuum_cube(cubeE, nwin)

# + pycharm={"name": "#%%\n"}
cube_cont.write(cube_data_path / "ngc346-median-cont-011.fits", savemask='nan')

# + pycharm={"name": "#%%\n"}
(cubeE - cube_cont).write(cube_data_path / "ngc346-csub-011.fits", savemask='nan')

# + [markdown] pycharm={"name": "#%% md\n"}
# That took about 1 min with an 11 pixel window

# + pycharm={"name": "#%%\n"}
# %timeit?

# + pycharm={"name": "#%%\n"}
# %%timeit -n 1 -r 1
nwin = 101
cube_cont = get_median_continuum_cube(cubeE, nwin)

# + pycharm={"name": "#%%\n"}
# %%timeit -n 1 -r 1
cube_cont.write(cube_data_path / "ngc346-median-cont-101.fits", savemask='nan')
(cubeE - cube_cont).write(cube_data_path / "ngc346-csub-101.fits", savemask='nan')

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Subtracting spectrum from the blank box
#
# After thinking about what the Peter spectrum look like, I have a hypothesis:
# > There is an extra unwanted additive component that has continuum and absorption lines. *I suspoect that this must be scattered moonlight, which has the solar photospheric absorption lines*
#
# The next question is what to do about it. One thing to try is to find a "blank" part of the nebula and extract that spectrum and subtract it from the rest. We need to find a region that is truly blank in all respects:
# 1. No ionized emission
# 2. No neutral emission
# 3. No stars, or at least no bright ones
# 4. No dust scattering
#
# Luckily, there is such a region in the upper middle part of the nebula, so we will try and use it.

# + pycharm={"name": "#%%\n"}
blank_region = rg.Regions.read(small_data_path / "ngc346-muse-blank-box.reg")[0]
blank_region

# + [markdown] pycharm={"name": "#%% md\n"}
# Convert to a PixelRegion with the WCS from the Peter cube

# + pycharm={"name": "#%%\n"}
blank_pixel_region_P = blank_region.to_pixel(wcs=cubeP.wcs.wcs.celestial)

# + [markdown] pycharm={"name": "#%% md\n"}
# Extract the one-dimensional spectrum from  the blank region.

# + pycharm={"name": "#%%\n"}
blank_spectrum_P = get_spectrum_from_region(
    cubeP, blank_pixel_region_P, reduction_method=np.nanmean,
)

# + [markdown] pycharm={"name": "#%% md\n"}
# Make a plot of the ratio between the blank spectrum  and the spectra from different regions. We wil use both the "source" and the "background" regions to get more variety

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex='all', sharey='all')
for label in source_regions:
    (blank_spectrum_P / spec_dict_P_neb[label]).plot(label=label, ax=axes[0])
    (blank_spectrum_P / spec_dict_P_neb[label + ' bg']).plot(label=label + ' bg', ax=axes[1])

for ax in axes:
    (blank_spectrum_P / spec_mean_P_neb).plot(label='mean neb', color='k', linewidth=2, ax=ax)
    (blank_spectrum_P / spec_mean_P_stars).plot(label='mean stars', color='k', linewidth=2, ax=ax, alpha=0.3)
    ax.set(ylim=[0, 1.5], ylabel='Fraction')
    ax.legend(ncol=4)
axes[0].set(xlabel='')
axes[0].set_title('Blank spectrum divided by region spectrum')
...;

# + [markdown] pycharm={"name": "#%% md\n"}
# Conclusions from the above figure:
# 1. The blank region continuum is roughly 30% of the mean continuum in the nebula pixels, and this is independent of wavelength
# 2. For the mean star pixels, the fraction is much lower: about 5%
# 3. Most of the source and bg regions behave very similarly to the mean nebula spectrum, with the exception of `YSO`, `FIL` and `FIL bg`
# 4. `YSO` is similar to the mean star spectrum
# 5. `FIL` and `FIL bg` have fractions close to unity because they are faint regions. For `FIL bg`, the fraction goes above 1 in the blue

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = (spec_mean_P_neb - blank_spectrum_P).subspec(w0, w0 + d_wave)
    section_mean_stars = (spec_mean_P_stars - blank_spectrum_P).subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_P_neb[label] - blank_spectrum_P
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg  #/ section_mean
        norm = np.median(ratio.data)
        ((ratio / norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean / norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
    # Mean spectrum of the stars
    norm = np.median(section_mean_stars.data)
    ((section_mean_stars / norm) - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.3)
    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    ax.plot([w1, w2], [0.7, 0.7], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, 0.7, f"  {v_sys:.0f} km/s", ha='left', va='center')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[0.6, 1.7],
        # yscale="log",
        ylabel="", xlabel="")
axes[0].set_title("Peter cube, blank-subtracted spectrum by nebular region, plus mean gas and mean stars")
axes[0].legend(ncol=7, fontsize="large", loc="lower left")
axes[-1].set(
    xlabel='Observed wavelength, Angstrom',
)
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-region-strip-spectra-sub-blank.pdf")

# + [markdown] pycharm={"name": "#%% md\n"}
# As expected, the largest effect of the subtraction is on the `FIL` sample, since it is only slightly brighter than the blank spectrum. Oddly, there is very little effect on the other samples.

# + [markdown] pycharm={"name": "#%% md\n"}
# #### BG region subtraction
#
# Now try subtracting the BG region of each source region

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = (spec_mean_P_neb - blank_spectrum_P).subspec(w0, w0 + d_wave)
    section_mean_stars = (spec_mean_P_stars - blank_spectrum_P).subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_P_neb[label] - spec_dict_P_neb[label + ' bg']
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg  #/ section_mean
        norm = np.median(ratio.data)
        ((ratio / norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean / norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
    # Mean spectrum of the stars
    norm = np.median(section_mean_stars.data)
    ((section_mean_stars / norm) - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.3)
    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    ax.plot([w1, w2], [0.7, 0.7], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, 0.7, f"  {v_sys:.0f} km/s", ha='left', va='center')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[0.6, 1.7],
        # yscale="log",
        ylabel="", xlabel="")
axes[0].set_title("Peter cube, BG-subtracted spectrum by nebular region, plus mean gas and mean stars")
axes[0].legend(ncol=7, fontsize="large", loc="lower left")
axes[-1].set(
    xlabel='Observed wavelength, Angstrom',
)
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-region-strip-spectra-sub-bg.pdf")

# + [markdown] pycharm={"name": "#%% md\n"}
# This does not work great either

# + [markdown] pycharm={"name": "#%% md\n"}
# ## Return to divide by mean nebula
#
# I now think that this shows more promise than I thought.  We will try it for pairs of source and bg.
#
# First, the `GLOB` region

# + pycharm={"name": "#%%\n"}
reg_id = "GLOB"

fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    spec = spec_dict_P_neb[reg_id].subspec(w0, w0 + d_wave) / section_mean
    spec_bg = spec_dict_P_neb[reg_id + ' bg'].subspec(w0, w0 + d_wave) / section_mean
    spec_med = np.median(spec.data)
    spec_bg_med = np.median(spec_bg.data)
    (spec - spec_med).plot(ax=ax, alpha=0.6, label=reg_id)
    (spec - spec_bg).plot(ax=ax, alpha=0.6)
    (spec_bg - spec_bg_med).plot(ax=ax, alpha=0.6, label=reg_id + ' bg')

    ax.set(
        ylim=[-0.1, 0.1],
        ylabel="", xlabel="")
axes[0].set_title("Peter Ratio: Region / Mean Nebula")
axes[0].legend()
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
reg_id = "FIL"

fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    spec = spec_dict_P_neb[reg_id].subspec(w0, w0 + d_wave) / section_mean
    spec_bg = spec_dict_P_neb[reg_id + ' bg'].subspec(w0, w0 + d_wave) / section_mean
    spec_med = np.median(spec.data)
    spec_bg_med = np.median(spec_bg.data)
    (spec - spec_med).plot(ax=ax, alpha=0.6, label=reg_id)
    (spec - spec_bg).plot(ax=ax, alpha=0.6)
    (spec_bg - spec_bg_med).plot(ax=ax, alpha=0.6, label=reg_id + ' bg')

    ax.set(
        ylim=[-0.1, 0.1],
        ylabel="", xlabel="")
axes[0].set_title("Peter Ratio: Region / Mean Nebula")
axes[0].legend()
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
reg_id = "YSO"

fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    spec = spec_dict_P_neb[reg_id].subspec(w0, w0 + d_wave) / section_mean
    spec_bg = spec_dict_P_neb[reg_id + ' bg'].subspec(w0, w0 + d_wave) / section_mean
    spec_med = np.median(spec.data)
    spec_bg_med = np.median(spec_bg.data)
    (spec - spec_med).plot(ax=ax, alpha=0.6, label=reg_id)
    (spec - spec_bg - (spec_med - spec_bg_med)).plot(ax=ax, alpha=0.6)
    (spec_bg - spec_bg_med).plot(ax=ax, alpha=0.6, label=reg_id + ' bg')

    ax.set(
        ylim=[-1.5, 1.5],
        ylabel="", xlabel="")
axes[0].set_title("Peter Ratio: Region / Mean Nebula")
axes[0].legend()
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
reg_id = "NEUT-B"

fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    spec = spec_dict_extra_P_neb[reg_id].subspec(w0, w0 + d_wave) / section_mean
    spec_bg = spec_dict_extra_P_neb[reg_id + ' bg'].subspec(w0, w0 + d_wave) / section_mean
    spec_med = np.median(spec.data)
    spec_bg_med = np.median(spec_bg.data)
    (spec - spec_med).plot(ax=ax, alpha=0.6, label=reg_id)
    (spec - spec_bg - (spec_med - spec_bg_med)).plot(ax=ax, alpha=0.6)
    (spec_bg - spec_bg_med).plot(ax=ax, alpha=0.6, label=reg_id + ' bg')

    ax.set(
        ylim=[-0.02, 0.02],
        ylabel="", xlabel="")
axes[0].set_title("Peter Ratio: Region / Mean Nebula")
axes[0].legend()
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    shift = 0.0
    for reg_id in source_regions:
        spec = spec_dict_P_neb[reg_id].subspec(w0, w0 + d_wave) / section_mean
        spec_bg = spec_dict_P_neb[reg_id + ' bg'].subspec(w0, w0 + d_wave) / section_mean
        spec_med = np.median(spec.data)
        spec_bg_med = np.median(spec_bg.data)
        spec_net = (spec - spec_med) - (spec_bg - spec_bg_med)
        norm = 5 if reg_id == 'YSO' else 0.1 if reg_id == 'FIL' else 0.5
        (shift + spec_net / norm).plot(ax=ax, alpha=0.6, linewidth=3, label=reg_id)
        shift += 0.1
    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    ax.plot([w1, w2], [0.5, 0.5], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, 0.5, f"  {v_sys:.0f} km/s", ha='left', va='center', fontsize='small')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[-0.1, 0.7],
        ylabel="", xlabel="")
axes[0].set_title("Peter Cube: Median-subtracted [(Source – BG) / Mean Nebula] for each region")
axes[0].legend(ncol=5, fontsize='xx-small', loc='upper left')
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-region-strip-ratio-diff-spectra.pdf")

# + pycharm={"name": "#%%\n"}
n_sections = 12
start_waves, stop_waves, d_wave = split_spectrum(cubeP, n_sections)

fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    shift = 0.0
    for reg_id in source_regions_extra:
        spec = spec_dict_extra_P_neb[reg_id].subspec(w0, w0 + d_wave) / section_mean
        spec_bg = spec_dict_extra_P_neb[reg_id + ' bg'].subspec(w0, w0 + d_wave) / section_mean
        spec_med = np.median(spec.data)
        spec_bg_med = np.median(spec_bg.data)
        spec_net = (spec - spec_med) - (spec_bg - spec_bg_med)
        norm = 0.5
        (shift + spec_net / norm).plot(ax=ax, alpha=0.6, linewidth=3, label=reg_id)
        shift += 0.05
    # Indicator for rest-to-observed wavelength transformation
    w1 = w0 + 0.9 * d_wave
    w2 = w1 * (1.0 + v_sys / 3e5)
    ax.plot([w1, w2], [0.5, 0.5], color='m', lw=5, solid_capstyle='butt')
    ax.text(w2, 0.5, f"  {v_sys:.0f} km/s", ha='left', va='center', fontsize='small')

    # Put a finer wavelength grid with lines every 10 AA
    ax.minorticks_on()
    ax.grid(which='minor', axis='x', color='y', lw=0.5, alpha=0.5)
    ax.set(
        ylim=[-0.1, 0.7],
        ylabel="", xlabel="")
axes[0].set_title("Peter Cube: Median-subtracted [(Source – BG) / Mean Nebula] for each region")
axes[0].legend(ncol=5, fontsize='xx-small', loc='upper left')
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-region-extras-strip-ratio-diff-spectra.pdf")

# + pycharm={"name": "#%%\n"}


# + pycharm={"name": "#%%\n"}


# + [markdown] pycharm={"name": "#%% md\n"}
# ## Older stuff
#
# This is things that no longer make much sense. Should probably delete. Leaving around for now in case anything can be repurposed maybe later

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Divide by mean nebula

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_P_neb[label]  #- spec_dict_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg / section_mean
        norm = np.median(ratio.data)
        lw = 1 if label == 'YSO' else 2
        ((ratio / norm) + shift).plot(ax=ax, label=label, linewidth=lw, alpha=0.6)
        shift += 0.03
    norm = np.median(section_mean.data)
    # ((section_mean/norm) + shift).plot(ax=ax, label="mean cube", linewidth=3, color="k")
    ax.set(
        ylim=[0.95, 1.2],
        # yscale="log",
        ylabel="", xlabel="")
axes[0].set_title("Peter Ratio: Region / Mean Nebula")
axes[0].legend(ncol=6, fontsize="small")
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-region-strip-ratio-spectra.pdf")

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12, 2 * n_sections))
for w0, ax in zip(start_waves, axes):
    section_mean = spec_mean_E_neb.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_E_neb[label]  #- spec_dict_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg / section_mean
        norm = np.median(ratio.data)
        lw = 1 if label == 'YSO' else 2
        ((ratio / norm) + shift).plot(ax=ax, label=label, linewidth=lw, alpha=0.6)
        shift += 0.03
    norm = np.median(section_mean.data)
    # ((section_mean/norm) + shift).plot(ax=ax, label="mean cube", linewidth=3, color="k")
    ax.set(
        ylim=[0.95, 1.2],
        # yscale="log",
        ylabel="", xlabel="")
axes[0].set_title("ESO Ratio: Region / Mean Nebula")
axes[0].legend(ncol=6, fontsize="small")
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("eso-region-strip-ratio-spectra.pdf")

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Even older, or at least less relevant

# + pycharm={"name": "#%%\n"}
spec_dict_P = {label: get_spectrum_from_region(cubeP, reg) for label, reg in region_dict.items()}
spec_dict_E = {label: get_spectrum_from_region(cubeE, reg) for label, reg in region_dict.items()}

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 8))
for label, spec in spec_dict_P.items():
    if label.endswith(" bg"):
        lw = 1
    else:
        lw = 2
    (spec / specP).plot(label=label, linewidth=lw)
ax.legend(ncol=5)
ax.set(ylim=[0, 0.02])
...;

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(12, 8))
for label, spec in spec_dict_E.items():
    if label.endswith(" bg"):
        lw = 1
    else:
        lw = 2
    (spec / specE).plot(label=label, linewidth=lw)
ax.legend(ncol=5)
ax.set(ylim=[0, 0.03])
...;

# + pycharm={"name": "#%%\n"}
cubeP.info()

# + pycharm={"name": "#%%\n"}
type(cubeP.wave)

# + pycharm={"name": "#%%\n"}
cubeP.unit, specP.unit

# + pycharm={"name": "#%%\n"}

