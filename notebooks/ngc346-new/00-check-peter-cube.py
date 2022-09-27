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
from  astropy.table import Table
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
# It is easier to fix that at the FITS header level before the `Cube` is created

# + pycharm={"name": "#%%\n"}
hdulist = fits.open(data_path / "ADP.2017-10-16T11_04_19.247.fits")
hdulist.info()

# + pycharm={"name": "#%%\n"}
WCS(hdulist['DATA'].header)

# + [markdown] pycharm={"name": "#%% md\n"}
# Read my hand-tweaked WCS that I exported from DS9

# + pycharm={"name": "#%%\n"}
eso_gaia_hdr = fits.Header(wcsfile.read(small_data_path / "ngc346-muse-deep-3d-gaia.wcs"))
eso_gaia_hdr

# + [markdown] pycharm={"name": "#%% md\n"}
# Update the header info in both the `DATA` and `STAT` HDUs

# + pycharm={"name": "#%%\n"}
hdulist[1].header.update(eso_gaia_hdr)
hdulist[2].header.update(eso_gaia_hdr)

# + pycharm={"name": "#%%\n"}
cubeE = Cube(hdulist=hdulist)

# + [markdown] pycharm={"name": "#%% md\n"}
# Find offset between the two frames

# + pycharm={"name": "#%%\n"}
(cubeE.wcs.wcs
 .array_index_to_world(0, 0)
 .separation(
    cubeP.wcs.wcs
    .array_index_to_world(0, 0)
).arcsec)

# + [markdown] pycharm={"name": "#%% md\n"}
# That is about 0.3 pixels,  which should be good enough. On the other hand, maybe it would be better to just copy the Peter WCS to the ESO cube. Another day maybe.

# + pycharm={"name": "#%%\n"}
cubeE.wave

# + pycharm={"name": "#%%\n"}
cubeP.wave

# + pycharm={"name": "#%%\n"}
300000 * (4599.97314453125 - 4599.94482421875) / 4599.97314453125

# + pycharm={"name": "#%%\n"}


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
# Choose 12 sections, which means that each section spans about 400 Angstrom.

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
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for section_P, section_E, ax in zip(spec_sections_P, spec_sections_E,  axes):
    section_P.plot(ax=ax)
    section_E.plot(ax=ax)
    ax.set(yscale="log", ylabel="", xlabel="")
fig.tight_layout()

# + [markdown] pycharm={"name": "#%% md\n"}
# So, the strongest discrepancy is the [O I] 5577 line, which will be extremely weak  from  the nebula, but is strong from the night sky. This is almost completely removed in the ESO cube but is still strong in the Peter cube.

# + [markdown] pycharm={"name": "#%% md\n"}
# The Peter cube is also higher at the position of real nebula lines, such as H beta and [O III].  In these cases, I think it is not true night sky emission, but is rather the outer parts of the nebula that are being picked up. These are over-subtracted in the ESO cube, resulting in the lines being apparently seen in absorption in fainter spaxels of the cube,  The Peter cube does not have this problem.

# + [markdown] pycharm={"name": "#%% md\n"}
# Next, we look at the ratio of the two cubes, zooming in on the range +/-20% so we can see the weak night sky lines and atmospheric absorption.

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for section_P, section_E, ax in zip(spec_sections_P, spec_sections_E,  axes):
    (section_P / section_E ).plot(ax=ax)
    ax.set(ylim=[0.8, 1.2], ylabel="", xlabel="")
axes[0].set_title("Ratio: Peter / ESO")
fig.tight_layout()

# + [markdown] pycharm={"name": "#%% md\n"}
# We see a lot of weak absorption features with EW of less than an Angstrom. Also, lots of emission lines. In most sections, it is possible to distinguish the two, but at longer wavelengths this is more difficult.
#
# We could look at the model sky spectra to try and help with this.  We probably need to be able to separate the sky emission from the absorption, since the first is additive, while the second is multiplicative. **Actually, this might be a way of doing it:** in principle the absorption lines should be spatially constant in a map of `cube/avspec`, whereas the emission lines should be constant in a map of `(cube - avspec)`

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Fitting the continuum
#
# I had an idea that we can use sigma clipping in the fit so that we do not have to bother with masking out the emission lines first.  I will try this on the entire cube first, since I  am not sure how fast it is going to be to apply it pixel by pixel.
#
# - [ ] TODO do the fits

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
ha_map_P.shape,   ha_map_E.shape

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
        reduction_method: callable=np.sum,
        debug: bool=False,
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
mask_yso_P

# + [markdown] pycharm={"name": "#%% md\n"}
# #### Calculate mean spectra for stars and for diffuse nebula
#
# Remember that the cube mask is `True` for invalid pixels

# + pycharm={"name": "#%%\n"}
cube_P_neb = cubeP.copy()
cube_P_neb.mask = cube_P_neb.mask | mask_stars_P | mask_yso_P
spec_mean_P_neb = cube_P_neb.mean(axis=(1, 2))

# + pycharm={"name": "#%%\n"}
cube_P_stars = cubeP.copy()
cube_P_stars.mask = cube_P_stars.mask | (~mask_stars_P)
spec_mean_P_stars = cube_P_stars.mean(axis=(1, 2))

# + pycharm={"name": "#%%\n"}
cube_E_neb = cubeE.copy()
cube_E_neb.mask = cube_E_neb.mask | mask_stars_E
spec_mean_E_neb = cube_E_neb.mean(axis=(1, 2))

# + pycharm={"name": "#%%\n"}
cube_E_stars = cubeE.copy()
cube_E_stars.mask = cube_E_stars.mask | (~mask_stars_E)
spec_mean_E_stars = cube_E_stars.mean(axis=(1, 2))

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

# + pycharm={"name": "#%%\n"}


# + pycharm={"name": "#%%\n"}


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
pec_dict_E_neb = {
    label: get_spectrum_from_region(cube_E_neb, reg, reduction_method=np.nanmean)
    for label, reg in region_dict.items()
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
    'cmr.ember', 10, cmap_range=(0.25, 0.95), return_fmt='hex')

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
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    section_mean_stars = spec_mean_P_stars.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_P_neb[label] #- spec_dict_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg #/ section_mean
        norm = np.median(ratio.data)
        ((ratio/norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean/norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
    # Mean spectrum of the stars
    norm = np.median(section_mean_stars.data)
    ((section_mean_stars/norm) - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.3)
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
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    section_mean_stars = spec_mean_P_stars.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions_extra:
        spec = spec_dict_extra_P_neb[label] #- spec_dict_extra_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg #/ section_mean
        norm = np.median(ratio.data)
        ((ratio/norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.05
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean/norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
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

fig, axes = plt.subplots(n_row, n_col, figsize=(12,  2 * n_row), sharey="all")
for w0, w1, ax in zip(start_waves,  stop_waves, axes.flat):
    section_mean = spec_mean_P_neb.subspec(w0, w1)
    section_mean_stars = spec_mean_P_stars.subspec(w0, w1)
    spec_list = []
    for label in all_source_regions:
        spec = all_spec_dict_P_neb[label] #- all_spec_dict_P_neb[label + " bg"]
        sub_spec = spec.subspec(w0, w1)
        norm = np.median(sub_spec.data)
        color = color_dict[label.split('-')[0]]
        spec_list.append(sub_spec.data.data / norm)
        (sub_spec / norm).plot(ax=ax, label=label, linewidth=0.5, alpha=0.6, color=color)
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean/norm)).plot(ax=ax, label="mean nebula", linewidth=0.5, color="k")
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
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = spec_mean_E_neb.subspec(w0, w0 + d_wave)
    section_mean_stars = spec_mean_E_stars.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_E_neb[label] #- spec_dict_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg #/ section_mean
        norm = np.median(ratio.data)
        ((ratio/norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean/norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
    # Mean spectrum of the stars
    norm = np.median(section_mean_stars.data)
    ((section_mean_stars/norm) - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.3)
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
blank_region = rg.Regions.read(small_data_path / "ngc346-muse-blank-box.reg" )[0]
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
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = (spec_mean_P_neb - blank_spectrum_P).subspec(w0, w0 + d_wave)
    section_mean_stars = (spec_mean_P_stars - blank_spectrum_P).subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_P_neb[label] - blank_spectrum_P
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg #/ section_mean
        norm = np.median(ratio.data)
        ((ratio/norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean/norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
    # Mean spectrum of the stars
    norm = np.median(section_mean_stars.data)
    ((section_mean_stars /norm) - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.3)
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
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = (spec_mean_P_neb - blank_spectrum_P).subspec(w0, w0 + d_wave)
    section_mean_stars = (spec_mean_P_stars - blank_spectrum_P).subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_P_neb[label] - spec_dict_P_neb[label + ' bg']
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg #/ section_mean
        norm = np.median(ratio.data)
        ((ratio/norm) + shift).plot(ax=ax, label=label, linewidth=2, alpha=0.6)
        shift += 0.1
    # Mean spectrum of the nebula
    norm = np.median(section_mean.data)
    ((section_mean/norm) + shift).plot(ax=ax, label="mean nebula", linewidth=3, color="k")
    # Mean spectrum of the stars
    norm = np.median(section_mean_stars.data)
    ((section_mean_stars /norm) - 0.1).plot(ax=ax, label="mean stars", linewidth=3, color="k", alpha=0.3)
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

fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
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

fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
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

fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
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

fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
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
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    shift = 0.0
    for reg_id in source_regions:
        spec = spec_dict_P_neb[reg_id].subspec(w0, w0 + d_wave) / section_mean
        spec_bg = spec_dict_P_neb[reg_id + ' bg'].subspec(w0, w0 + d_wave) / section_mean
        spec_med = np.median(spec.data)
        spec_bg_med = np.median(spec_bg.data)
        spec_net = (spec - spec_med) - (spec_bg - spec_bg_med)
        norm = 5 if reg_id == 'YSO' else 0.1 if reg_id == 'FIL' else 0.5
        (shift + spec_net/norm).plot(ax=ax, alpha=0.6, linewidth=3, label=reg_id)
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

fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    shift = 0.0
    for reg_id in source_regions_extra:
        spec = spec_dict_extra_P_neb[reg_id].subspec(w0, w0 + d_wave) / section_mean
        spec_bg = spec_dict_extra_P_neb[reg_id + ' bg'].subspec(w0, w0 + d_wave) / section_mean
        spec_med = np.median(spec.data)
        spec_bg_med = np.median(spec_bg.data)
        spec_net = (spec - spec_med) - (spec_bg - spec_bg_med)
        norm = 0.5
        (shift + spec_net/norm).plot(ax=ax, alpha=0.6, linewidth=3, label=reg_id)
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
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_P_neb[label] #- spec_dict_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg / section_mean
        norm = np.median(ratio.data)
        lw = 1 if label == 'YSO' else 2
        ((ratio/norm) + shift).plot(ax=ax, label=label, linewidth=lw, alpha=0.6)
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
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = spec_mean_E_neb.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_E_neb[label] #- spec_dict_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg / section_mean
        norm = np.median(ratio.data)
        lw = 1 if label == 'YSO' else 2
        ((ratio/norm) + shift).plot(ax=ax, label=label, linewidth=lw, alpha=0.6)
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

