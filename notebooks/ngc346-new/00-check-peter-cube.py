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

# # Compare Peter Zeidler's NGC 346 cube with ESO data archive cube of same observations

# + pycharm={"name": "#%%\n"}
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from  astropy.table import Table
import seaborn as sns
from mpdaf.obj import Cube, Image, Spectrum
import regions as rg
sns.set_context("talk")
# -

# ## Spectrum of the whole cube

datapath = Path.home() / "Work/Muse-Hii-Data/SMC-NGC-346"
cubeP = Cube(str(datapath / "PeterZeidler" / "DATACUBE_FINAL_fwhm_cor.fits"))

# Also load the standard ESO cube for comparison

cubeE = Cube(str(datapath / "ADP.2017-10-16T11_04_19.247.fits"))

# Extract full-cube spectra. This next cell takes several seconds to execute.

specP = cubeP.sum(axis=(1, 2))
specE = cubeE.sum(axis=(1, 2))

# ### Global comparison of full spectrum
#
# First look at the full wav range on a log intensity scale

fig, ax = plt.subplots(figsize=(12, 4))
specP.plot(ax=ax, label="Peter")
specE.plot(ax=ax, label="ESO")
ax.legend()
ax.set(
    yscale="log",
)
...;

# So Peter cube is consistently brighter, presumably due to better photometry. Also has more sky features visible in the spectrum.
#
# Next, look ast the ratio between the two on a linear scale:

fig, ax = plt.subplots(figsize=(12, 4))
(specP / specE).plot(ax=ax, label="Ratio")
ax.legend()
avratio = np.nanmedian(specP.data / specE.data)
ax.axhline(avratio, color='r', zorder=-100)
ax.set(
    ylim=[0.0, None],
)
...;

avratio

# ### Split into multiple wavelength sections
#
# Choose 12 sections, which means that each section spans about 400 Angstrom.

wave_min, wave_max = specP.wave.get_range()
n_sections = 12
d_wave = (wave_max - wave_min) / n_sections
start_waves = wave_min + np.arange(n_sections) * d_wave
start_waves

# Scale the ESO cube by the median ratio, so we can compare in more detail.

spec_sections_P = [specP.subspec(w0, w0 + d_wave) for w0 in start_waves]
spec_sections_E = [avratio * specE.subspec(w0, w0 + d_wave) for w0 in start_waves]

# First compare on a logarithmic flux scale.

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for section_P, section_E, ax in zip(spec_sections_P, spec_sections_E,  axes):
    section_P.plot(ax=ax)
    section_E.plot(ax=ax)
    ax.set(yscale="log", ylabel="", xlabel="")
fig.tight_layout()
# -

# So, the strongest discrepancy is the [O I] 5577 line, which will be extremely weak  from  the nebula, but is strong from the night sky. This is almost completely removed in the ESO cube but is still strong in the Peter cube.

# + [markdown] pycharm={"name": "#%% md\n"}
# The Peter cube is also higher at the position of real nebula lines, such as H beta and [O III].  In these cases, I think it is not true night sky emission, but is rather the outer parts of the nebula that are being picked up. These are over-subtracted in the ESO cube, resulting in the lines being apparently seen in absorption in fainter spaxels of the cube,  The Peter cube does not have this problem.
# -

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
# -

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
# -

# ## Look at particular regions
# To start with I will use the bowshock ones that I made for Jesus to use for the Spitzer spectra. Later, I should do some that concentrate more on the filaments.

# + pycharm={"name": "#%%\n"}
small_data_path = Path.cwd().parent.parent / "data"
region_file = "ngc346-jesus-muse-pixels.reg"
pixel_regions = rg.Regions.read(small_data_path / region_file)

# + [markdown] pycharm={"name": "#%% md\n"}
# I am using pixel coordinates rather than  celestial since the Peter cube has been aligned to MUSE, whereas the ESO cube has not.  There is a slight difference in cube shape, so we should check if that matters.

# + pycharm={"name": "#%%\n"}
region_dict = {reg.meta["label"]: reg for reg in pixel_regions}
Table(rows=region_dict.items())

# + [markdown] pycharm={"name": "#%% md\n"}
# Calculate the fraction of the full spectrum flux that is in the Ha line.  We will use this as a criterion for whether a given pixel is dominated by stars or nebular emission. For example, using a threshold of 0.04.

# + pycharm={"name": "#%%\n"}
ha_frac_P = ha_map_P / white_map_P


# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(8, 8))
ha_map_P.plot(ax=ax, zscale=True, cmap="gray_r", title="Positions of sample regions")
ax.contourf(ha_frac_P.data, levels=[0.00, 0.01, 0.02, 0.03, 0.04], cmap="Greens")
for label, reg in region_dict.items():
    is_bg = label.endswith("bg")
    linestyle = "dotted" if is_bg else "solid"
    fontweight = "normal" if is_bg else "bold"
    reg.plot(ax=ax, facecolor='none', edgecolor='red',
             lw=2,

             linestyle=linestyle,
             )
    ax.text(reg.center.x, reg.center.y, label,
            ha="center", va="center", fontsize="xx-small", fontweight=fontweight, color="w")
...;


# + [markdown] pycharm={"name": "#%% md\n"}
# The green blobs are have a low Ha fraction and are dominated by the stellar continuum. We could maybe mask them out.

# + pycharm={"name": "#%%\n"}
def get_spectrum_from_region(
        cube: Cube,
        region: rg.PixelRegion,
        reduction_method: callable=np.sum,
) -> Spectrum:
    region_mask = region.to_mask()
    nv, ny, nx = cube.shape
    # Slices into 2D arrays
    slices_large, slices_small = region_mask.get_overlap_slices((ny, nx))
    slices_cube = (slice(None, None),) + slices_large
    image_mask_large = region_mask.to_image((ny, nx))
    image_mask_small = image_mask_large[slices_large]
    cube_cutout = cube.data[slices_cube]
    spec = reduction_method(cube_cutout * image_mask_small[None, :, :], axis=(1, 2))
    return Spectrum(wave=cube.wave, data=spec, unit=cube.unit)



# + pycharm={"name": "#%%\n"}
type(np.sum)

# + pycharm={"name": "#%%\n"}


# + pycharm={"name": "#%%\n"}
mask_stars_P = ha_frac_P.data < 0.04

# + pycharm={"name": "#%%\n"}
cube_P_neb = cubeP.copy()
cube_P_neb.mask = cube_P_neb.mask | mask_stars_P
spec_mean_P_neb = cube_P_neb.mean(axis=(1, 2))

# + pycharm={"name": "#%%\n"}
spec_dict_P_neb = {
    label: get_spectrum_from_region(cube_P_neb, reg, reduction_method=np.nanmean)
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
source_regions = set(s.split()[0] for s in region_dict.keys())
source_regions

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(n_sections, 1, figsize=(12,  2 * n_sections))
for w0, ax in zip(start_waves,  axes):
    section_mean = spec_mean_P_neb.subspec(w0, w0 + d_wave)
    shift = 0.0
    for label in source_regions:
        spec = spec_dict_P_neb[label] #- spec_dict_P_neb[label + " bg"]
        section_reg = spec.subspec(w0, w0 + d_wave)
        ratio = section_reg #/ section_mean
        norm = np.median(ratio.data)
        ((ratio/norm) + shift).plot(ax=ax, label=label, linewidth=2)
        shift += 0.1
    norm = np.median(section_mean.data)
    ((section_mean/norm) + shift).plot(ax=ax, label="mean cube", linewidth=3, color="k")
    ax.set(
        ylim=[0.5, 1.7],
        # yscale="log",
        ylabel="", xlabel="")
axes[0].set_title("Ratio: Region / Full Cube")
axes[0].legend(ncol=6, fontsize="x-small")
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-region-strip-spectra.pdf")

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
        ((ratio/norm) + shift).plot(ax=ax, label=label, linewidth=2)
        shift += 0.1
    norm = np.median(section_mean.data)
    # ((section_mean/norm) + shift).plot(ax=ax, label="mean cube", linewidth=3, color="k")
    ax.set(
        ylim=[0.5, 1.7],
        # yscale="log",
        ylabel="", xlabel="")
axes[0].set_title("Ratio: Region / Full Cube")
axes[0].legend(ncol=6, fontsize="x-small")
fig.tight_layout()

# + pycharm={"name": "#%%\n"}
fig.savefig("peter-region-strip-ratio-spectra.pdf")


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
# -

# ### Divide  by the mean cube


