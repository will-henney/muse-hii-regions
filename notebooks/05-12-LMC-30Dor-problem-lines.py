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

# # Revisiting problematic lines
#
# Some of the line maps have discontinuities between the ABCD fields. In this notebook, we are going to work out what is wrong and fix it.

# ## [O I] 6300 line 
#
# This has a sky component, which is probably causing the problem. However, we ought to be able to remove it since the LMC is redshifted by over 200 km/s.
#
# Also, we can use the ratio of 6363/6300 in the faint parts to guide us. 

# +
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import cmasher as cm
from mpdaf.obj import Cube
import regions
import sys
import warnings
from whispy.linetools import EmissionLine, SpectralRange, VelocityScale
from whispy import moments
from whispy import extract
from whispy import sky

import astropy.units as u
import astropy.constants as const

sns.set_context("talk")
sns.set_color_codes()
# -

# Start off by using the subcubes, as we did in the 05-04 notebook

labels = "ABCD"
csub = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-62-71-contsub.fits") 
    for label in labels
}
corig = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-62-71.fits") 
    for label in labels
}

# Use the new classes in linetools

SpectralRange(6300.30).wavlim

# So, we should no longer need the following cell:

#C_KMS = const.c.to(u.km/u.s).value
#
#def wavlimits(wav0, vlim=[100.0, 300.0]):
#    wav1 = wav0 * (1.0 + vlim[0] / C_KMS)
#    wav2 = wav0 * (1.0 + vlim[1] / C_KMS)
#    return wav1, wav2
      

oi_window_median = {label: cube.select_lambda(6250, 6400).median(axis=(1, 2)) for label, cube in csub.items()}

# +
fig, ax = plt.subplots(figsize=(12, 6))
for label, spec in oi_window_median.items():
    spec.plot(label=label, linewidth=2)

#wav1, wav2 = wavlimits(6300.30)
for wav0 in 6300.30, 6363.78, 6312.06, 6347.11:
    ax.axvspan(*SpectralRange(wav0, vlim=[100, 300]).wavlim, alpha=0.2)
    ax.axvspan(*SpectralRange(wav0, vlim=[-50, 50]).wavlim, alpha=0.1, color="k")
ax.legend()
ax.set_title("Median spectra around [O I] 6300, 6363 lines")
...;
# -

# So, a few things:
#
# 1. This is the median of the BG-subtracted cube of each of the ABCD fields. I originally tried with the mean, but that was weird for the BG-subtracted ones (presumably there are a few pixels that are problematic, mainly where there is a bright star). 
# 2. It looks like in fields B and D the sky has been oversubtracted, whereas in A and C it has been undersubtracted. 
# 3. I am showing the window for velocity ranges of [-50, 50] (gray) and [100, 300] (blue) for the O I, S III and Si II lines. 
# 4. The latter is the default window that I have been using for line extraction, and it looks like it is missing a lot of the flux on the red side of the line profile.
#
#

# So the first solution to try is just to use a better extraction window. Maybe [150, 400]

# +
fig, ax = plt.subplots(figsize=(12, 6))
for label, spec in oi_window_median.items():
    spec.plot(label=label, linewidth=2)

newvlim = [125, 425]
for wav0 in 6300.30,:
    ax.axvspan(*SpectralRange(wav0, vlim=newvlim).wavlim, alpha=0.2)
    ax.axvspan(*SpectralRange(wav0, vlim=[-50, 50]).wavlim, alpha=0.1, color="k")
ax.legend()
ax.set_title(f"Improved extraction window for [O I] 6300: {newvlim}")
ax.set(
    xlim=[6290, 6330],
)
...;
# -

# ### Calculate moments with new extraction window

oimoms = {
    label:
    moments.find_moments(
        cube.select_lambda(*SpectralRange(wav0, vlim=newvlim).wavlim)
    )
    for label, cube in csub.items()
}

scale = 500
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oimoms["D"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 0])
axes[0, 0].set_title("D")
oimoms["C"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 1])
axes[0, 1].set_title("C")
oimoms["B"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 0])
axes[1, 0].set_title("B")
im = oimoms["A"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 1])
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# Yay, all positive now!  (With the exception of a few stars, which presumably have the absorption line).
#
# We will look at the velocities too to see if they are improved any:

kws = dict(vmin=6304, vmax=6308, cmap=cm.fusion_r)
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oimoms["D"][1].plot(ax=axes[0, 0], **kws)
axes[0, 0].set_title("D")
oimoms["C"][1].plot(ax=axes[0, 1], **kws)
axes[0, 1].set_title("C")
oimoms["B"][1].plot(ax=axes[1, 0], **kws)
axes[1, 0].set_title("B")
im = oimoms["A"][1].plot(ax=axes[1, 1], **kws)
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# Yes, that looks a lot better too. Note that the the units are still wavelength, not velocity.

medlam = np.mean(oimoms["A"][1].data)
medlam

# We can convert to velocity to compare with our new window.

const.c.to(u.km/u.s) * (medlam - wav0) / wav0, newvlim

# So that is pretty well centered.  There does still seem to be a little bit of an effect of the sky still. For instance, there looks like there is a discontinuity between C and D

np.mean(oimoms["A"][2].data)

kws = dict(vmin=0, vmax=2, cmap=cm.horizon)
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oimoms["D"][2].plot(ax=axes[0, 0], **kws)
axes[0, 0].set_title("D")
oimoms["C"][2].plot(ax=axes[0, 1], **kws)
axes[0, 1].set_title("C")
oimoms["B"][2].plot(ax=axes[1, 0], **kws)
axes[1, 0].set_title("B")
im = oimoms["A"][2].plot(ax=axes[1, 1], **kws)
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# ### Repeat for a window that isolates the sky

skyvlim = [-100, 100]
oiskymoms = {
    label:
    moments.find_moments(
        cube.select_lambda(*SpectralRange(wav0, vlim=skyvlim).wavlim)
    )
    for label, cube in csub.items()
}

scale = 2000
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oiskymoms["D"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 0])
axes[0, 0].set_title("D")
oiskymoms["C"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 1])
axes[0, 1].set_title("C")
oiskymoms["B"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 0])
axes[1, 0].set_title("B")
im = oiskymoms["A"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 1])
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# This all looks very uniform in each field.  One idea would be to fit a Gaussian to the median of each field and subtract it.

kws = dict(vmin=6298, vmax=6302, cmap=cm.fusion_r)
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oiskymoms["D"][1].plot(ax=axes[0, 0], **kws)
axes[0, 0].set_title("D")
oiskymoms["C"][1].plot(ax=axes[0, 1], **kws)
axes[0, 1].set_title("C")
oiskymoms["B"][1].plot(ax=axes[1, 0], **kws)
axes[1, 0].set_title("B")
im = oiskymoms["A"][1].plot(ax=axes[1, 1], **kws)
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# ### Fit the Sky component with a Gaussian

# New version that uses the `sky` module

em = EmissionLine("[O I] 6300", wav0)
em

# The default `sky_vlim` was `-50, 50`, which was not enough pixels to do the fit, so use `-100, 100` instead

oi_sky = {
    label: sky.Sky(cube, em, sky_vlim=[-100, 100], full_vlim=[-400, 600])
    for label, cube in csub.items()
}
oi_sky

# +
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=False)

for ax, [label, skyline] in zip(axes.flat, oi_sky.items()):
    finewavs = np.linspace(
        skyline.wavs.min(), 
        skyline.wavs.max(), 
        200,
    )
    ax.axvspan(*skyline.sky_range.wavlim, alpha=0.1)
    skyline.spec.plot(label=label, linewidth=2, ax=ax)
    ax.plot(finewavs, skyline.gauss(finewavs))
    residuals = skyline.spec.data - skyline.gauss(skyline.wavs)
    ax.scatter(skyline.wavs, residuals)
    ax.legend()
    ax.axvline(skyline.em.wav0, color="r")
    ax.axhline(0, color="r")
fig.tight_layout()
...;
# -

# That seems to have worked OK, very similar to the last time.

# The fits look like they worked very well. The residuals are shown by the blue dots.  We get positive residual on the right edge because of the nebular component starting to rise, which is totally expected. 
#
# We can look at the fitted models:

[_.gauss for _ in oi_sky.values()]

# So it looks like A was observed on a different date from the others, since the mean wavelength is significantly bluer. If we could be bothered, it would be good to calculate the heliocentric correction for each observing date to make sure that these values make sense. 
#
# Now compare with the extraction window again:

oi_skysub = {label: sky.remove_sky(skyline) for label, skyline in oi_sky.items()}

newvlim

# +
fig, ax = plt.subplots(figsize=(12, 6))
for label, skysub in oi_skysub.items():
    newspec = skysub.median(axis=(1, 2))
    norm = np.max(newspec.subspec(*SpectralRange(wav0, vlim=newvlim).wavlim).data)
    newspec /= norm
    newspec.plot(label=label, linewidth=2)

newvlim = [50, 500]
for wav0 in 6300.30,:
    ax.axvspan(*SpectralRange(wav0, vlim=newvlim).wavlim, alpha=0.1)
    ax.axvspan(*SpectralRange(wav0, vlim=[-50, 50]).wavlim, alpha=0.1, color="k")
ax.legend()
ax.set_title(f"Sky-subtracted profiles for [O I] 6300: vel window = {newvlim}")
ax.set(
    xlim=[6290, 6330],
    ylim=[-0.1, 1.1],
)
...;
# -

# ### Repeat the moment calculation with the sky-subtracted cube

oimomss = {
    label:
    moments.find_moments(cube.select_lambda(*SpectralRange(wav0, vlim=newvlim).wavlim))
    for label, cube in oi_skysub.items()
}

scale = 1000
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oimomss["D"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 0])
axes[0, 0].set_title("D")
oimomss["C"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 1])
axes[0, 1].set_title("C")
oimomss["B"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 0])
axes[1, 0].set_title("B")
im = oimomss["A"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 1])
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

kws = dict(vmin=6304, vmax=6308, cmap=cm.fusion_r)
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oimomss["D"][1].plot(ax=axes[0, 0], **kws)
axes[0, 0].set_title("D")
oimomss["C"][1].plot(ax=axes[0, 1], **kws)
axes[0, 1].set_title("C")
oimomss["B"][1].plot(ax=axes[1, 0], **kws)
axes[1, 0].set_title("B")
im = oimomss["A"][1].plot(ax=axes[1, 1], **kws)
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

kws = dict(vmin=0.5, vmax=2, cmap=cm.heat)
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oimomss["D"][2].plot(ax=axes[0, 0], **kws)
axes[0, 0].set_title("D")
oimomss["C"][2].plot(ax=axes[0, 1], **kws)
axes[0, 1].set_title("C")
oimomss["B"][2].plot(ax=axes[1, 0], **kws)
axes[1, 0].set_title("B")
im = oimomss["A"][2].plot(ax=axes[1, 1], **kws)
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# These are way better.  The mean and stddev are noisy in the faint parts of course, but we are no lomger seeing the suspicious behavior that we saw before. 

# Now we just need to find a way to generalise this.
#
# But first, have a look at subtracting off the sky mean wavelength to see if it eliminates some of the pattern noise. 

(oimomss["D"][1] - oiskymoms["D"][1]).plot(vmin=3.5, vmax=7.5, cmap=cm.fusion_r)

# Does that help? Maybe a bit, but probably not enough to be worth doing.

# ### Look at the corner plots of the line moments in the different panels

for label in "ABCD":
    moments.moments_corner_plot(
        oimomss[label], 
        rebin=8,
        hist_bins=50,
        image_bins=25,
        flabel=f"{label}-oi", 
        ilabel=f"Field {label} [O I]",
        label="6300",
        restwav=wav0,
        irange=[20, 8000],
        vrange=[180, 320],
        srange=[20, 120],
    )
...;

# It is surprising how varied these are, but it is mainly because they are sampling different parts of brightness distribution: A and C are mainly bright parts and are all quite close to the V=+250 km/s and with narrow widths. B and D are mainly fainter emission, which has more variety of V and with generally higher sigma. 
#
# Quite similar behavior is seen in the Ha profiles that I looke at for the Javier project. 

# ## Fluorescent O I lines
#
# The 8446 line had a similar problem but to a much smaller degree.

csub = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-78-87-contsub.fits") 
    for label in labels
}
corig = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-78-87.fits") 
    for label in labels
}
ccont = {
    label: 
    Cube(f"../big-data/lmc-30dor-{label}-subcube-78-87-cont.fits") 
    for label in labels
}

em = EmissionLine("oi-8446", 8446.36)
oi_sky = {
    label: sky.Sky(cube, em, sky_vlim=[-150, 100], full_vlim=[-400, 600])
    for label, cube in csub.items()
}
oi_sky

# Extract the moments for the interval that would correspond to the sky:

skyvlim = [50, 100]
oiskymoms = {
    label:
    moments.find_moments(
        cube.select_lambda(*SpectralRange(em.wav0, vlim=skyvlim).wavlim)
    )
    for label, cube in csub.items()
}

# And plot the zeroth moment

scale = 200
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oiskymoms["D"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 0])
axes[0, 0].set_title("D")
oiskymoms["C"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 1])
axes[0, 1].set_title("C")
oiskymoms["B"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 0])
axes[1, 0].set_title("B")
im = oiskymoms["A"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 1])
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# _This shows that it is not the sky after all!_
#
# Instead, it is the H line Pa 18 8438
#
# This just goes to show that cannot assume that a solution for one line will work for the others!

newvlim = 100, 400

# +
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=False)

for ax, [label, skyline] in zip(axes.flat, oi_sky.items()):
    finewavs = np.linspace(
        skyline.wavs.min(), 
        skyline.wavs.max(), 
        200,
    )
    ax.axvspan(*skyline.sky_range.wavlim, alpha=0.1)
    ax.axvspan(*SpectralRange(em.wav0, vlim=newvlim).wavlim, alpha=0.1, color="g")
    skyline.spec.plot(label=label, linewidth=2, ax=ax)
    ax.plot(finewavs, skyline.gauss(finewavs))
    residuals = skyline.spec.data - skyline.gauss(skyline.wavs)
    ax.scatter(skyline.wavs, residuals)
    ax.legend()
    ax.axvline(skyline.em.wav0, color="r")
    ax.axhline(0, color="r")
fig.tight_layout()
...;
# -

[_.gauss for _ in oi_sky.values()]

# So we have successfully managed to fit this component, even though it is not actually sky. But it is of no use, since it is not spatially uniform, so subtracting the median value makes no sense. 

oi_skysub = {label: sky.remove_sky(skyline) for label, skyline in oi_sky.items()}

newvlim

# We can calculate the moments, both with and without the subtraction of the supposed sky component.  

oimoms = {
    label:
    moments.find_moments(
        cube.select_lambda(*SpectralRange(em.wav0, vlim=newvlim).wavlim)
    )
    for label, cube in csub.items()
}
oimomss = {
    label:
    moments.find_moments(
        cube.select_lambda(*SpectralRange(em.wav0, vlim=newvlim).wavlim)
    )
    for label, cube in oi_skysub.items()
}

# But there is not much use in using the skysub version. I did look at it, and it is marginally worse, as well as making absolutely no sense.

scale = 1000
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oimoms["D"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 0])
axes[0, 0].set_title("D")
oimoms["C"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 1])
axes[0, 1].set_title("C")
oimoms["B"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 0])
axes[1, 0].set_title("B")
im = oimoms["A"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 1])
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

kws = dict(vmin=em.wav0 + 5, vmax=em.wav0 + 9, cmap=cm.fusion_r)
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oimoms["D"][1].plot(ax=axes[0, 0], **kws)
axes[0, 0].set_title("D")
oimoms["C"][1].plot(ax=axes[0, 1], **kws)
axes[0, 1].set_title("C")
oimoms["B"][1].plot(ax=axes[1, 0], **kws)
axes[1, 0].set_title("B")
im = oimoms["A"][1].plot(ax=axes[1, 1], **kws)
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

kws = dict(vmin=0.5, vmax=2, cmap=cm.heat)
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
oimoms["D"][2].plot(ax=axes[0, 0], **kws)
axes[0, 0].set_title("D")
oimoms["C"][2].plot(ax=axes[0, 1], **kws)
axes[0, 1].set_title("C")
oimoms["B"][2].plot(ax=axes[1, 0], **kws)
axes[1, 0].set_title("B")
im = oimomss["A"][2].plot(ax=axes[1, 1], **kws)
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# So these are the mom0 and mom1 maps without the sky subtraaction. They are mostly fine, with the exception of the B field, which has negative parts. 
#
# This is very similar to what we saw for this line in the 05-09 notebook.  I suspect it may have something to do with the continuum fitting. 

widerange = SpectralRange(em.wav0, vlim=[-4000, 4000])
margin = 30
hislice = slice(-100, None)
loslice = slice(None, 100)
specB = csub["B"][:, -margin:, loslice].select_lambda(*widerange.wavlim).mean(axis=(1, 2))
specD = csub["D"][:, :margin,  loslice].select_lambda(*widerange.wavlim).mean(axis=(1, 2))
specA = csub["A"][:, -margin:, hislice].select_lambda(*widerange.wavlim).mean(axis=(1, 2))
specC = csub["C"][:, :margin,  hislice].select_lambda(*widerange.wavlim).mean(axis=(1, 2))

# We will compare overlapping parts of the fields: B vs D and A vs C

fig, ax = plt.subplots(figsize=(12, 4))
specB.plot(label="B")
specD.plot(label="D")
ax.axvline(VelocityScale(em.wav0).vel2wav(250))
ax.axvline(em.wav0, linestyle="dashed")
ax.axhline(0.0)
ax.legend()
...;

fig, ax = plt.subplots(figsize=(12, 4))
specA.plot(label="A")
specC.plot(label="C")
ax.axvline(VelocityScale(em.wav0).vel2wav(250))
ax.axvline(em.wav0, linestyle="dashed")
ax.axhline(0.0)
ax.legend()
...;

# So it is clear there is some strange oversubtraction of the "sky" in field B that causes lots of the lines to appear in absorption.
#
# This is particularly apparent on the outer edge.
#
# For A and C, on the other hand, the agreement is very good.

for label in "ABCD":
    moments.moments_corner_plot(
        oimoms[label], 
        rebin=8,
        hist_bins=50,
        image_bins=25,
        flabel=f"{label}-oi", 
        ilabel=f"Field {label} O I",
        label="8446",
        restwav=em.wav0,
        irange=[20, 8000],
        vrange=[180, 320],
        srange=[20, 120],
    )
...;

# The bright ends of the A, C, and D distributions are similar to 6300, but the faint ends are all different, each in its own unique way. 
#
# This is almost certainly due to the effects of contaminating lines and/or bad sky.  In the case of field B, there is no bright end and everything is a mess. 

# ## Si I 8150.57 line or [Fe I] 8151.3424 
#
# Formerly a mystery line. This is the one with the most secure ID. _Or so I thought_ – it turns out that the wavelength is about an Angstrom redder than expected if it is really this line.
#
# **Now I think that [Fe I] 8151.3424 is a more likely explanation**. But we need to look for other components of the multiplet.
#
# * 7708.8389 – this is just to red of a terrestrial absorption feature,  which is probably over-subtracted.  See 05-07. There is something there,  but it looks like sky
# * others are also hidden

# +
em = EmissionLine("si-i-8151", 8151.7)

skyvlim = [100, 200]
skyvlim = [0, 150]
skyrange = SpectralRange(em.wav0, vlim=skyvlim)
siskymoms = {
    label:
    moments.find_moments(
        cube.select_lambda(*skyrange.wavlim)
    )
    for label, cube in csub.items()
}

newvlim = 200, 300
newvlim = 150, 350
nebrange = SpectralRange(em.wav0, vlim=newvlim)
simoms = {
    label:
    moments.find_moments(
        cube.select_lambda(*nebrange.wavlim)
    )
    for label, cube in csub.items()
}
# -

# So we are experimenting with different limits for the ranges.  We will iterate with plotting the 1D spectrum below.
#
# I have defined 
# * `skyrange` as a velocity range where we do not see the line (not necessarily actual sky), giving `siskymoms`
# * `nebrange` as a velocity range where we see the neutral clouds clearly, giving `simoms`

scale = 50
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
siskymoms["D"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 0])
axes[0, 0].set_title("D")
siskymoms["C"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 1])
axes[0, 1].set_title("C")
siskymoms["B"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 0])
axes[1, 0].set_title("B")
im = siskymoms["A"][0].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 1])
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# So it looks like some of this might be sky, especially in panels A and B. But in panel D there is a lot of spatial structure, possibly related to the WR star

# Now try taking the difference `simoms - siskymoms` and rebin at 4x4 to improve s/n

sidiff = {k: (simoms[k][0] - siskymoms[k][0]).rebin(4) for k in "ABCD"}

scale = 50
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(16, 12))
sidiff["D"].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 0])
axes[0, 0].set_title("D")
sidiff["C"].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[0, 1])
axes[0, 1].set_title("C")
sidiff["B"].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 0])
axes[1, 0].set_title("B")
im = sidiff["A"].plot(vmin=-scale, vmax=scale, cmap=cm.fusion, ax=axes[1, 1])
axes[1, 1].set_title("A")
fig.colorbar(im, ax=axes)
...;

# Do the gausssian fit to the _sky_ component and plot the median spectra. Although this is of very little use in this case.

si_sky = {
    label: sky.Sky(cube, em, sky_vlim=skyvlim, full_vlim=[-1000, 1000])
    for label, cube in csub.items()
}


# +
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=False)

for ax, [label, skyline] in zip(axes.flat, si_sky.items()):
    finewavs = np.linspace(
        skyline.wavs.min(), 
        skyline.wavs.max(), 
        200,
    )
    ax.axvspan(*skyline.sky_range.wavlim, alpha=0.1)
    ax.axvspan(*nebrange.wavlim, alpha=0.1, color="g")
    skyline.spec.plot(label=label, linewidth=2, ax=ax)
    ax.plot(finewavs, skyline.gauss(finewavs))
    residuals = skyline.spec.data - skyline.gauss(skyline.wavs)
    ax.scatter(skyline.wavs, residuals)
    ax.legend()
    ax.axvline(skyline.em.wav0, color="r")
    ax.axhline(0, color="r")
fig.tight_layout()
...;
# -

# But this is a bit of a disaster. I can no longer find where the line is supposed to be!!!
#
# But we do know that it is there, we just need to select the right spatial region.  Start off with choosing rectangular boxes in field A since that has the brightest emission:

# I made model sky spectra using the [skycalc](https://www.eso.org/observing/etc/bin/simu/skycalc) webapp. I did separate ones for fields A and C, based on their observation  dates, although there does not seem to be much difference between them. 

# +
from astropy.table import Table

skytabA = Table.read(
    Path.cwd().parent / "data" / "sky" / "skytable-30dorA.fits"
)
skytabC = Table.read(
    Path.cwd().parent / "data" / "sky" / "skytable-30dorC.fits"
)

# +
#widerange = SpectralRange(em.wav0, vlim=[-2700, 2700])
widerange = SpectralRange(em.wav0, vlim=[-4000, 4000])


slice_neutral_A = slice(None), slice(80, 240), slice(240, None)
slice_ionized_A = slice(None), slice(200, None), slice(160, 240)
slice_neutral_C = slice(None), slice(140, 180), slice(80, 120)
slice_ionized_C = slice(None), slice(200, None), slice(80, 160)


spec_neutral_A = (
    csub["A"][slice_neutral_A]
    .select_lambda(*widerange.wavlim)
    .mean(axis=(1, 2))
)
spec_ionized_A = (
    csub["A"][slice_ionized_A]
    .select_lambda(*widerange.wavlim)
    .mean(axis=(1, 2))
)

spec_neutral_C = (
    csub["C"][slice_neutral_C]
    .select_lambda(*widerange.wavlim)
    .mean(axis=(1, 2))
)
spec_ionized_C = (
    csub["C"][slice_ionized_C]
    .select_lambda(*widerange.wavlim)
    .mean(axis=(1, 2))
)
# -

marklines = [
    EmissionLine("[Fe I] 8151", 8151.3424),
    EmissionLine("He I 8156", 8155.71),
    EmissionLine("Ca I 8125", 8125.3),
    EmissionLine("N I 8185", 8184.862),
    EmissionLine("N I 8188", 8188.012),
    EmissionLine("N I 8216", 8216.336),
]
vel0 = 250.0

# +
fig, ax = plt.subplots(figsize=(12, 5))
ax.axvspan(*skyrange.wavlim, alpha=0.1)
ax.axvspan(*nebrange.wavlim, alpha=0.1, color="g")
ax.axhline(0, color="k", lw=1)

factor_A = 0.15
spec_neutral_A.plot(label="neutral region")
spec_ionized_A.plot(alpha=0.3, color="k", label="ionized region")
(spec_ionized_A * factor_A).plot(label=f"ionized x {factor_A}")

skytab_offset = -10
skytab_offset_abs = -14
ax.fill_between(
    skytabA['lam'] * u.nm.to(u.Angstrom),
    (skytabA['flux_ael'] / 1e2) + skytab_offset,
    skytab_offset,
    alpha=0.5,
    color="r",
)
ax.fill_between(
    skytabA['lam'] * u.nm.to(u.Angstrom),
    skytabA['trans_ma'] * 3 + skytab_offset_abs,
    skytab_offset_abs,
    alpha=0.5,
    color="c",
)

for mline in marklines:
    ax.axvline(
        VelocityScale(mline.wav0).vel2wav(vel0),
        ymin=0.5, 
        ymax=0.7,
        color="k",
        linewidth=3,
    )
ax.set(
    ylim=[-15, 35],
#    ylim=[-50, 700],
)
ax.set_title("Field A: ionized vs neutral around 8152 Å")
ax.legend(fontsize="x-small", ncol=3)
fig.tight_layout()
fig.savefig("explore-8152-A.pdf")
...;
# -

mline = marklines[0]
VelocityScale(mline.wav0).vel2wav(vel0)

# The neutral cloud is in blue. The peak ionized emission is in gray, and a scaled version is in orange (scaled so as to match the neutral cloud for the He I 8156  line).
#
# I have also plotted the predicted sky spectrum
#
# So, with this we can identify various types of line:
#
# * If the gray and blue coincide, then this is consistent with a sky line, which should show no spatial variation.
#     * This seems to be confirmed for the two lines either side of 8100, but other cases are not so clear. 
# * If the blue and orange coincide, then it is consistent with an ionized line.
#     * 8156 looks like He I
#     * 8149 may be a line, but I do not know what it is
#     * 8175 is another one

# +
fig, ax = plt.subplots(figsize=(12, 5))
ax.axvspan(*skyrange.wavlim, alpha=0.1)
ax.axvspan(*nebrange.wavlim, alpha=0.1, color="g")
ax.axhline(0, color="k", lw=1)
factor_C = 0.35
spec_neutral_C.plot(label="neutral region")
spec_ionized_C.plot(alpha=0.3, color="k", label="ionized region")
(spec_ionized_C * factor_C).plot(label=f"ionized x {factor_C}")

skytab_offset = -10
skytab_offset_abs = -14
ax.fill_between(
    skytabC['lam'] * u.nm.to(u.Angstrom),
    (skytabC['flux_ael'] / 1e2) + skytab_offset,
    skytab_offset,
    alpha=0.5,
    color="r",
)
ax.fill_between(
    skytabC['lam'] * u.nm.to(u.Angstrom),
    skytabC['trans_ma'] * 3 + skytab_offset_abs,
    skytab_offset_abs,
    alpha=0.5,
    color="c",
)

for mline in marklines:
    ax.axvline(
        VelocityScale(mline.wav0).vel2wav(vel0),
        ymin=0.5, 
        ymax=0.7,
        color="k",
        linewidth=3,
    )

ax.set(
    ylim=[-15, 35],
)
ax.set_title("Field C: ionized vs neutral around 8152 Å")
ax.legend(fontsize="x-small", ncol=3)
fig.tight_layout()
fig.savefig("explore-8152-C.pdf")
...;
# -

[csub[k].data_header["TITLE"] for k in "ABCD"]



# ## More on the deep neutral lines

vel0 = 250.0
VelocityScale(7708.8389).vel2wav(vel0)

# ### I have found a new neutral line around 8790
#
# It is similar to 8151 and the others. So we now have 4 lines.

VelocityScale(8790).vel2wav(vel0)

# ### And another new line at 8851
#
# This one has a slightly different morphology. It looks like it comes from less deep in the PDR, although that may  be dut to the blending with a higher ionization line.

# + tags=[]
VelocityScale(8851).vel2wav(vel0)
# -

# ### And yet another one at 8889
#
# This weaker,  but is definitely there

# + tags=[]
VelocityScale(8889).vel2wav(vel0)
# -

# ### And 8946
#
# Intermediate strength between the weak and the strong. 

# + tags=[]
VelocityScale(8946).vel2wav(vel0)
# -

# ### And 9020
#
# This one is quite strong, although 
# it has an ionized line just to its blue

# + tags=[]
VelocityScale(9020).vel2wav(vel0)
# -

# ### And 9029
#
# Intermediate strength between the weak and the strong. 

# + tags=[]
VelocityScale(9029).vel2wav(vel0)
# -

# Now I have got round to 9114, which we already had. There was also a weak one at 9100 and another at 9057 but that is badly blended with an H line or something. And an extremely weak one maybe at 9035.
#
# Then  after 9114, there is a weak one (NEW) at 9129 before a stronger one at 9146 (already known). There are  no more after that, The sky lines take over. 

# And I found one more too. Very weak at 8284. It is on the red flank of an ionized line (possibly H I).  It seems to have a slightly different distribution from the others. Slightly less deep, similar to 8851.

VelocityScale(9146).vel2wav(vel0)
VelocityScale(9129).vel2wav(vel0)
VelocityScale(8660).vel2wav(vel0)
VelocityScale(8726).vel2wav(vel0)
VelocityScale(8284).vel2wav(vel0)

# ### Ah, no there are more!  I wrapped around and there are some I missed: 8660 is pretty strong
#
# Then we get to C I 8727 and there are no more until I hit 8790 again. 
#
# 8151 8284 8660 8790 8851 8889 8946 9020 9029 9057 9100 9114 9129 9146
#
# That is 14 lines!!! plus C I 8727
#
# And in addition, there is a possible very weak lines at 8810, 8920
#
# Nte that 9057 is also super weak and on the flank of an ionized line. 
#
# Finally, there are intermediate lines, such as 8696, which look deeper than Fe II but shallower than the deep neutral lines. We had flagged it as [Co II] in the Orion project, but I doubt that it really is


