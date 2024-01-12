# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
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

# + pycharm={"name": "#%%\n"}
from astropy.table import Table
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import yaml
from astropy.modeling.models import BlackBody
import astropy.units as u
from astropy import constants


# + pycharm={"name": "#%%\n"}
sns.set_context("talk")
sns.set_color_codes()

# + pycharm={"name": "#%%\n"}
datapath = Path.cwd().parent.parent / "data-jesus"

# + pycharm={"name": "#%%\n"}
regions = "BS", "YSO", "bkg", "bkg2"
bands = "SL1", "SL2", "LL1", "LL2"
def _load_spectrum(region, band):
    try:
        return Table.read(
            str(datapath / "NGC346_spec" / f"ngc346_{region}_final_{band}.tbl"),
            format="ascii.ipac"
        )
    except FileNotFoundError:
        return None

data = {
    region: {
        band: _load_spectrum(region, band)
        for band in bands
    }
    for region in regions
}

# + pycharm={"name": "#%%\n"}
data["bkg2"]

# + pycharm={"name": "#%%\n"}
line_list = yaml.safe_load(open(datapath / "spitzer-lines.yaml"))
line_list


# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(
    figsize=(20, 8),
)
for region in regions:
    for band in bands:
        if data[region][band]:
            ax.plot("WAVELENGTH", "FLUX", data=data[region][band])
for linedata in line_list:
    if 'elow' in linedata and linedata['elow'] > 1e4:
        # Highly excited configuration
        continue
    elif linedata['label'].split()[-1].startswith("V") or linedata['label'] == "He II":
        # Highly ionized
        linewidth = 0.5
        y0 = 2.0
        fontsize = 'xx-small'
    elif linedata['label'].startswith('PAH'):
        # PAH band
        linewidth = 2.0
        y0 = 300.0
        fontsize = 'small'
    elif linedata['label'].startswith('H'):
        # Hydrogen, Helium, or H_2 line
        linewidth = 1.5
        y0 = 100.0
        fontsize = 'x-small'
    else:
        # Ground configuration
        linewidth = 1.0
        y0 = 10.0
        fontsize = 'x-small'
    ax.axvline(linedata['wave'], linestyle='dashed', color='r', linewidth=linewidth)
    ax.text(
        linedata['wave'], y0, f"{linedata['label']} {linedata['wave']:.3f}",
        ha="center", va="center", rotation="vertical", fontsize=fontsize,
    )

ax.set(
    yscale="log",
    ylim=[0.3, 500],
    ylabel="Flux",
    xlabel="Wavelength, micron",
)
fig.savefig('jesus-spitzer-spectra.pdf', bbox_inches='tight')
...;

# + [markdown] pycharm={"name": "#%% md\n"}
#     ### Take ratios with respect to the background

# + pycharm={"name": "#%%\n"}
for region in regions:
    for band in bands:
        data[region][band]["RATIO"] = data[region][band]["FLUX"] / data["bkg2"][band]["FLUX"]

# + pycharm={"name": "#%%\n"}
fig, [ax, axx] = plt.subplots(
    2, 1, sharex=True,
    figsize=(20, 8),
)
for band in bands:
    ax.plot("WAVELENGTH", "RATIO", data=data["BS"][band])
    axx.plot("WAVELENGTH", "RATIO", data=data["YSO"][band])
for a in ax, axx:
    a.axhline(1.0, color='k', linestyle='dashed')
ax.set(
    ylim=[0., 2.9],
    ylabel="BS / BG",
)
axx.set(
    ylim=[0., 16],
    ylabel="YSO / BG",
    xlabel="Wavelength, micron",
)
fig.savefig('jesus-spitzer-ratios.pdf', bbox_inches='tight')
...;

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Subtract BG
#
# Mixture of the two BG regions.  Aim is to get the [S III] lines to disappear on the BG-subtracted BS spectrum

# + pycharm={"name": "#%%\n"}
mix2 = 0.7
mix1 = 1 - mix2

for region in regions:
    for band in bands:
        bg = mix1 * data["bkg"][band]["FLUX"] + mix2 * data["bkg2"][band]["FLUX"]
        data[region][band]["BGSUB"] = data[region][band]["FLUX"] - bg

# + [markdown] pycharm={"name": "#%% md\n"}
# Compare with  a black body
#

# + pycharm={"name": "#%%\n"}
((u.cm / u.s) * u.MJy / u.sr / u.micron ** 2).to(u.erg / u.s / u.cm**2 / u.sr / u.micron)
#(u.MJy / u.sr).to(u.erg / u.s / u.cm**2 / u.sr / u.Hz)

# + pycharm={"name": "#%%\n"}
waves = np.linspace(5.0, 40.0, 200) * u.micron
bb = BlackBody(
    temperature=135 * u.K,
    scale=5e-15 * constants.c * u.MJy / u.sr / u.micron ** 2,
)

# + pycharm={"name": "#%%\n"}
fig, [ax, axx, axxx] = plt.subplots(
    3, 1, sharex=True,
    figsize=(10, 12),
)
for band in bands:
    d = data["BS"][band]
    ax.plot(d["WAVELENGTH"], d["BGSUB"] / d["WAVELENGTH"]**2, color='r')
    d = data["YSO"][band]
    axx.plot(d["WAVELENGTH"], d["BGSUB"] / d["WAVELENGTH"]**2, color='m')
    d = data["bkg"][band]
    axxx.plot(d["WAVELENGTH"], d["FLUX"] / d["WAVELENGTH"]**2, color='g')
    d = data["bkg2"][band]
    axxx.plot(d["WAVELENGTH"], d["FLUX"] / d["WAVELENGTH"]**2, color='y')
ax.plot(waves, bb(waves) / constants.c, color="k", linestyle="dotted")
for a in ax, axx, axxx:
    a.axhline(0.0, color='k', linestyle='dashed')
ax.set(
    ylim=[-0.05, None],
    ylabel=r"$F_\lambda$ (BS $-$ BG)",
)
axx.set(
    ylim=[None, None],
    ylabel=r"$F_\lambda$ (YSO $-$ BG)",
)
axxx.set(

    ylim=[None, None],
    ylabel=r"$F_\lambda$ (BG)",
    xlabel="Wavelength, micron",
)
fig.savefig('jesus-spitzer-bgsub.pdf', bbox_inches='tight')
...;

# + pycharm={"name": "#%%\n"}
bb.lambda_max.to(u.micron)

# + pycharm={"name": "#%%\n"}
(constants.c / bb.nu_max).to(u.micron)

# + [markdown] pycharm={"name": "#%% md\n"}
# Use modified BB using Cloudy opacities

# + pycharm={"name": "#%%\n"}
opac_tab = Table.read("../../data/xsec-infrared-dust-silicate_ism_10.ecsv")
wavgrid = opac_tab["Wavelength"].data * u.micron
kappa_s = opac_tab["Opacity"].data
opac_tab = Table.read("../../data/xsec-infrared-dust-graphite_ism_10.ecsv")
fac = 1.5
kappa_c = fac * opac_tab["Opacity"].data

# + pycharm={"name": "#%%\n"}
bb_s = BlackBody(
    temperature=140 * u.K,
    scale=2.8e-15 * constants.c * u.MJy / u.sr / u.micron ** 2,
)
bb_ss = BlackBody(
    temperature=350 * u.K,
    scale=2.8e-15 * constants.c * u.MJy / u.sr / u.micron ** 2,
)
bb_c = BlackBody(
    temperature=120 * u.K,
    scale=2.8e-15 * constants.c * u.MJy / u.sr / u.micron ** 2,
)
mbb_sed_s = wavgrid * kappa_s * bb_s(wavgrid) / constants.c
mbb_sed_ss = wavgrid * (kappa_s + 4 * kappa_c) * bb_ss(wavgrid) / constants.c
hfac = 0.0008
mbb_sed_ss *= hfac
mbb_sed_c = wavgrid * kappa_c * bb_c(wavgrid) / constants.c
fac = 1.5
mbb_sed_c *= fac
mbb_sed = mbb_sed_s + mbb_sed_ss + mbb_sed_c
mbb_sed_c *= 7 / mbb_sed.max()
mbb_sed_s *= 7 / mbb_sed.max()
mbb_sed_ss *= 7 / mbb_sed.max()
mbb_sed *= 7 / mbb_sed.max()

# + [markdown] pycharm={"name": "#%% md\n"}
# Or the same, but the SED

# + pycharm={"name": "#%%\n"}
fig, [ax, axx, axxx] = plt.subplots(
    3, 1, sharex=True,
    figsize=(8, 9),
)
for band in bands:
    d = data["BS"][band]
    ax.plot(d["WAVELENGTH"], d["BGSUB"] / d["WAVELENGTH"], color='r')
    d = data["YSO"][band]
    axx.plot(d["WAVELENGTH"], d["BGSUB"] / d["WAVELENGTH"], color='m')
    d1 = data["bkg"][band]
    d2 = data["bkg2"][band]
    axxx.plot(d["WAVELENGTH"], (mix1 * d1["FLUX"] + mix2 * d2["FLUX"]) / d["WAVELENGTH"], color='g')
ax.plot(waves, waves * bb(waves) / constants.c, color="k", linestyle="dotted")
# axx.plot(wavgrid, mbb_sed, color="k", linestyle="dotted")
# axx.plot(wavgrid, mbb_sed_c, color="r", linestyle="dotted")
# axx.plot(wavgrid, mbb_sed_s, color="c", linestyle="dotted")
# axx.plot(wavgrid, mbb_sed_ss, color="y", linestyle="dotted")

# axxx.plot(waves, 1e3 * waves * bb2(waves) / constants.c, color="k", linestyle="dotted")
ax.set(
    ylim=[-0.5, None],
    ylabel=r"$\lambda F_\lambda$ (BS $-$ BG)",
)
axx.set(
    ylim=[-0.5, None],
    ylabel=r"$\lambda F_\lambda$ (YSO $-$ BG)",
)
axxx.set(
    xlim=[3.5, 41.5],
    ylim=[-0.5, None],
    yticks=[0, 2, 4, 6],
    ylabel=r"$\lambda F_\lambda$ (BG)",
    xlabel=r"Wavelength $\lambda$, micron",
)
# for a in ax, axx, axxx:
#     a.set_yscale("symlog", linthresh=0.5, linscale=0.1)
#     a.set_xscale("log")
fig.savefig('jesus-spitzer-bgsub-sed.pdf', bbox_inches='tight')
...;

# + [markdown] pycharm={"name": "#%% md\n"}
# In the end I did not need any if the modified bb curves for the YSO spectrum, since none of them were much good.  It is possible to get the dip between the two silicate peaks fitted pretty well, but not at the same time as the long and short wavelengths
#

# + [markdown] pycharm={"name": "#%% md\n"}
# ### Look at all the dust opacity files

# + pycharm={"name": "#%%\n"}
opac_files = (Path.cwd().parent.parent / "data").glob('xsec-infrared-dust-*.ecsv')
opacities = {p.stem.split('-')[-1]: Table.read(str(p)) for p in opac_files}
opacities

# + pycharm={"name": "#%%\n"}
dtypes = "silicate_ism_10", "graphite_ism_10"
dweights = 0.1, 1.0
fig, ax = plt.subplots(figsize=(12, 10))
opac_sum = np.zeros_like(kappa)
for label, weight in zip(dtypes, dweights):
    ax.plot(wavgrid, weight * opacities[label]["Opacity"], label=f"{label} x {weight:.1f}")
    opac_sum += weight * opacities[label]["Opacity"]
ax.plot(wavgrid, opac_sum, label="Total")

ax.legend()
ax.set(xlim=[4, 42], ylim=[0.0, 3], yscale='linear')
...;

# + pycharm={"name": "#%%\n"}

