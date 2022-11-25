from mpdaf.obj import Spectrum
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import pandas as pd
from pathlib import Path
import yaml
import typer
import astropy.constants as const  # type: ignore
import astropy.units as u  # type: ignore

LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value

STAR_ZONES = ("MYSO")
def main(
        id_label: str,
        wave_range: tuple[float, float],
        separation: float=5.0,
        star_scale: float=5.0,
        v_sys: float=171.1,
        zone_file: str="zones.yaml",
):
    """Plot of spectra from zones in a given wave range"""
    fig, ax = plt.subplots(figsize=(6, 6.4))
    figfile = f"spec1d-{id_label}-{int(wave_range[0])}-{int(wave_range[1])}.pdf"

    with open(zone_file) as f:
        zones = yaml.safe_load(f)
    df = pd.read_csv(f"all-lines-{id_label}/uil-final-table.csv").set_index("Index")
    doppler =  1.0 + v_sys / LIGHT_SPEED_KMS
    nzones = len(zones)
    dy = separation
    offset = 0.0
    for zone in reversed(zones):
        specfile = f"{zone['label']}-{id_label}-spec1d.fits"
        spec = Spectrum(f"zone_spectra/{specfile}")
        # Now select the wavelength range
        spec = spec.subspec(*wave_range)
        # MYSO needs reducing in height
        scale = star_scale if zone["label"].endswith(STAR_ZONES) else 1.0
        # Choose a slightly darker version of the key color
        color = sns.dark_palette(
            tuple(zone["husl"]),
            input="husl",
            as_cmap=True,
        )(0.6)
        ax.plot(
            1e10 * spec.wave.coord() / doppler,
            spec.data / scale + offset,
            drawstyle="steps-mid",
            linewidth=0.7,
            color=color,
        )
        # (spec / scale + offset).plot(label=zone["label"], linewidth=1, color=color)
        ax.axhline(offset, linewidth=0.5, color=color)
        label = zone["label"].split("-")[-1]
        if zone["label"].endswith(STAR_ZONES):
            label = label + f" / {int(scale)}"
        ax.text(
            wave_range[1] + 5, offset,
            label,
            ha="left", va="center", color=color,
        )
        offset += dy
    ax.text(
        wave_range[1] + 5, offset,
        "Zone",
        ha="left", va="center", color="k",
    )
    for data in df.itertuples():
        if wave_range[0] / doppler < data.wave0 < wave_range[1] / doppler:
            linestyle = "dotted" if data.blend else "solid"
            if data.Type.startswith('Deep'):
                ax.axvline(data.wave0, 0.85, 0.95, color="0.9", lw=1.5, linestyle=linestyle)
            else:
                ax.axvline(data.wave0, 0.85, 0.9, color="0.9", lw=1, linestyle=linestyle)
    ax.minorticks_on()
    ax.yaxis.set_tick_params(which='minor', left=False)
    ax.grid(which="major", linewidth=0.5)
    ax.set_ylim(-3 * separation, 9 * separation)
    ax.set_xlim(None, wave_range[1])
    ax.set_xlabel("STP wavelength in rest frame of nebula, Å")
    ax.set_ylabel(
        "Continuum-subtracted mean brightness plus offset\n"
        r"$10^{-20}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$ pix$^{-1}$"
    )
    sns.despine()
    fig.savefig(figfile, bbox_inches="tight")
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
