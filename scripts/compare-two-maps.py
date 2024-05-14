from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from pathlib import Path
import typer
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.wcs import WCS

SAVEPATH = Path("maps-compare")


def combo_folder(combo: str):
    cubeid, winid = combo.split("-")
    cubename = {"P": "PZ", "E": "muse"}[cubeid]
    return Path(f"maps-n346-{cubename}-2pass-csub-{winid}")


def line_path(combo: str, line: str):
    candidates = list(combo_folder(combo).glob(f"*/*{line}-ABC.fits"))
    assert len(candidates) == 1
    return candidates[0]


def get_data(line_path: Path, suffix: str = "ABC", return_header: bool = False):
    """Get the image data from a line path"""
    hdu = fits.open(line_path.with_stem(line_path.stem.replace("ABC", suffix)))[0]
    if return_header:
        return hdu.data, hdu.header
    else:
        return hdu.data


def get_zone_spectra(
    combo: str,
    zones: list[str] = ["0", "I", "II", "III", "IV", "S", "MYSO", "BG"],
):
    """Get the 1-d spectrum of each zone from a given cube"""
    specdict = {}
    for zone in zones:
        spec_file = f"zone_spectra/zone-{zone}-{combo}-mean-spec1d.fits"
        hdu = fits.open(spec_file)[0]
        specdict[zone] = hdu.data
        if not "wave" in specdict:
            specdict["wave"] = (
                WCS(hdu.header)
                .array_index_to_world(np.arange(hdu.data.size))
                .to_value("Angstrom")
            )
    # Add some combination zones
    # Average diffuse
    specdict["I,III"] = (specdict["I"] + specdict["III"]) / 2
    # Average filaments
    specdict["0,II"] = (specdict["0"] + specdict["II"]) / 2
    return specdict


def get_zone_masks(
    zones_folder: Path = Path.cwd().parent / "n346-lines",
    zones: list[str] = ["0", "I", "II", "III", "IV", "S", "MYSO", "BG"],
):
    maskdict = {
        zone: fits.open(zones_folder / f"zone-{zone}-mask.fits")[0].data
        for zone in zones
    }
    # Add the combination zones
    maskdict["I,III"] = maskdict["I"] | maskdict["III"]
    maskdict["0,II"] = maskdict["0"] | maskdict["II"]
    return maskdict


def rgb_ABC(line_path: Path):
    rgb = []
    for chan in "CBA":
        data = get_data(line_path, suffix=chan)
        rgb.append(data)
    return np.stack(rgb, axis=-1)


def scale_image(data, vmin, vmax):
    """Linear rescaling from [vmin, vmax] -> [0, 1]"""
    return (data - vmin) / (vmax - vmin)


def auto_scale_channels(rgb, p=1.0, mask=None):
    """Rescale each channel to the 1-99 percentile range."""
    _rgb = np.empty_like(rgb)
    for i in range(3):
        _rgb[..., i] = scale_image(
            rgb[..., i],
            *np.nanpercentile(
                rgb[..., i] if mask is None else rgb[..., i][mask],
                [p, 100 - p],
            ),
        )
    return _rgb


def split_line_string(line: str):
    """Extract the ion and the central wavelength from a line string."""
    parts = line.split("-")
    assert len(parts) in (3, 4)
    species = " ".join(parts[:-3])
    wave0 = float(".".join(parts[-2:]))
    return species, wave0


def main(
    acombo: str = "P-007",
    bcombo: str = "E-007",
    line: str = "h-i-6562-79",
    histogram_gamma: float = 2.0,
    smooth: float = 0.0,
    mask_out_stars: bool = False,
    star_mask_threshold: float = 10.0,
    star_map_path: Path = Path.cwd().parent / "n346-lines" / "zone-S-bright-map.fits",
    zones_folder: Path = Path.cwd().parent / "n346-lines",
    subtract_bg: bool = False,
):
    species, wave0 = split_line_string(line)
    line_path_a = line_path(acombo, line)
    line_path_b = line_path(bcombo, line)
    abc_a, hdr = get_data(line_path_a, return_header=True)
    abc_b = get_data(line_path_b)
    rgb_a = rgb_ABC(line_path_a)
    rgb_b = rgb_ABC(line_path_b)
    wave_obs = float(hdr["lambda_obs"])
    index0 = int(hdr["Index"])
    # Optionally smooth the images
    if smooth > 0.0:
        kernel = Gaussian2DKernel(smooth)
        for i in range(3):
            rgb_a[..., i] = convolve_fft(rgb_a[..., i], kernel)
            rgb_b[..., i] = convolve_fft(rgb_b[..., i], kernel)
            abc_a = convolve_fft(abc_a, kernel)
            abc_b = convolve_fft(abc_b, kernel)
    # Optionally mask out the stars
    if mask_out_stars:
        star_map = fits.open(star_map_path)[0].data
        if smooth > 0.0:
            star_map = convolve_fft(star_map, kernel)
        # star_mask is True when when we have no star
        star_mask = star_map < star_mask_threshold
    else:
        star_mask = np.ones_like(abc_a, dtype=bool)
    star_mask_rgb = np.stack([star_mask] * 3, axis=-1)

    # Get the zone masks
    zone_masks = get_zone_masks(zones_folder)

    # Load the spectra for each zone into dicts
    spec_a = get_zone_spectra(acombo)
    spec_b = get_zone_spectra(bcombo)
    # Optionally subtract the background spectrum from each zone
    if subtract_bg:
        for zone in spec_a.keys():
            if zone not in ("wave", "BG"):
                spec_a[zone] -= spec_a["BG"]
                spec_b[zone] -= spec_b["BG"]

    # Optionally apply the BG subtraction to the images
    if subtract_bg:
        # Subtract from each channel of the RGB images
        for i in range(3):
            # index in wave array for the current channel
            wav_index = index0 + 1 - i
            rgb_a[..., i] -= spec_a["BG"][wav_index]
            rgb_b[..., i] -= spec_b["BG"][wav_index]
        # Subtract from the summed ABC images
        abc_a -= np.sum(spec_a["BG"][index0 - 1 : index0 + 2])
        abc_b -= np.sum(spec_b["BG"][index0 - 1 : index0 + 2])

    # Calculate suitable limits for the plots and histograms
    amin, amax = np.nanpercentile(rgb_a[star_mask_rgb], [1, 99])
    aspan = amax - amin
    amin -= 0.1 * aspan
    amax += 0.1 * aspan
    bmin, bmax = np.nanpercentile(rgb_b[star_mask_rgb], [1, 99])
    bspan = bmax - bmin
    bmin -= 0.1 * bspan
    bmax += 0.1 * bspan
    abmax = max(amax, bmax)
    abmin = min(min(amin, bmin), 0.0)

    # Set up the figure
    sns.set_color_codes("muted")
    fig, ax = plt.subplots(3, 3, figsize=(10, 7))

    # RGB images of the ABC channels from the two cubes
    ax_ima, ax_imb = ax[0, 0], ax[1, 0]
    ax_ima.imshow(auto_scale_channels(rgb_a, mask=star_mask), origin="lower")
    ax_imb.imshow(auto_scale_channels(rgb_b, mask=star_mask), origin="lower")
    if mask_out_stars:
        ax_ima.contour(
            star_map, levels=[star_mask_threshold], colors="r", linewidths=0.5
        )
        ax_imb.contour(
            star_map, levels=[star_mask_threshold], colors="r", linewidths=0.5
        )
    ax_ima.set_title(f"{acombo} {line}")
    ax_imb.set_title(f"{bcombo} {line}")

    # Zones plotted on the map
    ax_map = ax[2, 0]
    for zone, color, alpha in [
        ("0,II", "m", 0.7),
        ("I,III", "c", 0.7),
        ("IV", "y", 0.9),
        ("BG", "0.9", 1.0),
   ]:
        ax_map.contourf(
            zone_masks[zone],
            levels=[0.5, 1.5],
            colors=color,
            alpha=alpha,
        )

    # And over-plot brightness contours
    nlevels = 10
    map = (abc_a + abc_b) / 2.0
    if mask_out_stars:
        map[~star_mask] = np.nan
    levels = np.linspace(0.2, 1.0, nlevels) * np.nanpercentile(map, 99.5) * 1.1
    # levels = np.geomspace(0.01, 1.0, nlevels) * np.nanpercentile(abc_a, 99.5) * 1.1
    linewidths = np.linspace(0.1, 0.7, nlevels)
    ax_map.contour(
        map,
        levels=levels,
        # norm=colors.PowerNorm(gamma=5, vmin=None, vmax=abmax, clip=True),
        # cmap="gray_r",
        origin="lower",
        alpha=1.0,
        linewidths=linewidths,
        colors="k",
    )
    ax_map.set_aspect("equal")

    ## 1D spectra in right column of plots
    ##
    wslice = slice(index0 - 7, index0 + 8)
    smaxima = []
    sminima = []
    for axx, zone in [
        [ax[0, -1], "IV"],
        [ax[1, -1], "I,III"],
        [ax[2, -1], "0,II"],
        # [ax[2, 0], "BG"],
    ]:
        axx.plot(
            spec_a["wave"][wslice],
            spec_a[zone][wslice],
            label=f"{acombo} zone {zone}",
            ds="steps-mid",
        )
        axx.plot(
            spec_b["wave"][wslice],
            spec_b[zone][wslice],
            label=f"{bcombo} zone {zone}",
            ds="steps-mid",
        )
        # Save the highest and lowest values for the y-axis limits
        if zone != "BG":
            smaxima.append(np.nanmax(spec_a[zone][wslice]))
            smaxima.append(np.nanmax(spec_b[zone][wslice]))
            sminima.append(np.nanmin(spec_a[zone][wslice]))
            sminima.append(np.nanmin(spec_b[zone][wslice]))

        # Add the zero line
        axx.axhline(0.0, color="k", lw=0.5)
        # add RGB rectangles for the ABC channels
        dwave = np.diff(spec_a["wave"])[0]
        for i, color in ([index0 - 1, "b"], [index0, "g"], [index0 + 1, "r"]):
            axx.axvspan(
                spec_a["wave"][i] - dwave / 2,
                spec_a["wave"][i] + dwave / 2,
                facecolor=color,
                alpha=0.2,
                lw=0,
            )
        axx.set_title(f"Zone {zone} ", loc="right", y=0.8)

    # Set common plot limits for all spectra
    smax = max(smaxima)
    smin = min(sminima)
    sspan = smax - smin
    smin -= 0.1 * sspan
    smax += 0.1 * sspan
    for axx in ax[:, -1]:
        axx.set_ylim(smin, smax)

    ## Correlations in middle column
    # Correlations between the two cubes, channel by channel
    ax_corr_rgb = ax[2, 1]
    nbins = 100
    H_rgb = np.empty((nbins, nbins, 3))
    for ichan in range(3):
        x = rgb_a[..., ichan]
        y = rgb_b[..., ichan]
        mask = np.isfinite(x) & np.isfinite(y)
        mask &= (x != 0.0) & (y != 0.0)
        if mask_out_stars:
            mask &= star_mask
        H, xedges, yedges = np.histogram2d(
            x[mask],
            y[mask],
            weights=x[mask] - abmin,
            bins=nbins,
            range=[[abmin, abmax], [abmin, abmax]],
        )
        # if ichan == 1:
        #     print(xedges)
        #     print(yedges)
        #     print(H)
        H_rgb[..., ichan] = H.T

    ax_corr_rgb.imshow(
        auto_scale_channels(H_rgb, p=0) ** (1 / histogram_gamma),
        origin="lower",
        extent=[abmin, abmax, abmin, abmax],
    )
    # Add marginal rgb distributions along the axes
    for ichan, color in enumerate(["r", "g", "b"]):
        Hx = np.sum(H_rgb[..., ichan], axis=0) / np.sum(H_rgb[..., ichan])
        Hy = np.sum(H_rgb[..., ichan], axis=1) / np.sum(H_rgb[..., ichan])
        ax_corr_rgb.stairs(
            values=Hx * 2 * (abmax - abmin) + abmin,
            edges=xedges,
            orientation="vertical",
            baseline=abmin,
            fill=True,
            alpha=0.5,
            color=color,
            lw=0.5,
        )
        ax_corr_rgb.stairs(
            values=Hy * 2 * (abmax - abmin) + abmin,
            edges=yedges,
            orientation="horizontal",
            baseline=abmin,
            fill=True,
            alpha=0.5,
            color=color,
            lw=0.5,
        )

    ax_corr_rgb.axhline(0.0, color="w", lw=0.5, linestyle="--")
    ax_corr_rgb.axvline(0.0, color="w", lw=0.5, linestyle="--")
    ax_corr_rgb.plot([abmin, abmax], [abmin, abmax], color="w", lw=0.5, linestyle="--")
    ax_corr_rgb.set(
        xlabel=f"{acombo} {line}",
        ylabel=f"{bcombo} {line}",
    )

    # Recalculate the moments, since we may have smoothed the arrays
    m1_a = (rgb_a[..., 0] - rgb_a[..., 2]) / abc_a
    m1_b = (rgb_b[..., 0] - rgb_b[..., 2]) / abc_b
    m2_a = (rgb_a[..., 0] + rgb_a[..., 2]) / abc_a
    m2_b = (rgb_b[..., 0] + rgb_b[..., 2]) / abc_b

    # Correlations between the velocity moments
    for axx, mlabel, mrange, mlines in [
        [ax[0, 1], "m1", (-0.7, 0.7), (-0.5, 0.0, 0.5)],
        [ax[1, 1], "m2", (-0.2, 1.2), (0, 2 / 3)],
    ]:
        # Joint histogram of velocity moments
        if mlabel == "m1":
            m_a = m1_a
            m_b = m1_b
        else:
            m_a = m2_a
            m_b = m2_b
        # m_a = get_data(line_path_a, suffix=mlabel)
        # m_b = get_data(line_path_b, suffix=mlabel)
        mask = np.isfinite(m_a) & np.isfinite(m_b)
        mask &= abc_a > 0.0
        if mask_out_stars:
            mask &= star_mask
        mmin, mmax = mrange
        H, xedges, yedges = np.histogram2d(
            m_a[mask],
            m_b[mask],
            weights=abc_a[mask],
            bins=nbins,
            range=[[mmin, mmax], [mmin, mmax]],
        )
        axx.imshow(
            H.T ** (1 / histogram_gamma),
            origin="lower",
            extent=[mmin, mmax, mmin, mmax],
            cmap="gray_r",
        )
        axx.plot([mmin, mmax], [mmin, mmax], color="r", lw=0.5, linestyle="--")
        for mline in mlines:
            axx.axhline(mline, color="r", lw=0.5, linestyle="--")
            axx.axvline(mline, color="r", lw=0.5, linestyle="--")
            axx.set(
                xlabel=f"{acombo} {mlabel}",
                ylabel=f"{bcombo} {mlabel}",
            )
        # Add marginal distributions along the axes
        Hx = np.sum(H.T, axis=0) / np.sum(H)
        Hy = np.sum(H.T, axis=1) / np.sum(H)
        axx.stairs(
            values=Hx * 1 * (mmax - mmin) + mmin,
            edges=xedges,
            orientation="vertical",
            baseline=mmin,
            fill=True,
            alpha=0.5,
            color="r",
            lw=0.5,
        )
        axx.stairs(
            values=Hy * 1 * (mmax - mmin) + mmin,
            edges=yedges,
            orientation="horizontal",
            baseline=mmin,
            fill=True,
            alpha=0.5,
            color="r",
            lw=0.5,
        )

    sns.despine()
    figfile = SAVEPATH / f"{acombo}-{bcombo}-{line}.pdf"
    fig.savefig(figfile, bbox_inches="tight")

    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
