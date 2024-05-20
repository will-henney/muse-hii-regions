from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from pathlib import Path
import typer
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.wcs import WCS
import astropy.constants as const  # type: ignore
import astropy.units as u  # type: ignore

SAVEPATH = Path("maps-compare")
LIGHT_SPEED_KMS = const.c.to(u.km / u.s).value

# husl tuples (hue, saturation, lightness) for each zone
KEY_COLORS = {
    "IV": (50, 90, 50),  # orange-yellow
    "I,III": (200, 85, 60),  # blue-green
    "0,II": (275, 80, 40),  # purple
}

ZONE_TEXTS = {
    "IV": "Bow shock",
    "I,III": "Diffuse",
    "0,II": "Filaments",
}

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
    vsys: float,
    zones: list[str] = ["0", "I", "II", "III", "IV", "S", "MYSO", "BG"],
):
    """Get the 1-d spectrum of each zone from a given cube"""
    specdict = {}
    for zone in zones:
        spec_file = f"zone_spectra/zone-{zone}-{combo}-mean-spec1d.fits"
        hdu = fits.open(spec_file)[0]
        if not "wave" in specdict:
            specdict["wave"] = (
                WCS(hdu.header)
                .array_index_to_world(np.arange(hdu.data.size))
                .to_value("Angstrom")
            )
        specdict[zone] = hdu.data

    # Add some combination zones
    # Average diffuse
    specdict["I,III"] = (specdict["I"] + specdict["III"]) / 2
    # Average filaments
    specdict["0,II"] = (specdict["0"] + specdict["II"]) / 2
    # Average non-filaments - to use for correcting hyper-local continuum
    specdict["I,III,IV"] = (specdict["I,III"] + specdict["IV"]) / 2
    # Correct the velocities to nebular frame
    specdict["wave"] /= 1 + vsys / LIGHT_SPEED_KMS
    return specdict


def get_moments_from_zone_spectra(specdict: dict[str, np.ndarray], ipeak: int):
    """Get the moments of the spectra from a dictionary of zone spectra"""
    moments = {}
    for zone, spec in specdict.items():
        if zone == "wave":
            continue
        M0 = np.sum(spec[ipeak - 1 : ipeak + 2])
        moments[zone] = {
            "sum": M0,
            "m1": (spec[ipeak + 1] - spec[ipeak - 1]) / M0,
            "m2": (spec[ipeak + 1] + spec[ipeak - 1]) / M0,
        }
    return moments


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


def get_trim_mask(im, margin=5):
    """Return mask that trims the edges off an image by a given margin"""
    mask = np.zeros_like(im, dtype=bool)
    mask[margin:-margin, margin:-margin] = True
    return mask


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
    trim_edges: int = 0,
    vsys: float = 171.1,
    fix_continuum: bool = False,
    hyper_local_blue: int = 3,
    hyper_local_red: int = 3,
    mark_moments_bow_shock: bool = True,
    mark_moments_nebula: bool = True,
    mark_moments_filaments: bool = True,
    debug: bool = False,
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

    # Optionally trim the edges
    if trim_edges > 0:
        trim_mask = get_trim_mask(abc_a, margin=trim_edges)
        trim_mask_rgb = np.stack([trim_mask] * 3, axis=-1)
        # Apply the trim mask to the images
        abc_a[~trim_mask] = np.nan
        abc_b[~trim_mask] = np.nan
        rgb_a[~trim_mask_rgb] = np.nan
        rgb_b[~trim_mask_rgb] = np.nan

    # Optionally smooth the images
    if smooth > 0.0:
        kernel = Gaussian2DKernel(smooth)
        for i in range(3):
            rgb_a[..., i] = convolve_fft(rgb_a[..., i], kernel, preserve_nan=True)
            rgb_b[..., i] = convolve_fft(rgb_b[..., i], kernel, preserve_nan=True)
        abc_a = convolve_fft(abc_a, kernel, preserve_nan=True)
        abc_b = convolve_fft(abc_b, kernel, preserve_nan=True)
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

    # Get the mask for each of the spatial zones
    zone_masks = get_zone_masks(zones_folder)

    # Load the spectra for each zone into dicts
    spec_a = get_zone_spectra(acombo, vsys)
    spec_b = get_zone_spectra(bcombo, vsys)
    # Optionally subtract the background spectrum from each zone
    if subtract_bg:
        for zone in spec_a.keys():
            if zone not in ("wave", "BG"):
                spec_a[zone] -= spec_a["BG"]
                spec_b[zone] -= spec_b["BG"]

    # Optionally fix the hyper-local continuum level using average spectrum from diffuse zones
    if fix_continuum:
        # We use the pixels that are 3 to the left/right of the peak
        hyper_local_a = (
            spec_a["I,III,IV"][index0 - hyper_local_blue]
            + spec_a["I,III,IV"][index0 + hyper_local_red]
        ) / 2
        hyper_local_b = (
            spec_b["I,III,IV"][index0 - hyper_local_blue]
            + spec_b["I,III,IV"][index0 + hyper_local_red]
        ) / 2
        if debug:
            print(f"Hyper-local continuum: {hyper_local_a:.2f}, {hyper_local_b:.2f}")
        for zone in spec_a.keys():
            # Very important not to subtract continuum from BG zone,
            # otherwise the BG subtraction will cancel out the
            # continuum subtraction for the other zones (2024-05-17
            # spent whole day figuring this out!!!)
            if zone not in ("wave", "BG"):
                spec_a[zone] -= hyper_local_a
                spec_b[zone] -= hyper_local_b

    # Get the moments of the spectra
    moments_spec_a = get_moments_from_zone_spectra(spec_a, index0)
    moments_spec_b = get_moments_from_zone_spectra(spec_b, index0)

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

    # Optionally apply the hyper-local continuum fix to the images
    if debug:
        rgb_a_orig = rgb_a.copy()
        abc_a_orig = abc_a.copy()
    if fix_continuum:
        if debug:
            print(f"Subtracting hyper-local continuum from RGB and ABC images")
        rgb_a -= hyper_local_a
        rgb_b -= hyper_local_b
        abc_a -= 3 * hyper_local_a
        abc_b -= 3 * hyper_local_b

    # Sanity check that the summed image is still equal to the sum of the ABC channels
    mm = np.isfinite(abc_a)
    check1 = abc_a[mm]
    check2 = np.nansum(rgb_a, axis=-1)[mm]
    # assert np.allclose(check1, check2, atol=0.1, rtol=1e-3), np.sort(np.abs(check1 - check2))[-3:]

    # Recalculate the moments, since we may have smoothed the arrays
    # and/or subtracted the background and hyper-local continuum
    m1_a = (rgb_a[..., 0] - rgb_a[..., 2]) / abc_a
    m1_b = (rgb_b[..., 0] - rgb_b[..., 2]) / abc_b
    m2_a = (rgb_a[..., 0] + rgb_a[..., 2]) / abc_a
    m2_b = (rgb_b[..., 0] + rgb_b[..., 2]) / abc_b

    if debug:
        m2_a_orig = (rgb_a_orig[..., 0] + rgb_a_orig[..., 2]) / abc_a_orig
        print(
            "Median m2 =", np.nanmedian(m2_a), "original m2 =", np.nanmedian(m2_a_orig)
        )

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

    # Save the recalibrated images and spectra

    # Set up the figure
    sns.set_color_codes("muted")
    fig, ax = plt.subplots(3, 3, figsize=(10, 8.5))

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
        if zone in KEY_COLORS:
            color = sns.light_palette(KEY_COLORS[zone], input="husl", as_cmap=True)(0.5)
        ax_map.contourf(
            zone_masks[zone],
            levels=[0.5, 1.5],
            colors=[color, color],
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
        cmap = sns.light_palette(KEY_COLORS[zone], input="husl", as_cmap=True)
        axx.plot(
            spec_a["wave"][wslice],
            spec_a[zone][wslice],
            label=f"{acombo}",
            ds="steps-mid",
            color=cmap(0.5),
            lw=1.5,
        )
        axx.plot(
            spec_b["wave"][wslice],
            spec_b[zone][wslice],
            label=f"{bcombo}",
            ds="steps-mid",
            color=cmap(1.0),
            lw=1.0,
        )
        # Save the highest and lowest values for the y-axis limits
        if zone != "BG":
            # Consider only the ABC pixels for the maximum
            smaxima.append(np.nanmax(spec_a[zone][index0 - 1 : index0 + 2]))
            smaxima.append(np.nanmax(spec_b[zone][index0 - 1 : index0 + 2]))
            # Consider all pixels for the minimum
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
        # Add an indication of the rest wavelength
        axx.axvline(wave0, color="k", lw=0.5, ls="dotted")
        axx.legend(loc="upper left", fontsize="small")
        axx.set_title(f"{ZONE_TEXTS[zone]} ", loc="right", y=0.8, color=cmap(1.0))
        axx.minorticks_on()
        # Center the viewport on the rest wavelength
        axx.set_xlim(wave0 - 10, wave0 + 10)

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

    # Correlations between the velocity moments
    ax_m1, ax_m2 = ax[0, 1], ax[1, 1]
    for axx, mlabel, mrange, mlines, ticks in [
        [ax_m1, "m1", (-0.7, 0.7), (-0.5, 0.0, 0.5), (-0.5, 0.0, 0.5)],
        [ax_m2, "m2", (-0.2, 1.2), (0, 2 / 3), (0, 0.5, 1)],
    ]:
        # Joint histogram of velocity moments
        if mlabel == "m1":
            m_a = m1_a
            m_b = m1_b
        else:
            m_a = m2_a
            m_b = m2_b
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
        # Optionally mark the moments of the average spectrum for each zone
        for zone, mark_moments in [
            ["I,III", mark_moments_nebula],
            ["0,II", mark_moments_filaments],
            ["IV", mark_moments_bow_shock],
        ]:
            if mark_moments:
                mzone_a = moments_spec_a[zone][mlabel]
                mzone_b = moments_spec_b[zone][mlabel]
                color = color = sns.light_palette(
                    KEY_COLORS[zone], input="husl", as_cmap=True
                )(0.7)
                axx.plot(
                    mzone_a, mzone_b, marker="+", color=color, mew=0.5, ms=10, alpha=1.0
                )

        # Ensure same tick marks on x, y axes
        axx.set_xticks(ticks)
        axx.set_yticks(ticks)
        axx.minorticks_on()
        axx.minorticks_on()

    sns.despine()
    figfile = SAVEPATH / f"{acombo}-{bcombo}-{line}.pdf"

    # Set spacing between the panels. We need it to be consistent
    # between different runs, so we cannot use tight_layout
    plt.subplots_adjust(
        left=0.07,
        bottom=0.07,
        right=0.97,
        top=0.95,
        wspace=0.3,
        hspace=0.3,
    )
    # fig.tight_layout()
    fig.savefig(figfile)

    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
