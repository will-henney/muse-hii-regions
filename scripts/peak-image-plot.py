from  typing import Union
from mpdaf.obj import Spectrum, Cube, Image
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import to_hex
from matplotlib.patches import BoxStyle
import matplotlib
import seaborn as sns
import typer

def corners_text(
        ax: matplotlib.axes.Axes,
        texts: tuple[str, str, str, str],
        pad: float=0.05,
        **kwds,
) -> None:
    """Write 4 texts to 4 corners of graph"""
    ax.text(pad, 1 - pad, texts[0],
            ha="left", va="top", transform=ax.transAxes, **kwds)
    ax.text(1 - pad, 1 - pad, texts[1],
            ha="right", va="top", transform=ax.transAxes, **kwds)
    ax.text(pad, pad, texts[2],
            ha="left", va="bottom", transform=ax.transAxes, **kwds)
    ax.text(1 - pad, pad, texts[3],
            ha="right", va="bottom", transform=ax.transAxes, **kwds)


def main(
        cube_file: str,
        peak_file: str,
        star_mask_file: Union[str, None]=None,
        vsys: float=170.0,
        ncolumns: int=15,
        subtract_base: bool=False,
        wavelength_window_pad: float=1.5,
        scale_by_percentile: bool=True,
        use_rainbow_colors: bool=True,
        rainbow_saturation: float=95.0,
        rainbow_lightness: float=70.0,
        rainbow_blue_red: tuple[float, float]=(270.0, -45.0),
        extra_suffix: Union[str, None]=None,
):
    """Plot of images for each peak"""

    tab = Table.read(
        peak_file,
        format="ascii.ecsv",
    )
    cube = Cube(cube_file)
    if star_mask_file is not None:
        star_mask = Image(star_mask_file).data > 0.0

    nlines = len(tab)
    nrows = int(np.ceil(nlines / ncolumns))
    pane_size = 2
    fig, axes = plt.subplots(
        nrows, ncolumns,
        figsize=(ncolumns * pane_size, nrows * pane_size),
    )
    if extra_suffix is not None:
        suffix = peak_file.replace(".csv", f"-{extra_suffix}.pdf")
    else:
        suffix = peak_file.replace(".csv", ".pdf")
    figfile = f"peak-images-{suffix}"

    wavmin, wavmax = cube.wave.get_range()
    for row, ax in zip(tab, axes.flat):
        # Take the half-maximum wave range and expand it by the padding on both sides
        wav1 = cube.wave.coord(row["left_ips"]) - wavelength_window_pad
        wav2 = cube.wave.coord(row["right_ips"]) + wavelength_window_pad
        cube_window = cube.select_lambda(wav1, wav2)
        if subtract_base:
            # Subtract off average of channels from left and right bases
            ib1, ib2 = row["left_bases"], row["right_bases"],
            av_bases_image = 0.5 * (cube.data[ib1, :, :] + cube.data[ib2, :, :])
            cube_window.data -= av_bases_image[None, :, :]
        # Sum the window to get the extracted line image
        im = cube_window.sum(axis=0)
        if star_mask_file is not None:
            im.mask = im.mask | star_mask
        scale = row["prominences"]

        # For weak lines, rebin the pixels to get better s/n
        if scale < 0.5:
            im = im.rebin(8)
        elif scale < 2.0:
            im = im.rebin(4)
        elif scale < 8.0:
            im = im.rebin(2)

        # Brightness scaling
        if scale_by_percentile:
            # First, take 5th to 95th percentile span
            vmin, vmax = np.percentile(im.data[~im.mask], [5, 95])
            vspan = vmax - vmin
            # Then extend to by a certain fraction above and below
            vmin -= 0.0 * vspan
            vmax += 0.5 * vspan
        else:
            # If not using percentiles, just use the prominence data from find_peaks()
            vmin, vmax = -0.5*scale, 5*scale

        # Color map
        if use_rainbow_colors:
            # Fractional distance between blue and red ends of spectrum
            xwav = (row["Wavelength"] - wavmin) /  (wavmax - wavmin)
            # Convert to hue angle
            blue, red = rainbow_blue_red
            rainbow_hue = blue + (red - blue) * xwav
            # Make a nice color map using this Hue as a key color
            cmap = sns.light_palette(
                (rainbow_hue, rainbow_saturation, 100 - rainbow_lightness),
                input="husl",
                as_cmap=True,
            )
        else:
            cmap = "gray_r"

        im.plot(ax=ax, vmin=vmin, vmax=vmax, cmap=cmap.with_extremes(bad="0.5"))
        wav0 = 0.5 * (wav1 + wav2) / (1 + vsys / 300000)
        labels = (
            f"λ{wav0:.2f}",
            f"{vmax:.4g}", 
            f"#{row['Pixel']:04d}",
            f"δλ{row['widths']:.1f}",
        )
        corners_text(ax, labels, pad=0.02,
                     color="k", fontweight="bold",
                     bbox=dict(
                         facecolor="w",
                         boxstyle=BoxStyle.Round(pad=0.1),
                         alpha=0.6,
                     ),
                     )

    for ax in axes.flat:
        ax.set(xticks=[], yticks=[])
    sns.despine(left=True, bottom=True)
    fig.tight_layout(h_pad=0.2, w_pad=0.2)
    fig.savefig(figfile)
    print(figfile, end="")




if __name__ == "__main__":
    typer.run(main)
