from matplotlib import pyplot as plt
import matplotlib as mpl
from astropy.io import fits
import numpy as np
from typing import Union
import typer

DPI = 150

def main(
        fits_file: str,
        min_value: float=0.0,
        max_value: float=6.0,
        cmap: str="RdBu",
        nan_color: str="0.5",
        fig_file: str="",
):
    hdu = fits.open(fits_file)[0]
    ny, nx = hdu.data.shape
    fig, ax = plt.subplots(figsize=(nx/DPI, ny/DPI))
    ax.imshow(
        hdu.data,
        origin="lower",
        interpolation="none",
        vmin=min_value,
        vmax=max_value,
        cmap=mpl.colormaps[cmap].with_extremes(bad=nan_color),
    )
    ax.set(xticks=[], yticks=[])
    fig.subplots_adjust(0.0, 0.0, 1.0, 1.0)
    if not fig_file:
        fig_file = fits_file.replace(".fits", ".png")
    fig.savefig(fig_file, dpi=DPI)
    print(fig_file, end="")


if __name__ == "__main__":
    typer.run(main)
