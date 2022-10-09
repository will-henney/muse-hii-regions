from  typing import Union
from mpdaf.obj import Spectrum, Cube, Image
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import typer

def main(
        cube_file: str,
        peak_file: str,
        star_mask_file: Union[str, None]=None,
        vsys: float=170.0,
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
    ncolumns = 15
    nrows = int(np.ceil(nlines / ncolumns))
    pane_size = 2
    fig, axes = plt.subplots(
        nrows, ncolumns,
        figsize=(ncolumns * pane_size, nrows * pane_size),
    )
    suffix = peak_file.replace(".csv", ".pdf")
    figfile = f"peak-images-{suffix}"

    for row, ax in zip(tab, axes.flat):
        wav1 = cube.wave.coord(row["left_ips"])
        wav2 = cube.wave.coord(row["right_ips"])
        im = cube.select_lambda(wav1, wav2).sum(axis=0)
        if star_mask_file is not None:
            im.mask = im.mask | star_mask
        scale = row["prominences"]
        if scale < 2.0:
            im = im.rebin(4)
        elif scale < 8.0:
            im = im.rebin(2)
        im.plot(ax=ax, vmin=-0.5*scale, vmax=5*scale)
        wav0 = 0.5 * (wav1 + wav2) / (1 + vsys / 300000)
        title = f"{wav0:.2f} {scale:.4g}"
        ax.set_title(title)

    fig.tight_layout()
    fig.savefig(figfile)
    print(figfile, end="")


if __name__ == "__main__":
    typer.run(main)
