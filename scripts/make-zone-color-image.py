from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns
from astropy.io import fits
import numpy as np
import yaml
from typing import Union
import typer

DPI = 150

def main(
        zone_file: str="zones.yaml",
        saturation: int=95,
        lightness: int=30,
):
    with open(zone_file) as f:
        zones = yaml.safe_load(f)

    fig = None
    for zone in zones:
        # Load key brightness map for this zone
        hdu = fits.open(f"{zone['label']}-bright-map.fits")[0]
        # Load zone mask
        mhdu = fits.open(f"{zone['label']}-mask.fits")[0]
        # Set pixels outside the mask to NaN
        image = np.where(mhdu.data, hdu.data, np.nan)
        if fig is None:
            ny, nx = hdu.data.shape
            fig, ax = plt.subplots(figsize=(nx/DPI, ny/DPI))

        cmap = sns.light_palette(
            tuple(zone["husl"]),
            input="husl",
            as_cmap=True,
        )
        ax.imshow( image, origin="lower", interpolation="none",
                   vmin=-1.0, vmax=zone["max_bright"], cmap=cmap, )
    ax.set(xticks=[], yticks=[])
    fig.subplots_adjust(0.0, 0.0, 1.0, 1.0)
    fig_file = "zone-color-image.png"
    fig.savefig(fig_file, dpi=DPI)
    print(fig_file, end="")


if __name__ == "__main__":
    typer.run(main)
