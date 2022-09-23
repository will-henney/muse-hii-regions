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

# + [markdown] pycharm={"name": "#%% md\n"}
# # How region masks work
#
# These are the regions defiend in the `regions` library. They can be read from DS9 files or created by hand. There are versions in pixel and celestial sky coordinates, but here I will concentrate only on the pixel regions.

# + pycharm={"name": "#%%\n"}
from pathlib import Path

import regions as rg
import seaborn as sns
from matplotlib import pyplot as plt
from mpdaf.obj import Image

sns.set_context("talk")

# + [markdown] pycharm={"name": "#%% md\n"}
# ## TLDR recommendations
#
# For 2d images that are not too large, then we can just use
# 1. `reg.to_mask().to_image(image)` which returns a logical mask of the same shape as the image
# 2. or `reg.to_mask().multiply(image)` which directly applies the multiplication
#
# For 3d cubes, we can use `slices_large` which is the first of two results returned from `reg.to_mask().get_overlap_slices(shape)`. This is the slices into the image dimensions of the cube that extracts the bounding box of the region. Then we can extend that along the spectral axis. For non-rotated rectangular regions, this is enough. But for more general shapes, we can also apply it to the result of `.to_image()`

# + pycharm={"name": "#%%\n"}
small_data_path = Path.cwd().parent.parent / "data"
region_file = "ngc346-jesus-muse-pixels.reg"
pixel_regions = rg.Regions.read(small_data_path / region_file)

# + pycharm={"name": "#%%\n"}
line_map = Image(str(small_data_path / "ngc346-rgbchan-xxx-8151-sum.fits"))

# + pycharm={"name": "#%%\n"}
region_dict = {reg.meta["label"]: reg for reg in pixel_regions}

# + pycharm={"name": "#%%\n"}
vmin, vmax = 0.0, 100.0

# + pycharm={"name": "#%%\n"}
fig, ax = plt.subplots(figsize=(8, 8))
line_map.plot(vmin=vmin, vmax=vmax, ax=ax, cmap="gray_r", title="Positions of sample regions")
for label, reg in region_dict.items():
    is_bg = label.endswith(" bg")
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
# Look at one particular region in detail and create a region mask object:

# + pycharm={"name": "#%%\n"}
region = region_dict["FIL"]
region_mask = region.to_mask()

# + [markdown] pycharm={"name": "#%% md\n"}
# We can get two sets of slices. The first is the more useful, which is the slices into the full  image array that corresponds to the bounding box of the mask. The second is the same, but into the small mask cutout array, which is only necessary if the mask goes outside of the image.

# + pycharm={"name": "#%%\n"}
slices_large, slices_small = region_mask.get_overlap_slices(line_map.data.shape)
slices_large, slices_small

# + [markdown] pycharm={"name": "#%% md\n"}
# We will check that `slices_small` is just the full size of the cutout rectangle in this case:

# + pycharm={"name": "#%%\n"}
region_mask.cutout(line_map.data).shape

# + [markdown] pycharm={"name": "#%% md\n"}
# Now we will look at the different methods on the mask.

# + pycharm={"name": "#%%\n"}
fig, axes = plt.subplots(1, 3, sharex='all', sharey='all',  figsize=(12, 3))
axes[0].imshow(
    region_mask.cutout(line_map.data),
    vmin=vmin, vmax=vmax, origin="lower", cmap="gray_r")
axes[0].set_title("mask.cutout()")
axes[1].imshow(
    region_mask.to_image(line_map.data.shape)[slices_large],
    vmin=0, vmax=1.0, origin="lower", cmap="gray_r")
axes[1].set_title("mask.to_image()[slices_large]")
axes[2].imshow(
    region_mask.multiply(line_map.data),
    vmin=vmin, vmax=vmax, origin="lower", cmap="gray_r")
axes[2].set_title("mask.multiply()")
...;

# + [markdown] pycharm={"name": "#%% md\n"}
# How to extend this to cubes

# + pycharm={"name": "#%%\n"}
(slice(None, None),) + slices_large

# + pycharm={"name": "#%%\n"}

