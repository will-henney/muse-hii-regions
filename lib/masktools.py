"""Functions to manipulate masks of MUSE cubes and images

Author: Will Henney, IRyA-UNAM, 2021
"""
from mpdaf.obj import Cube, Image

def trim_edges(im, m):
    """Trim in-place m pixels of each edge of image by setting mask"""
    im.mask[:m, :] = True
    im.mask[-m:, :] = True
    im.mask[:, :m] = True
    im.mask[:, -m:] = True
    return None
