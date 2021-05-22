# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light,md
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

# # NGC 346 stars from Sabbi et al (2007)
#
# In my never-ending quest to get a catalog of all the stars in NGC 346, I am now going to try using the astroquery Vizier interface to get the tables from: 
#
# ```bibtex
# @article{Sabbi:2007h,
# 	author = {{Sabbi}, E. and {Sirianni}, M. and {Nota}, A. and {Tosi}, M. and {Gallagher}, J. and {Meixner}, M. and {Oey}, M.~S. and {Walterbos}, R. and {Pasquali}, A. and {Smith}, L.~J. and {Angeretti}, L.},
# 	journal = {\aj},
# 	month = jan,
# 	number = {1},
# 	pages = {44-57},
# 	title = {{Past and Present Star Formation in the SMC: NGC 346 and its Neighborhood}},
# 	volume = {133},
# 	year = 2007}
#
# ```

# ## Find which catalog we want
#
# I had to install the `astroquery` library from conda-forge.

from astroquery.vizier import Vizier

# I am following a similar method to the examples in the [astroquery docs](https://astroquery.readthedocs.io/en/latest/vizier/vizier.html):

catalog_list = Vizier.find_catalogs('Sabbi 2007')

# We now have a list of matching catalogs:

for k,v in catalog_list.items():
    print(k)
    print(v.description, end="\n\n")

# This gave more than I was expecting, but the one we want is `J/AJ/133/44`:

catalog_id = "J/AJ/133/44"

# ## Download the tables
#
# We have to set `Vizier.ROW_LIMIT = -1` to make sure we get the complete table.  Otherwise, the default is to return only 50 rows.

Vizier.ROW_LIMIT = -1
catalogs = Vizier.get_catalogs(catalog_id)
catalogs

# This gives two tables.  The first is a list of all 79960 stars that they have found in NGC 346.  The second is a list of the 16 sub-clusters that they have identified:

stars, subclusters = catalogs
subclusters

# ## Look at the star catalog

stars

# We want to convert the coordinate columns into degrees for ease of plotting, so we need some astropy libraries:

from astropy.coordinates import SkyCoord, Angle, Latitude, Longitude
import astropy.units as u

# I had originally thought of using SkyCoord, but this is not necessary

# +
# SkyCoord(stars["RAJ2000"], stars["DEJ2000"], unit=(u.hourangle, u.deg))
# -

# Easier to just convert the RA and DEC columns separately (note, these cells take a while to run – 10s of seconds on my laptop):

ra = Longitude(stars["RAJ2000"], unit=u.hourangle)

dec = Latitude(stars["DEJ2000"], unit=u.deg)

# Check that the values in degrees look reasonable:

ra.to(u.deg), dec

# Add them in as new columns to the table and set display precision to about 10 mas. 

stars["ra"] = ra.to(u.deg)
stars["dec"] = dec
stars["ra"].format = "{:.6f}"
stars["dec"].format = "{:.6f}"
stars


