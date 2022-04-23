# -*- coding: utf-8 -*-
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

from astropy import units as u
from astroquery.atomic import AtomicLineList

# # Look at the supposed Ca I lines
#
# ## The 7890 line
#
# This looks a lot like the [Fe III] lines

wavelength_range = (7890 * u.Angstrom, 0.2 * u.Angstrom)
#wavelength_range = (6500 * u.Angstrom, 6600 * u.Angstrom)
#AtomicLineList.FORM_URL = "https://www.pa.uky.edu/~peter/newpage/"
res = AtomicLineList.query_object(
    wavelength_range=wavelength_range,
    wavelength_type='Air',
    minimal_abundance="6",
    depl_factor="0",
    #element_spectrum='H I',
    #get_query_payload=True,
)

res

# So it looks like the [Ni III] 7889.9 line is a very good candidate.  It is a nebular forbidden line, so would be expected to be strong. 


