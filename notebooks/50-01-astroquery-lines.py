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
import astroquery.atomic as a
from astropy.table import conf

# # Lines around the neutral 8152 line
#
# There are certainly some He lines there.

wavelength_range = (8152 * u.Angstrom, 70 * u.Angstrom)
res = AtomicLineList.query_object(
    wavelength_range=wavelength_range,
    wavelength_type='Air',
    minimal_abundance="6",
    depl_factor="0",
    element_spectrum="\n".join(['H', 'He']),
    #get_query_payload=True,
)

res

wavelength_range = (8152 * u.Angstrom, 70 * u.Angstrom)
ionlist = " I-IV\n".join(
    [
        "C", 
        "N", 
        "O", 
        "S",
        "Si",
        "Cl",
        "Fl",
        "Fe",
        "Ne",
        "Ar",
        "Ni",
        "Ca",
        "Mg",
    ]
)
res = AtomicLineList.query_object(
    wavelength_range=wavelength_range,
    wavelength_type='Air',
    lower_level_energy_range=3 * u.eV,
    #transitions=aa.Transition.nebular,
    element_spectrum=ionlist,
    #get_query_payload=True,
)

res

# +
wavelength_range = (8152 * u.Angstrom, 70 * u.Angstrom)
ionlist = "\n".join([
    "N I", "N II", 
    "C I", "C II", 
    "O I", "O II", 
    "Si I", "Si II", 
    "Ca I", "Ca II", 
])

res = AtomicLineList.query_object(
    wavelength_range=wavelength_range,
    wavelength_type='Air',
    #lower_level_energy_range=3 * u.eV,
    #transitions=aa.Transition.nebular,
    element_spectrum=ionlist,
    #get_query_payload=True,
)
res.show_in_notebook()
# -

# # All the [Fe I] lines (and others)

# First look for all the neutral lines of common elements that fall in the wavelength range where we see the lines.
#
# Start off with just the forbidden lines, so we do not get too many. Also, impose an upper energy level cutoff of less than 4 eV, which would make sense if they wer excited collisionally.

wavelength_range = (
    8100 * u.Angstrom, 9200 * u.Angstrom
)
elements = "\n".join(
    [
        "C", "N", "O",
        "Mg", "Ca", "Na", 
        "S", "Si", "P",
        "F", "Cl",
        "Ar", "Xe",
        "Fe", "Ni", "Co",
    ]
)
res = AtomicLineList.query_object(
    wavelength_range=wavelength_range,
    wavelength_type='Air',
    upper_level_energy_range=10 * u.eV,
    element_spectrum=elements,
    transitions=aa.Transition.nebular,
)

# Add a new column that gives the ion stage: I, II, III, IV, etc

res["stage"] = [s.strip("[]").split()[-1] for s in res["SPECTRUM"]]

# And use it to selecct just the neutral ions, since the API does not allow this to be done on the query. 

res[res["stage"]=="I"].show_in_notebook()

# So there are some partial overlaps with the observed wavelengths here. 
#
# For instance, 9020 might be [Co I] 9019.65
#
# But most of the wavelengths are misses.
#
#

(1 * u.eV).to(1/u.cm, equivalencies=u.equivalencies.spectral())

# Next thing to try is to expand search to include permitted lines as well as forbidden ones. If they are excited by fluorescence of recombination, then they should have excitation energies less than the ionization potential of the neutral atom, which for these elements is about 6 eV.  

res = AtomicLineList.query_object(
    wavelength_range=wavelength_range,
    wavelength_type='Air',
    upper_level_energy_range=7 * u.eV,
    element_spectrum="Fe",
    #transitions=aa.Transition.nebular,
)
res["stage"] = [s.strip("[]").split()[-1] for s in res["SPECTRUM"]]

res[res["stage"]=="I"].show_in_notebook()



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

# +
# AtomicLineList.query_object_async??
# -

aa.Transition.nebular

simple_transitions = [aa.Transition.all, aa.Transition.nebular]

simple_transitions

"Neb" in simple_transitions

aa.utils.is_valid_transitions_param("Neb")

aa.utils.is_valid_transitions_param(aa.Transition.nebular)

from astropy import units as u


# +
# u.spectral??

# +
# u.Equivalency?
# -


