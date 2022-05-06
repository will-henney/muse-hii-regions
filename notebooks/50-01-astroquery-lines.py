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

# # Possible IDs for the deep neutral lines
#
# These are the lines, whose spatial distribution is similar to C I 8727 and the H_2 2.12 micron line. 
#
# ## Lines around the neutral 8152 line
#
# This is the shortest wavelength of these lines, which we are lookin at first.  We are also interested in other lines in the same wavelength range, even if they are form the ionized gas. 
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

# ## Trying to find IDs for the lines

# First look for all the neutral lines of common elements that fall in the wavelength range where we see the lines.
#
#
# ### Forbidden neutral lines of all the common elements
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

# And use it to select just the neutral ions, since the API does not allow this to be done on the query. 

res[res["stage"]=="I"].show_in_notebook()

# #### Confident ID: [Fe I] 8151.3424
#
# I am pretty sure of this one. 

# #### Other potential matches
#
# So there are some partial overlaps with the observed wavelengths here. 
#
# For instance, 9020 might be [Co I] 9019.65
#
# But most of the wavelengths are misses.
#
#

(1 * u.eV).to(1/u.cm, equivalencies=u.equivalencies.spectral())

# ### Permitted lines of selected neutral atoms
#
# Next thing to try is to expand search to include permitted lines as well as forbidden ones. If they are not from autoionizing levels, then they should have upper level energies less than the ionization potential of the neutral atom, which for these elements is about 6 eV.
#
# I have put a limit of 4 eV to start with, so I do not get too many hits. 
#
# The excitation could be collisions, fluorescence, or recombination.  In the case of collisions, the upper level cannot be very high since w have a very low T, probably hundreds of K rather than thousands. 
#
# But for the C I line, recombinations are thought to dominate. If the same is true of the other lines, then the brightness will be proportional to the abundance, more or less I should think. 
#
# The C I 8727 line is not a resonance line, and is 10 times weaker than C I 9850 (which is), according to Cesarsky 1982A&A...113L...7C
#
# This means there is room for observable recombination lines from less abundant elements, but they will be most likely resonant lines since those are the ones that capture the greatest fraction of the total recombination rate. 

res = AtomicLineList.query_object(
    wavelength_range=wavelength_range,
    wavelength_type='Air',
    upper_level_energy_range=4 * u.eV,
    element_spectrum="Fe",
    #transitions=aa.Transition.nebular,
)
res["stage"] = [s.strip("[]").split()[-1] for s in res["SPECTRUM"]]

res[res["stage"]=="I"].show_in_notebook()

# So that is with an upper limit of 4 eV for the energy, which gives 39 candidate lines of Fe I. If I use 7 eV instead, then I get 10 times as many lines. But then, the question is: why do we see only some of those (perhaps they vary greatly in strength, but the database does not have A values for the majority).
#
# This gives a potential match for the 9146 line of Fe I 9146.1285, but not for anything else. 

# We can try it again for the other likely elements.
#
# Try looking first at resonance lines of neutral atoms. I am restricting the lower energy level to less than 0.1 eV, so that will include the fine structure splitting of the ground configuration. 
#

res = AtomicLineList.query_object(
    wavelength_range=wavelength_range,
    wavelength_type='Air',
    lower_level_energy_range=0.1 * u.eV,
    #element_spectrum="Fe",
    #transitions=aa.Transition.nebular,
)
res["stage"] = [s.strip("[]").split()[-1] for s in res["SPECTRUM"]]

# + tags=[]
res[res["stage"]=="I"].show_in_notebook()
# -

# So that gives no good matches, apart from the Fe I line that we already knew about. 

# Just to make sure, we eill look at singly ionized species too, just in case there is one with an extremely low ionization potential.

res[res["stage"]=="II"].show_in_notebook()

# Certainly no candidates for our neutral lines, although we do pick up the [Cl II] lines at 8578.69 and 9123.6, which we certainly observe, although the shorte one is blended with a He I line.

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

# ## The 6104 high ionization line

wavelength_range = (6102 * u.Angstrom, 1 * u.Angstrom)
#wavelength_range = (6500 * u.Angstrom, 6600 * u.Angstrom)
#AtomicLineList.FORM_URL = "https://www.pa.uky.edu/~peter/newpage/"
res = AtomicLineList.query_object(
    wavelength_range=wavelength_range,
    wavelength_type='Air',
    transitions=aa.Transition.nebular,
)

res.show_in_notebook()



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


