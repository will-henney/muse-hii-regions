"""
Routines taken from:

Hubble Source Catalog API Notebook: SMC Color-Magnitude Diagram August
2019, Rick White

https://archive.stsci.edu/hst/hsc/help/api/hscv3_smc_api.html

This doc string by William Henney 2021-05-21.  All the code is from
Rick White's notebook (see above), but I have tidied it up with black. 

Dependencies: requests.py
"""


import requests, json
from astropy.table import Table

hscapiurl = "https://catalogs.mast.stsci.edu/api/v0.1/hsc"


def hsccone(
    ra,
    dec,
    radius,
    table="summary",
    release="v3",
    format="csv",
    magtype="magaper2",
    columns=None,
    baseurl=hscapiurl,
    verbose=False,
    **kw
):
    """Do a cone search of the HSC catalog

    Parameters
    ----------
    ra (float): (degrees) J2000 Right Ascension
    dec (float): (degrees) J2000 Declination
    radius (float): (degrees) Search radius (<= 0.5 degrees)
    table (string): summary, detailed, propermotions, or sourcepositions
    release (string): v3 or v2
    magtype (string): magaper2 or magauto (only applies to summary table)
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'numimages.gte':2)
    """

    data = kw.copy()
    data["ra"] = ra
    data["dec"] = dec
    data["radius"] = radius
    return hscsearch(
        table=table,
        release=release,
        format=format,
        magtype=magtype,
        columns=columns,
        baseurl=baseurl,
        verbose=verbose,
        **data
    )


def hscsearch(
    table="summary",
    release="v3",
    magtype="magaper2",
    format="csv",
    columns=None,
    baseurl=hscapiurl,
    verbose=False,
    **kw
):
    """Do a general search of the HSC catalog (possibly without ra/dec/radius)

    Parameters
    ----------
    table (string): summary, detailed, propermotions, or sourcepositions
    release (string): v3 or v2
    magtype (string): magaper2 or magauto (only applies to summary table)
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'numimages.gte':2).  Note this is required!
    """

    data = kw.copy()
    if not data:
        raise ValueError("You must specify some parameters for search")
    if format not in ("csv", "votable", "json"):
        raise ValueError("Bad value for format")
    url = "{}.{}".format(cat2url(table, release, magtype, baseurl=baseurl), format)
    if columns:
        # check that column values are legal
        # create a dictionary to speed this up
        dcols = {}
        for col in hscmetadata(table, release, magtype)["name"]:
            dcols[col.lower()] = 1
        badcols = []
        for col in columns:
            if col.lower().strip() not in dcols:
                badcols.append(col)
        if badcols:
            raise ValueError(
                "Some columns not found in table: {}".format(", ".join(badcols))
            )
        # two different ways to specify a list of column values in the API
        # data['columns'] = columns
        data["columns"] = "[{}]".format(",".join(columns))

    # either get or post works
    # r = requests.post(url, data=data)
    r = requests.get(url, params=data)

    if verbose:
        print(r.url)
    r.raise_for_status()
    if format == "json":
        return r.json()
    else:
        return r.text


def hscmetadata(table="summary", release="v3", magtype="magaper2", baseurl=hscapiurl):
    """Return metadata for the specified catalog and table

    Parameters
    ----------
    table (string): summary, detailed, propermotions, or sourcepositions
    release (string): v3 or v2
    magtype (string): magaper2 or magauto (only applies to summary table)
    baseurl: base URL for the request

    Returns an astropy table with columns name, type, description
    """
    url = "{}/metadata".format(cat2url(table, release, magtype, baseurl=baseurl))
    r = requests.get(url)
    r.raise_for_status()
    v = r.json()
    # convert to astropy table
    tab = Table(
        rows=[(x["name"], x["type"], x["description"]) for x in v],
        names=("name", "type", "description"),
    )
    return tab


def cat2url(table="summary", release="v3", magtype="magaper2", baseurl=hscapiurl):
    """Return URL for the specified catalog and table

    Parameters
    ----------
    table (string): summary, detailed, propermotions, or sourcepositions
    release (string): v3 or v2
    magtype (string): magaper2 or magauto (only applies to summary table)
    baseurl: base URL for the request

    Returns a string with the base URL for this request
    """
    checklegal(table, release, magtype)
    if table == "summary":
        url = "{baseurl}/{release}/{table}/{magtype}".format(**locals())
    else:
        url = "{baseurl}/{release}/{table}".format(**locals())
    return url


def checklegal(table, release, magtype):
    """Checks if this combination of table, release and magtype is acceptable

    Raises a ValueError exception if there is problem
    """

    releaselist = ("v2", "v3")
    if release not in releaselist:
        raise ValueError(
            "Bad value for release (must be one of {})".format(", ".join(releaselist))
        )
    if release == "v2":
        tablelist = ("summary", "detailed")
    else:
        tablelist = ("summary", "detailed", "propermotions", "sourcepositions")
    if table not in tablelist:
        raise ValueError(
            "Bad value for table (for {} must be one of {})".format(
                release, ", ".join(tablelist)
            )
        )
    if table == "summary":
        magtypelist = ("magaper2", "magauto")
        if magtype not in magtypelist:
            raise ValueError(
                "Bad value for magtype (must be one of {})".format(
                    ", ".join(magtypelist)
                )
            )


def mastQuery(request, url="https://mast.stsci.edu/api/v0/invoke"):
    """Perform a MAST query.

    Parameters
    ----------
    request (dictionary): The MAST request json object
    url (string): The service URL

    Returns the returned data content
    """

    # Encoding the request as a json string
    requestString = json.dumps(request)
    r = requests.post(url, data={"request": requestString})
    r.raise_for_status()
    return r.text


def resolve(name):
    """Get the RA and Dec for an object using the MAST name resolver

    Parameters
    ----------
    name (str): Name of object

    Returns RA, Dec tuple with position
    """

    resolverRequest = {
        "service": "Mast.Name.Lookup",
        "params": {"input": name, "format": "json"},
    }
    resolvedObjectString = mastQuery(resolverRequest)
    resolvedObject = json.loads(resolvedObjectString)
    # The resolver returns a variety of information about the resolved object,
    # however for our purposes all we need are the RA and Dec
    try:
        objRa = resolvedObject["resolvedCoordinate"][0]["ra"]
        objDec = resolvedObject["resolvedCoordinate"][0]["decl"]
    except IndexError as e:
        raise ValueError("Unknown object '{}'".format(name))
    return (objRa, objDec)
