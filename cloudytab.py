from astropy.table import Table  # type: ignore
from astropy.io.ascii import InconsistentTableError  # type: ignore
from pathlib import Path

# File extensions that might be present, but which are NOT Cloudy save files
IGNORE_EXTENSIONS = ["pdf", "png", "jpg"]

# Input and output files, which are ingested literally
IO_EXTENSIONS = ["in", "out"]


class CloudyModel:
    """
    Lightweight wrapper for output from Cloudy run

    For example:

    >>> from cloudytab import CloudyModel
    >>> m = CloudyModel("myfolder/mymodel")

    `m.files` contains a list of all the files that were found:
              `['myfolder/mymodel.in', 'myfolder/mymodel.ovr', ETC]`

    `m.data` contains dict of astropy.Table's, one for each save file:
              `{'ovr': <Table length=289> ..., ETC}`

    `m.io['in']` and `m.io['out']` contain the input and output streams

    `m.skipped` contains a dict of each extension that was skipped with its reason.
    """

    def __init__(self, prefix: str):
        self.filepaths = Path(".").glob(f"{prefix}.*")
        self.data = {}
        self.io = {}
        self.skipped = {}
        for filepath in self.filepaths:
            saveid = filepath.suffix.strip(".")
            if saveid in IGNORE_EXTENSIONS:
                # Figure files, etc need to be skipped
                self.skipped[saveid] = "Extension is listed in IGNORE_EXTENSIONS"
            elif saveid in IO_EXTENSIONS:
                # Special case of input and output files
                with open(filepath) as f:
                    # Just save the whole file as a string
                    self.io[saveid] = f.read()
            else:
                # Assume all else are save files
                try:
                    self.data[saveid] = Table.read(
                        str(filepath),
                        delimiter="\t",
                        guess=False,
                        fast_reader=False,
                        format="ascii.commented_header",
                    )
                except UnicodeDecodeError:
                    # Binary files can raise this error - ignore them
                    self.skipped[saveid] = "Table.read() raised UnicodeDecodeError"
                except InconsistentTableError:
                    # The "save heating" files can raise this error - skip them
                    self.skipped[saveid] = "Table.read() raised InconsistentTableError"
