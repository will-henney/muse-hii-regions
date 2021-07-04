from astropy.table import Table
from astropy.io.ascii import InconsistentTableError
import glob

# File extensions that might be present, but which are NOT Cloudy save files
IGNORE_EXTS = ["pdf", "png", "jpg"]


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

    `m.skipped` contains a dict of each extensions that was skipped with its reason.
    """

    def __init__(self, prefix: str):
        self.files = glob.glob(prefix + ".*")
        self.data = {}
        self.io = {}
        self.skipped = {}
        for file_ in self.files:
            saveid = file_.split(".")[-1]
            if saveid in IGNORE_EXTS:
                # Figure files, etc need to be skipped
                self.skipped[saveid] = "Extension is listed in IGNORE_EXTS"
            elif saveid in ["in", "out"]:
                # Special case of input and output files
                with open(file_) as f:
                    # Just save the whole file as a string
                    self.io[saveid] = f.read()
            else:
                # Assume all else are save files
                try:
                    self.data[saveid] = Table.read(
                        file_,
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
