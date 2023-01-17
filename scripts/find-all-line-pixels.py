import numpy as np
import yaml
from pathlib import Path
import typer

INDEX_PATTERN = "[0-9]" * 4
IGNORED_TYPES = ["noise", "stellar"]

def main(
        data_dir: str="all-lines-orig",
):
    # Get all the lines
    line_files = sorted(Path(data_dir).glob(f"{INDEX_PATTERN}.yaml"))
    data = [
        yaml.safe_load(path.open()) for path in line_files
    ]
    # Make a list of all indices that we wish to exclude from BG
    # estimation
    indices = []
    indices_core = []
    for line in data:
        if not line["Type"]:
            # Skip case that Type is null
            continue
        if line["Type"].lower().lstrip("?") in IGNORED_TYPES:
            # Skip other Types that we do not want, e.g., Noise
            continue                
        idx = line["Index"]
        # Use a 3-pixel window around each line
        indices.extend([idx - 1, idx, idx + 1])
        indices_core.append(idx)

    # Save to a file as a list of integers
    np.savetxt(
        Path(data_dir) / "line-indices.txt",
        indices,
        fmt="%04d",
    )
    np.savetxt(
        Path(data_dir) / "line-indices-core.txt",
        indices_core,
        fmt="%04d",
    )

if __name__ == "__main__":
    typer.run(main)
