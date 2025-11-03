from mpdaf.obj import Image
from pathlib import Path

ROOT = Path.cwd().parent.parent


def get_image_raw(
    lineid: str,
    variant: str = "csub",
    channel: str = "ABC",
) -> Image:
    p = ROOT / "data" / "n346-bow-lines" / f"maps-n346-PZ-2pass-{variant}-007"
    candidates = list(p.glob(f"**/*-{lineid}-{channel}.fits"))
    assert candidates
    assert len(candidates) == 1, f"Multiple candidates: {candidates}"
    imfile = candidates[0]
    im = Image(str(imfile))
    im.data = im.data.astype(float)
    return im
