"""Utilities for keeping regenerated figure assets on the manuscript canvas."""
from pathlib import Path

from PIL import Image


def resize_png_canvas(path, width, height, *, background="white", dpi=300):
    """Fit a rendered PNG to a fixed pixel canvas without distorting its content.

    The image is uniformly resized to fit and centered on a white canvas. This
    keeps spatial axes and markers at their correct aspect ratio while making
    the published PNG dimensions deterministic.
    """
    path = Path(path)
    with Image.open(path) as source:
        image = source.convert("RGB")
    scale = min(width / image.width, height / image.height)
    size = (max(1, round(image.width * scale)), max(1, round(image.height * scale)))
    if size != image.size:
        image = image.resize(size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", (width, height), background)
    offset = ((width - size[0]) // 2, (height - size[1]) // 2)
    canvas.paste(image, offset)
    temp = path.with_name(path.stem + ".canvas.tmp.png")
    canvas.save(temp, format="PNG", dpi=(dpi, dpi), optimize=True)
    temp.replace(path)
    return path
