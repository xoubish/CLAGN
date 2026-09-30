"""Embed only rendered figure pixels, never a clipped full literature page."""
from io import BytesIO
from pathlib import Path
from tempfile import TemporaryDirectory
import math
import subprocess

from PIL import Image
from pypdf import PdfReader, PdfWriter
from pypdf.generic import (
    DecodedStreamObject, DictionaryObject, NameObject, NumberObject,
)


def cropped_figure(page, box, dpi=600):
    """Return a lossless image-only PDF page; box is top-down page points.

    Render only the requested rectangle before constructing a fresh PDF.
    No source text, fonts, annotations, or off-crop image pixels are copied.
    """
    left, top, right, bottom = box
    width, height = right-left, bottom-top
    assert width > 0 and height > 0
    factor = dpi/72
    x, y = math.floor(left*factor), math.floor(top*factor)
    w, h = math.ceil(right*factor)-x, math.ceil(bottom*factor)-y
    with TemporaryDirectory(prefix='jwst-figure-crop-') as tmp:
        tmp = Path(tmp)
        source = PdfWriter()
        source.add_page(page)
        source.write(tmp/'source.pdf')
        subprocess.run([
            'pdftoppm', '-r', str(dpi), '-x', str(x), '-y', str(y),
            '-W', str(w), '-H', str(h), '-png', '-singlefile',
            str(tmp/'source.pdf'), str(tmp/'crop'),
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        with Image.open(tmp/'crop.png') as image:
            image = image.convert('RGB')
            pixels = DecodedStreamObject()
            pixels.set_data(image.tobytes())
            pixels.update({
                NameObject('/Type'): NameObject('/XObject'),
                NameObject('/Subtype'): NameObject('/Image'),
                NameObject('/Width'): NumberObject(image.width),
                NameObject('/Height'): NumberObject(image.height),
                NameObject('/ColorSpace'): NameObject('/DeviceRGB'),
                NameObject('/BitsPerComponent'): NumberObject(8),
            })
    writer = PdfWriter()
    result = writer.add_blank_page(width=width, height=height)
    result[NameObject('/Resources')] = DictionaryObject({
        NameObject('/XObject'): DictionaryObject({
            NameObject('/FigureCrop'): writer._add_object(pixels.flate_encode()),
        }),
    })
    drawing = DecodedStreamObject()
    drawing.set_data(f'q {width} 0 0 {height} 0 0 cm /FigureCrop Do Q'.encode())
    result[NameObject('/Contents')] = writer._add_object(drawing)
    output = BytesIO()
    writer.write(output)
    output.seek(0)
    return PdfReader(output).pages[0]
