#!/usr/bin/env python3
"""
strip_xlsx_drawings.py — remove orphan drawing relationships from openxlsx-written
workbooks so the .xlsx validates and opens without an Excel "repair" prompt.

openxlsx writes, per worksheet, a relationships entry pointing at a drawing part
(xl/drawings/drawingN.xml) and a legacy VML drawing (vmlDrawingN.vml) that it never
actually creates, and declares the drawing part in [Content_Types].xml. The sheet
bodies do not reference these, so removing the orphan relationships, the drawing
content-type overrides, the vml default, and any stray drawing parts leaves a clean,
standards-valid workbook with identical data.

Usage:  python3 strip_xlsx_drawings.py file1.xlsx [file2.xlsx ...]
Idempotent: a workbook with no orphan drawing refs is rewritten unchanged.
"""
import sys, re, os, zipfile, tempfile

REL_DRAW = re.compile(r'<Relationship[^>]*Type="[^"]*/(?:drawing|vmlDrawing)"[^>]*/>')
CT_OVER  = re.compile(r'<Override[^>]*PartName="/xl/drawings/[^"]*"[^>]*/>')
CT_VML   = re.compile(r'<Default[^>]*Extension="vml"[^>]*/>')
DIM      = re.compile(r'<dimension ref="[^"]*"/>')
CELL_REF = re.compile(r'<c r="([A-Z]+)(\d+)"')


def _col_to_num(col):
    n = 0
    for ch in col:
        n = n * 26 + (ord(ch) - 64)
    return n


def _num_to_col(n):
    s = ""
    while n > 0:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def fix_dimension(xml):
    # openxlsx writes <dimension ref="A1"/> regardless of the true used range, which makes
    # readers that trust the dimension tag (e.g. openpyxl read_only, pandas via it) report a
    # single 1x1 cell. Recompute the true range from the actual cell references so every reader
    # sees all columns/rows. Cell data itself is untouched.
    refs = CELL_REF.findall(xml)
    if not refs:
        return xml
    cols = [_col_to_num(c) for c, _ in refs]
    rows = [int(r) for _, r in refs]
    ref = f"{_num_to_col(min(cols))}{min(rows)}:{_num_to_col(max(cols))}{max(rows)}"
    return DIM.sub(f'<dimension ref="{ref}"/>', xml)

def clean(path):
    with zipfile.ZipFile(path) as z:
        items = [(i, z.read(i.filename)) for i in z.infolist()]
    changed = False
    out = []
    for info, data in items:
        name = info.filename
        # drop any real drawing/vml parts (defensive; usually none exist)
        if name.startswith("xl/drawings/"):
            changed = True
            continue
        if name.startswith("xl/worksheets/_rels/") and name.endswith(".rels"):
            txt = data.decode("utf-8")
            new = REL_DRAW.sub("", txt)
            if new != txt:
                changed = True; data = new.encode("utf-8")
        elif re.match(r"xl/worksheets/sheet\d+\.xml$", name):
            txt = data.decode("utf-8")
            new = fix_dimension(txt)
            if new != txt:
                changed = True; data = new.encode("utf-8")
        elif name == "[Content_Types].xml":
            txt = data.decode("utf-8")
            new = CT_VML.sub("", CT_OVER.sub("", txt))
            if new != txt:
                changed = True; data = new.encode("utf-8")
        out.append((info, data))
    if not changed:
        print(f"  {os.path.basename(path)}: already clean")
        return
    fd, tmp = tempfile.mkstemp(suffix=".xlsx", dir=os.path.dirname(os.path.abspath(path)))
    os.close(fd)
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for info, data in out:
            z.writestr(info, data)
    os.replace(tmp, path)
    print(f"  {os.path.basename(path)}: stripped orphan drawing refs")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: strip_xlsx_drawings.py file.xlsx [...]")
    for p in sys.argv[1:]:
        clean(p)
