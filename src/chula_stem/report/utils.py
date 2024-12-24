#!/usr/bin/env ipython

def style_cells(
    start: tuple,
    ncols: int = 0,
    nrows: int = 0,
    fontname: str = "",
    fontsize: str = "",
    align: str = "",
    background: str = "",
    valign: str = "",
    textcolor=None,
    underline: tuple = (),
) -> list:
    """
    Coordinates for table style are given as (column, row)
    """
    parameter_map: dict = {
        fontname: "FONTNAME",
        fontsize: "FONTSIZE",
        textcolor: "TEXTCOLOR",
        background: "BACKGROUND",
        align: "ALIGN",
        valign: "VALIGN",
        underline: "LINEBELOW",
    }
    if ncols and nrows:
        end: tuple = start[0] + ncols - 1, start[1] + nrows - 1
    else:
        end = (-1, -1)

    def style_helper(format: str, value) -> tuple:
        if isinstance(value, tuple):
            return (format, start, end) + value
        else:
            return (format, start, end, value)

    styles: list = []
    for param, name in parameter_map.items():
        if param:
            styles.append(style_helper(name, param))
    return styles
