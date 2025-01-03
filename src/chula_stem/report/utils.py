#!/usr/bin/env ipython


from reportlab.lib.colors import HexColor


def style_cells(
    start: tuple,
    ncols: int = 0,
    nrows: int = 0,
    fontname: str = "",
    fontsize: str = "",
    align: str = "",
    background: str = "",
    end: tuple = None,
    valign: str = "",
    textcolor=None,
    underline: tuple = (),
    lineafter: tuple = (),
    box: tuple = (),
    grid: tuple = (),
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
        lineafter: "LINEAFTER",
        valign: "VALIGN",
        box: "BOX",
        grid: "GRID",
        underline: "LINEBELOW",
    }
    if ncols and nrows and (not end):
        end_actual: tuple = start[0] + ncols - 1, start[1] + nrows - 1
    elif not end:
        end_actual = (-1, -1)
    else:
        end_actual = end

    def style_helper(format: str, value) -> tuple:
        if isinstance(value, tuple):
            return (format, start, end_actual) + value
        else:
            return (format, start, end_actual, value)

    styles: list = []
    for param, name in parameter_map.items():
        if param:
            styles.append(style_helper(name, param))
    return styles


def alternating_bg(ncols: int, c1: str, c2: str, offset=1, **kwargs) -> list:
    """Produce a series of alternating background colors for table cells

    All cells in the 1st column will have bg "c1", cells in the 2nd column bg "c2",
    which repeats
    """
    style_list = []
    for i in range(ncols):
        color = c1 if i % 2 == 0 else c2
        style_list.extend(
            style_cells(
                start=(i, offset), end=(i, -1), background=HexColor(color), **kwargs
            )
        )
    return style_list
