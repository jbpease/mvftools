#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8

# Author: James B. Pease
"""Simple script to take argparse ArgumentParser object and make an .rst file
   This is a replacement for sphinx-argparse to implement custom
   formatting options.
"""


import sys
import os
import re
from argparse import _StoreAction, _HelpAction
from argparse import _StoreTrueAction, _StoreFalseAction
import importlib


def main(filepath, function_name, outfile=None):
    """Takes a filepath and function_name that returns an ArgumentParser and
       returns  text for an .rst file"""
    sys.path.append(os.path.abspath(os.path.dirname(filepath)))
    outfile = open(outfile, 'w') if outfile is not None else sys.stdout
    filename = os.path.basename(
        filepath[:-3] if filepath.endswith(".py") else filepath)
    x = importlib.import_module(filename)
    met = getattr(x, function_name)
    met.selfdoc = True
    for fxn in met.__dict__:
        if fxn[0] == '_' or fxn == "selfdoc":
            continue
        print("Processing", fxn)
        xmet = getattr(met, fxn)
        parser = xmet(met)
        description = xmet.__doc__
        description = re.sub(r'\n\s*', '', description)
        print(("---\n\n"
               "## {}\n"
               "***{}***\n\n"
               "**Parameters**\n\n").format(
                   fxn,
                   description,
               ), file=outfile)
        entries = []
        for y in parser._actions:
            # print(type(y).__name__, y)
            if (isinstance(y, _StoreTrueAction) is True) or (
                isinstance(y, _StoreFalseAction) is True) or (
                    isinstance(y, _HelpAction) is True) or (
                    type(y).__name__ == "_IntRangeAction"):
                optstring = '/'.join(y.option_strings)
                sortkey = y.option_strings[0].strip('-')
                if isinstance(y, _HelpAction) is True:
                    sortkey = '00000' + sortkey
                elif sortkey[0].islower():
                    sortkey = '0' + sortkey
                entries.append([
                    1 if isinstance(y, _HelpAction) else 2,
                    sortkey,
                    ("\n"
                     "``{}`` = {} (flag, default={})\n\n"
                     ).format(
                         optstring,
                         y.help,
                         y.default)
                    ])
            elif isinstance(y, _StoreAction) is True:
                if len(y.option_strings) == 0:
                    optstring = y.dest
                    sortkey = y.dest
                    priority = 0
                else:
                    optstring = "``{}``{}".format(
                        '/'.join(y.option_strings),
                        " (required)" if y.required is True else "")
                    sortkey = y.option_strings[0].strip('-')
                    priority = 1 if y.required is True else 2
                if sortkey[0].islower():
                    sortkey = '0' + sortkey
                ytype = str(y.type)
                if ytype == "<built-in function open>":
                    ytype = "file path"
                elif "function abspath" in ytype:
                    ytype = "file path"
                elif ytype == "<class 'float'>":
                    ytype = "float"
                elif ytype == "<class 'int'>":
                    ytype = "integer"
                elif ytype == "<class 'str'>":
                    ytype = "string"
                entries.append([
                    priority,
                    sortkey,
                    (""
                     "{} = {} (type={}, default={})\n"
                     "{}").format(
                          optstring,
                          y.help,
                          ytype,
                          y.default,
                          ("" if y.choices is None
                           else "Choices: {}\n".format(
                             y.choices)))
                ])
        for entry in sorted(entries):
            if "==SUPPRESS==" in entry[2]:
                continue
            print(entry[2], file=outfile)
    return ''


if __name__ == '__main__':
    # main(sys.argv[1], sys.argv[2])
    main("../mvftools", "MVFcall", 'prog-desc.md')
