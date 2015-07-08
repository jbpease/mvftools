#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools

MVF Analysis Base
@author: James B. Pease
@author: Ben K. Rosenzweig


This file is part of MVFtools.

MVFtools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MVFtools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MVFtools.  If not, see <http://www.gnu.org/licenses/>.
"""

import os


class Counter(dict):
    """dict subclass for counter with integer values"""

    def add(self, key, val=1):
        """Adds integer value to key, default = 1)"""
        self[key] = self.get(key, 0) + val

    def subtract(self, key, val=1):
        """Subtracts integer value to key, default = 1)"""
        self[key] = self.get(key, 0) - val
        return ''

    def iter_sorted(self, percentage=True, n_values=None, sort='desc'):
        """Iterates values sorting in a list"""
        total = float(sum(self.values()))
        if not total:
            yield "Nothing counted!"
        entries = sorted([(v, k) for (k, v) in self.items()],
                         reverse=(sort != 'asc'))
        for (val, k) in entries[:n_values]:
            if percentage:
                yield "{!s} {!s} {!s}%\n".format(
                    k, val, round(float(val) * 100/total, 2))
            else:
                yield "{!s} {!s}\n".format(k, val)


class OutputFile(object):
    """Set up Output File
        Params:
            path: file path
            headers: list of header elements
    """

    def __init__(self, path, headers):
        self.headers = headers
        self.path = os.path.abspath(path)
        self.write_headers()

    def write_headers(self):
        """Write headers to file"""
        with open(self.path, 'w') as outfile:
            outfile.write('#' + '\t'.join(self.headers) + "\n")
        return ''

    def write_entry(self, entry):
        """Writes entry to file
            Arguments:
                entry: dict of values with keys matching header
        """
        with open(self.path, 'a') as outfile:
            outfile.write("\t".join([str(k not in entry and '.' or entry[k])
                                     for k in self.headers]) + "\n")
        return ''


class AnalysisModule(object):
    """General Functions for Analysis Modules
        Params:
            data = data dict
            params = parameters dict
            **kwargs = all other keywords are added to params as key=value

    """
    def __init__(self, data=None, params=None, **kwargs):
        self.data = data or {}
        self.params = params or {}
        self.params.update([(k, v) for k, v in kwargs.items()
                            if k not in ['data', 'params']])
