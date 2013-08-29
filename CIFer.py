#!/usr/bin/env python

class CIF(object):

    def __init__(self):
        self._data = {}

    def add_data(**kwargs):
        for key, val in kwargs.items():
            self._data[key].setdefault([])
            self._data.append(val)
