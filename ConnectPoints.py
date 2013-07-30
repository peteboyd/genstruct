#!/usr/bin/env python
import numpy as np

class ConnectPoint(object):
    
    def __init__(self):
        """Origin describes the point of intersection of two parameters,
        z describes the vector pointing along the bond (parallel),
        y describes a vector perpendicular to z for alignment purposes.
        """
        null = np.array([0., 0., 0., 1.])
        self.identifier = None
        self.origin = null.copy()
        self.y = null.copy()
        self.z = null.copy()
        # flag to determine if the point has been attached
        self.connected = False
        # tuple (sbu_order, connect.identifier)
        self.sbu_bond = None
        
    def from_config(self, line):
        """ Obtain the connectivity information from the config .ini file."""
        line = line.strip().split()
        self.identifier = int(line[0])
        # obtain the coordinate information.
        self.origin[:3] = np.array([float(x) for x in line[1:4]])
        self.z[:3] = np.array([float(x) for x in line[4:7]])
        self.y[:3] = np.array([float(x) for x in line[7:10]])
        try:
            self.special = int(line[11])
        except IndexError:
            self.special = None
        except ValueError:
            self.special = None
    
    def rotate(self, R):
        self.origin = np.dot(R, self.origin)
        self.y[:3] = np.dot(R[:3,:3], self.y[:3])
        self.z[:3] = np.dot(R[:3,:3], self.z[:3])
    
    def translate(self, vector):
        self.origin[:3] += vector
    
    def __neg__(self):
        self.z[:3] = -self.z[:3]
        return self
        
    @property
    def normal(self):
        try:
            return self._normal
        except AttributeError:
            normal = np.cross(self.y[:3], self.z[:3])
            self._normal = normal / np.linalg.norm(normal)
            return self._normal
    
        
        