#!/usr/bin/env python
_version_=3.0
from optparse import OptionParser
import optparse
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import os
try:
        from StringIO import StringIO ## for Python 2
except ImportError:
        from io import StringIO ## for Python 3
import sys
import re
import copy
from ast import literal_eval
from logging import info, debug, warning, error, critical

class Options(object):
    version = 3.001
    def __init__(self):

        self._command_options()
        self.job = configparser.SafeConfigParser()
        self._set_paths()
        self._load_defaults()
        self._load_job()
        self._set_attr()
    
    def _set_paths(self):
        if __name__ != '__main__':
            self.script_dir = os.path.dirname(__file__)
        else:
            self.script_dir = os.path.abspath(sys.path[0])
        self.job_dir = os.getcwd()
        self.jobname = os.path.splitext(os.path.basename(
                                       self.input_file))[0]

    def _command_options(self):
        """Load data from the command line."""

        usage = "%prog [options] input_file"
        version = "%prog " + "%f"%(_version_)
        parser = OptionParser(usage=usage, version=version)
        group = optparse.OptionGroup(parser, "Verbosity Options")
        group.add_option("-s", "--silent", action="store_true",
                          dest="silent",
                          help="Print nothing to the console.")
        group.add_option("-q", "--quiet", action="store_true",
                          dest="quiet",
                          help="Print only warnings and errors.")
        group.add_option("-v", "--verbose", action="store_true",
                          dest="verbose",
                          help="Print everything to the console.")
        parser.add_option_group(group)
        (self.cmd_options, local_args) = parser.parse_args()

        if not local_args:
            parser.print_help()
            sys.exit(1)
        elif len(local_args) != 1:
            error("Only one argument required, the input file")
            sys.exit(1)
        else:
            self.input_file = os.path.abspath(local_args[0])

        
    def _load_defaults(self):
        default_path = os.path.join(self.script_dir, 'defaults.ini')        
        try:
            filetemp = open(default_path, 'r')
            default = filetemp.read()
            filetemp.close()
            if not '[defaults]' in default.lower():
                default = '[defaults]\n' + default
            default = StringIO(default)
        except IOError:
            error("Error loading defaults.ini")
            default = StringIO('[defaults]\n')
        self.job.readfp(default)
        
    def _load_job(self):
        """Load data from the local job name."""
        if self.input_file is not None:
            try:
                filetemp = open(self.input_file, 'r')
                job = filetemp.read()
                filetemp.close()
                if not '[job]' in job.lower():
                    job = '[job]\n' + job
                job = StringIO(job)
            except IOError:
                job = StringIO('[job]\n')
        else:
            job = StringIO('[job]\n')
        self.job.readfp(job)
        
    def _set_attr(self):
        """Sets attributes to the base class. default options are over-written
        by job-specific options.

        """
        for key, value in self.job.items('defaults'):
            value = self.get_val('defaults', key)
            setattr(self, key, value)
        for key, value in self.job.items('job'):
            value = self.get_val('job', key)
            setattr(self, key, value)
        for key, value in self.cmd_options.__dict__.items():
            setattr(self, key, value)

    def get_val(self, section, key):
        """Returns the proper type based on the key used."""
        # known booleans
        booleans = ['verbose', 'quiet', 'silent', 'exhaustive',
                    'gen_single_bonding_set', 'calc_sbu_surface_area',
                    'calc_max_sbu_span', 'create_sbu_input_files', 'find_symmetry',
                    'relaxed_topology', 'build_from_tree']
        floats = ['distance_tolerance', 'overlap_tolerance',
                  'bond_angle_tolerance',
                  'cell_angle_cutoff', 'symmetry_precision']
        integers = ['organic_sbu_per_structure',
                    'metal_sbu_per_structure',
                    'structure_sbu_length',
                    'max_trials',
                    'max_structures']
        lists = ['topologies', 'sbu_files', 'organic_sbus',
                 'metal_sbus']
        tuple_of_tuples = ['sbu_combinations']
        
        if key in booleans:
            try:
                val = self.job.getboolean(section, key)
            except ValueError:
                val = False
        # known integers
        elif key in integers:
            try:
                val = self.job.getint(section, key)
            except ValueError:
                val = 0;
        # known floats
        elif key in floats:
            try:
                val = self.job.getfloat(section, key)
            except ValueError:
                val = 0.
        # known lists
        elif key in lists:
            p = re.compile('[,;\s]+')
            val = p.split(self.job.get(section, key))
            try:
                val = [int(i) for i in val if i]
            except ValueError:
                val = [i for i in val if i]
                
        # tuple of tuples.
        elif key in tuple_of_tuples:
            val = literal_eval(self.job.get(section, key)) \
                    if self.job.get(section, key) else None
            # failsafe if only one tuple is presented, need to embed it.
            if val is not None:
                if isinstance(val[0], int) or isinstance(val[0], float):
                    val = [val]
        else:
            val = self.job.get(section, key)
        return val

def Terminate():
    info("Genstruct terminated normally")
    sys.exit(1)
