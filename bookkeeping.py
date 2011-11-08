#!/usr/bin/env python

"""
reads and writes files, sets up logging
"""

import logging
import sys
import textwrap

class Log:

    def __init__(self, iter=None):

        if iter is None:
            self.iter = 1
        else:
            self.iter = iter

        self.logfile = logging.FileHandler(filename="genstruct.out")
        self.console = logging.StreamHandler(sys.stdout)
        self.init_logger()
        # set up writing to file and terminal

    def init_logger(self):
        # This is the logging handler
        self.logger = logging.getLogger('SBU building iteration %-6i'
                           %self.iter)

        # for now, keep as verbose
        stdout_level = logging.DEBUG
        file_level = logging.DEBUG

        self.long_format = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s: %(message)s',
                datefmt='%m/%d/%Y %I:%M:%S %p')
        self.time_format = logging.Formatter(
                '%(asctime)s - %(levelname)s: %(message)s',
                datefmt='%m/%d/%Y %I:%M:%S %p')
        self.reduced_format = logging.Formatter(
                '%(levelname)s: %(message)s')
        
        self.message = logging.Formatter(
                '%(message)s')
        
        # set up verbosity level
        self.console.setLevel(stdout_level)
        self.logfile.setLevel(file_level)
        self.logger.setLevel(stdout_level)

        # initialize the format as time_format 
        self.console.setFormatter(self.time_format)
        self.logfile.setFormatter(self.time_format)

        # add both the stdout and log file to the logger
        self.logger.addHandler(self.console)
        self.logger.addHandler(self.logfile)

    def debug(self, msg, format = None):
        """Send DEBUG to the logging handlers"""
        msg = textwrap.wrap(msg)
        if format is not None:
            self.console.setFormatter(format)
            self.logfile.setFormatter(format)
        for line in msg:
            self.logger.debug(line)

    def info(self, msg, format = None):
        """Send INFO to the logging handlers"""
        msg = textwrap.wrap(msg)
        if format is not None:
            self.console.setFormatter(format)
            self.logfile.setFormatter(format)
        for line in msg:
            self.logger.info(line)

    def warn(self, msg, format = None):
        """Send WARN to the logging handlers"""
        msg = textwrap.wrap(msg)
        if format is not None:
            self.console.setFormatter(format)
            self.logfile.setFormatter(format)
        for line in msg:
            self.logger.warning(line)

    def error(self, msg, format = None):
        """Send ERROR to the logging handlers"""
        msg = textwrap.wrap(msg)
        if format is not None:
            self.console.setFormatter(format)
            self.logfile.setFormatter(format)
        for line in msg:
            self.logger.error(line)

        sys.exit(0)

