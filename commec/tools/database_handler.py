#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Abstract base class for a database handler. 
Customize to screen desired database functionality.
"""
import os
import glob
from dataclasses import dataclass
import subprocess
import logging

@dataclass
class DatabaseVersion():
    """ Container class for outputting version related information from a database."""
    def __init__(self, input_version : str = "x.x.x", input_date : str = "Null", input_comment :str = ""):
        self.version_string = input_version
        self.version_date = input_date
        self.additional_comment = input_comment
    version_string : str
    version_date : str
    additional_comment : str

class DatabaseHandler():
    """ 
    Abstract Class for holding the directory, and file of a database, 
    including input and output files for screening, as well as screening arguments. 
    Override the Screen() method for a custom database implementation.
    """
    # Start Database Handler API
    def screen(self):
        """ Virtual function to be called by any child database implementation"""
        raise NotImplementedError(
            "This class must override the .screen() to use DatabaseHandler."
        )

    def get_version_information(self) -> DatabaseVersion:
        """ Override to call version retrieval information on a database."""
        new_version_info = DatabaseVersion("", "", "Version retrieval not yet supported for this database")
        return new_version_info

    def check_output(self):
        """ 
        Simply checks for the existance of the output file, 
        indicating that the database screening ran. Can be overrided if
        more complex checks for a particular db is desired.
        """
        if os.path.isfile(self.out_file): # Consider adding checks from File_Tools such that empty outputs are invalid?
            return True
        return False

    def __init__(self, directory : str, database_file : str, input_file : str, out_file : str):
        self.db_directory = directory
        self.db_file = database_file
        self.out_file = out_file
        self.input_file = input_file
        self.temp_log_file = f"{self.out_file}.log.tmp"
        self.arguments_dictionary = {}
        self.validate_directory()

    # End override API

    def check_input(self):
        """ 
        Simply checks for the existance of the input file, This is separated 
        from database location, as sometimes the input file is generated at a 
        previous step, and may not exist during DatabaseHandler instantiation.
        """
        if os.path.isfile(self.input_file):
            return True
        return False
    
    def is_empty(self) -> bool:
        """
        is_empty

        usage: check that a file is empty
        input:
        - name of file
        """
        try:
            abspath = os.path.abspath(os.path.expanduser(self.out_file))
            return os.path.getsize(abspath) == 0
        except OSError:
            # If there is an error (including FileNotFoundError) consider it empty
            return True
    
    def has_hits(self, filepath: str = None) -> bool:
        """
        has_hits
        usage: check to see if the file contains any hits (lines that don't start with #)
        Override for custom hit behaviour for a database, 
        currently returns true for any file with lines not starting with a comment #
        input:
        - path to file
        """
        file_to_check = filepath
        if file_to_check is None:
            file_to_check = self.input_file

        try:
            with open(file_to_check, "r", encoding="utf-8") as file:
                for line in file:
                    # Strip leading and trailing whitespace and check the first character
                    if not line.strip().startswith("#"):
                        # Found a hit!
                        return True
            return False
        except FileNotFoundError:
            # The file does not exist
            return False

    def get_arguments(self) -> list:
        """ 
        convert the arguments dictionary into a list, 
        structurally ready for appending to a command list of strs.
        """
        my_list = []
        for key, value in self.arguments_dictionary.items():
            my_list.append(str(key))
            if isinstance(value, list):
                my_list.append(" ".join(value))  # Extend the list with all elements in the array
            else:
                my_list.append(str(value))  # Append the value directly if it's not a list
        return my_list

    def validate_directory(self):
        """ 
        Validates that the directory, 
        and database file exists. Called on init.
        """
        if not os.path.isdir(self.db_directory):
            raise FileNotFoundError(f"Mandatory screening directory {self.db_directory} not found.")
        if not os.path.isfile(self.db_file):

            # The specified file doesn't exist, try to find a pattern of multiple files
            filename, extension = os.path.splitext(self.db_file)
            search_file = os.path.join(self.db_directory, "*" + os.path.basename(filename) + "*" + extension)
            files = glob.glob(search_file)
            if len(files) == 0:
                logging.error(f"{self.__class__.__name__} : Bad database file!")
                raise FileNotFoundError(f"Mandatory screening directory {self.db_file} not found.")

    def run_as_subprocess(self, command, out_file, raise_errors=False):
        """
        Run a command using subprocess.run, piping stdout and stderr to `out_file`.
        """
        logging.debug("SUBPROCESS: " + " ".join(command))

        with open(out_file, "a", encoding="utf-8") as f:

            result = subprocess.run(
                command, stdout=f, stderr=subprocess.STDOUT, check=raise_errors
            )

            if not self.is_succesful_result(result):
                command_str = ' '.join(command)
                logging.info("\t ERROR: command %s failed", command_str)
                raise RuntimeError(
                    f"subprocess.run of command '{command_str}' encountered error."
                    f" Check {out_file} for logs."
                )
            
    def is_succesful_result(self , result : subprocess.CompletedProcess[bytes]):
        """ Override for custom return code behaviour"""
        if result.returncode != 0:
            return False
        return True
    
    def __del__(self):
        if os.path.exists(self.temp_log_file):
            os.remove(self.temp_log_file)

