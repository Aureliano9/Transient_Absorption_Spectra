# Copyright 2020 Patrick C. Tapping
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Support for reading and writing Ultrafast Systems .ufs files.
"""

import struct

import numpy as np

def _read_string(ufsdata, cursor):
    """
    Read a string from ``ufsdata`` at position given by the cursor.
    
    An unsigned 32 bit integer gives length of string in bytes, then string data follows.
    Returns the string and the incremented cursor value as a tuple.
    """
    string_length = struct.unpack_from('>I', ufsdata, offset=cursor)
    cursor += 4
    string = struct.unpack_from('>{}s'.format(string_length[0]), ufsdata, offset=cursor)
    cursor += string_length[0]
    return (string[0].decode('utf8'), cursor)

def _read_uint(ufsdata, cursor):
    """
    Read an unsigned 32 bit integer from ``ufsdata`` at position given by the cursor.
    
    Returns the number and the incremented cursor value as a tuple.
    """
    number = struct.unpack_from('>I', ufsdata, offset=cursor)
    cursor += 4
    return (number[0], cursor)
    
def _read_double(ufsdata, cursor):
    """
    Read an 64 bit floating point number from ``ufsdata`` at position given by the cursor.
    
    Returns the number and the incremented cursor value as a tuple.
    """
    number = struct.unpack_from('>d', ufsdata, offset=cursor)
    cursor += 8
    return (number[0], cursor)

def _read_doubles(ufsdata, cursor, count):
    """
    Read a series of 64 bit floating point numbers from ``ufsdata`` at position given by the cursor.
    
    Returns the list of numbers and the incremented cursor value as a tuple.
    """
    numbers = struct.unpack_from('>{}d'.format(count), ufsdata, offset=cursor)
    cursor += 8*count
    return (np.array(numbers), cursor)

def _write_string(ufsstream, string):
    """
    Write a string to ``ufsstream`` following the UFS formatting protocol.
    """
    ufsstream.write(struct.pack('>I', len(string)))
    ufsstream.write(struct.pack('>{}s'.format(len(string)), string))

def _write_uint(ufsstream, value):
    """
    Write an unsigned integer to ``ufsstream`` following the UFS formatting protocol.
    """
    ufsstream.write(struct.pack('>I', value))

def _write_double(ufsstream, value):
    """
    Write a double-precision floating point number to ``ufsstream`` following the UFS formatting protocol.
    """
    ufsstream.write(struct.pack('>d', value))

def _write_doubles(ufsstream, values):
    """
    Write a list of double-precision floating point number to ``ufsstream`` following the UFS formatting protocol.
    """
    ufsstream.write(struct.pack('>{}d'.format(len(values)), *values))

class UFSData():
    def __init__(self, ufsdata=None, ufsfile=None):
        """
        Read the byte stream given by ``ufsdata`` and initialise a UFSData instance.

        :param ufsdata: Bytes-like object from which to decode UFS data.
        """
        #: Version string for the UFS file format.
        self.version = None
        #: Label string for the first axis.
        self.axis1_label = None
        #: Unit string for the first axis.
        self.axis1_units = None
        #: Values for the data points on the first axis.
        self.axis1_data = None
        #: Label string for the second axis.
        self.axis2_label = None
        #: Unit string for the second axis.
        self.axis2_units = None
        #: Values for the data points on the second axis.
        self.axis2_data = None
        #: Unit string for the matrix data points.
        self.data_units = None
        #: Data matrix as a 2D array.
        self.data_matrix = None
        #: String containing metadata information.
        self.metadata = None
        if not ufsdata is None:
            self.read_data(ufsdata)
        elif not ufsfile is None:
            self.read_file(ufsfile)
    
def read_data(self, ufsdata):
        """
        Read the byte stream given by ``ufsdata`` and store decoded data.

        :param ufsdata: Bytes-like object from which to decode UFS data.
        """

        # Keep track of our current location through the file
        cursor = 0x0

        (self.version, cursor) = _read_string(ufsdata, cursor)

        (self.axis1_label, cursor) = _read_string(ufsdata, cursor)
        (self.axis1_units, cursor) = _read_string(ufsdata, cursor)
        (axis1_count, cursor) = _read_uint(ufsdata, cursor)
        (self.axis1_data, cursor) = _read_doubles(ufsdata, cursor, axis1_count)

        (self.axis2_label, cursor) = _read_string(ufsdata, cursor)
        (self.axis2_units, cursor) = _read_string(ufsdata, cursor)
        (axis2_count, cursor) = _read_uint(ufsdata, cursor)
        (self.axis2_data, cursor) = _read_doubles(ufsdata, cursor, axis2_count)

        (self.data_units, cursor) = _read_string(ufsdata, cursor)
        (data_size0, cursor) = _read_uint(ufsdata, cursor)
        # Don't actually know how they handle 3D data... might be data_size0 += 1
        if data_size0 == 0: data_size0 = 1
        (data_size1, cursor) = _read_uint(ufsdata, cursor)
        (data_size2, cursor) = _read_uint(ufsdata, cursor)

        (self.data_matrix, cursor) = _read_doubles(ufsdata, cursor, (data_size0*data_size1*data_size2))
        self.data_matrix =self.data_matrix.reshape((data_size0, data_size1, data_size2))
        
        (self.metadata, cursor) = _read_string(ufsdata, cursor)


        # print("Version string: " + self.version)
        # print("Axis1: {}, {} values, {} to {} {}".format(self.axis1_label, axis1_count, self.axis1_data[0], self.axis1_data[-1], self.axis1_units))
        # print("Axis2: {}, {} values, {} to {} {}".format(self.axis2_label, axis2_count, self.axis2_data[0], self.axis2_data[-1], self.axis2_units))
        # print("Data: {} {} x {}".format(self.data_units, data_size1, data_size2))
        # print(self.metadata)
    
def read_file(self, filename):
        """
        Read a ``.ufs`` file.

        :param filename: Filename to read.
        """
        with open(filename, "rb") as f:
            filedata = f.read()
            self.read_data(filedata)

    
def write_data(self, destination):
        """
        Write the raw UFS data out to a stream.

        :param destination: Stream to write UFS data to.
        """
        _write_string(destination, self.version)
            
        _write_string(destination, self.axis1_label)
        _write_string(destination, self.axis1_units)
        _write_uint(destination, self.axis1_data.shape[0])
        _write_doubles(destination, self.axis1_data)
        
        _write_string(destination, self.axis2_label)
        _write_string(destination, self.axis2_units)
        _write_uint(destination, self.axis2_data.shape[0])
        _write_doubles(destination, self.axis2_data)
        
        _write_string(destination, self.data_units)
        _write_uint(destination, self.data_matrix.shape[0] - 1) # ??? or 1 if zero?
        _write_uint(destination, self.data_matrix.shape[1])
        _write_uint(destination, self.data_matrix.shape[2])
        _write_doubles(destination, self.data_matrix.flat)
        
        _write_string(destination, self.metadata)


def write_file(self, filename):
        """
        Write the raw UFS data out to a file.

        :param filename: Filename to write UFS data to.
        """
        with open(filename, "wb") as f:
            self.write_data(f)

obj = UFSData()
# obj.read_file('continuum spectrum_channel1 raw Pump off.ufs')