import os
import time
import numpy as np
import pandas as pd

from ..ebsd.project import ScanData, _rename_columns

__all__ = ['load_ang_file', 'load_scandata']


def _parse_info_header(line, pattern, dtype=str):
    info = dtype(line.split(pattern)[-1].strip())
    print(line.strip())
    return info


def load_ang_file(fname):
    """
    Loads a TSL ang file. 
    The fields of each line in the body of the TSL ang file are as follows:

        phi1 Phi phi2 x y IQ CI ph intensity fit

    where:
    phi1, Phi, phi2 : Euler angles (in radians) in Bunge's notation for 
        describing the lattice orientations and are given in radians.
        A value of 4 is given to each Euler angle when an EBSP could 
        not be indexed. These points receive negative confidence index 
        values when read into an OIM dataset.
    x, y : The horizontal and vertical coordinates of the points in the
        scan, in micrometers. The origin (0,0) is defined as the top-left
        corner of the scan.
    IQ : The image quality parameter that characterizes the contrast of 
        the EBSP associated with each measurement point.
    CI : The confidence index that describes how confident the software is
        that it has correctly indexed the EBSP, i.e., confidence that the
        angles are correct.
    ph : The material phase identifier. This field is 0 for single phase 
        OIM scans or 1,2,3... for multi-phase scans.
    intensity : An integer describing the intensity from whichever detector
        was hooked up to the OIM system at the time of data collection, 
        typically a forward scatter detector.
    fit : The fit metric that describes how the indexing solution matches 
        the bands detected by the Hough transform or manually by the user.
    footers:
        In addition there also may be extra sections - such as EDS counts 
        data.

    Parameters
    ----------
    fname : string
        Path to the ang file

    Returns
    -------
    scan : ScanData object

    """
    t0 = time.time()
    print('Reading file \"{}\"...'.format(fname))

    grid = None
    dx = None
    dy = None
    ncols_odd = None
    ncols_even = None
    nrows = None
    # Read and parse header
    with open(fname) as f:
        header = []
        for line in f:
            # If header
            if line[0] == '#' or line[0] == '\n':
                header.append(line)

                pattern = '# GRID:'
                if pattern in line:
                    grid = _parse_info_header(line, pattern)
                    continue
                pattern = '# XSTEP:'
                if pattern in line:
                    dx = _parse_info_header(line, pattern, float)
                    continue
                pattern = '# YSTEP:'
                if pattern in line:
                    dy = _parse_info_header(line, pattern, float)
                    continue
                pattern = '# NCOLS_ODD:'
                if pattern in line:
                    ncols_odd = _parse_info_header(line, pattern, int)
                    continue
                pattern = '# NCOLS_EVEN:'
                if pattern in line:
                    ncols_even = _parse_info_header(line, pattern, int)
                    continue
                pattern = '# NROWS:'
                if pattern in line:
                    nrows = _parse_info_header(line, pattern, int)
                    continue
            else:
                break

    if grid is None:
        raise Exception('Missing grid info')

    # Uses pandas to read ang file. pd.read_csv returns a pandas DataFrame
    data = pd.read_csv(fname, header=None, comment='#', delim_whitespace=True)

    # Rename the columns
    columns = list(data.columns)
    columns[:10] = ['phi1', 'Phi', 'phi2', 'x',
                    'y', 'IQ', 'CI', 'ph', 'intensity', 'fit']
    data.columns = columns

    print('\n{} points read in {:.2f} s'.format(len(data), time.time() - t0))

    return ScanData(data, grid, dx, dy, ncols_odd, ncols_even, nrows, header)

def load_ctf_file(fname):
    """
    Loads a Oxford ctf file. 
    The fields of each line in the body of the TSL ang file are as follows:

        phi1 Phi phi2 x y IQ CI ph intensity fit

    where:
    phi1, Phi, phi2 : Euler angles (in radians) in Bunge's notation for 
        describing the lattice orientations and are given in radians.
        A value of 4 is given to each Euler angle when an EBSP could 
        not be indexed. These points receive negative confidence index 
        values when read into an OIM dataset.
    x, y : The horizontal and vertical coordinates of the points in the
        scan, in micrometers. The origin (0,0) is defined as the top-left
        corner of the scan.
    IQ : The image quality parameter that characterizes the contrast of 
        the EBSP associated with each measurement point.
    CI : The confidence index that describes how confident the software is
        that it has correctly indexed the EBSP, i.e., confidence that the
        angles are correct.
    ph : The material phase identifier. This field is 0 for single phase 
        OIM scans or 1,2,3... for multi-phase scans.
    intensity : An integer describing the intensity from whichever detector
        was hooked up to the OIM system at the time of data collection, 
        typically a forward scatter detector.
    fit : The fit metric that describes how the indexing solution matches 
        the bands detected by the Hough transform or manually by the user.
    footers:
        In addition there also may be extra sections - such as EDS counts 
        data.

    Parameters
    ----------
    fname : string
        Path to the ang file

    Returns
    -------
    scan : ScanData object

    """
    t0 = time.time()
    print('Reading file \"{}\"...'.format(fname))

    # Read and parse header
    with open(fname) as f:
        header = []
        nmatches = 0
        for line in f:
            
            # If header
            if 'MAD' in line:
                break

            else:
                header.append(line)

                pattern = 'JobMode'
                if pattern in line:
                    grid = 'SqrGrid'
                    nmatches += 1
                    continue
                pattern = 'XStep'
                if pattern in line:
                    dx = _parse_info_header(line, pattern, float)
                    nmatches += 1
                    continue
                pattern = 'YStep'
                if pattern in line:
                    dy = _parse_info_header(line, pattern, float)
                    nmatches += 1
                    continue
                pattern = 'XCells'
                if pattern in line:
                    ncols_odd = _parse_info_header(line, pattern, int)
                    nmatches += 1
                    continue
                pattern = 'YCells'
                if pattern in line:
                    ncols_even = _parse_info_header(line, pattern, int)
                    nmatches += 1
                    continue

    if nmatches != 5:
        raise Exception('Info about scandata is missing in the file header.')

    # Uses pandas to read ang file. pd.read_csv returns a pandas DataFrame
    data = pd.read_csv(fname, header=None, comment='#',
                       delim_whitespace=True, skiprows=(len(header) + 1))
    # Rename the columns
    columns = list(data.columns)
    
    columns[:11] = ['ph', 'x', 'y', 'bands',
                    'error', 'phi1', 'Phi', 'phi2', 'fit', 'IQ', 'intensity']
    data.columns = columns

    data.phi1 = np.deg2rad(data.phi1)
    data.Phi = np.deg2rad(data.Phi)
    data.phi2 = np.deg2rad(data.phi2)

    print('\n{} points read in {:.2f} s'.format(len(data), time.time() - t0))

    return ScanData(data, grid, dx, dy, ncols_odd, ncols_odd, ncols_even, header)


def load_scandata(fname):
    """
    Load EBSD scan data

    Parameters
    ----------
    fname : string
        Path to the file to be loaded

    Returns
    -------
    scan : ScanData object

    """
    ext = os.path.splitext(fname)[-1]

    if ext == '.ang':
        scan = load_ang_file(fname)
    elif ext == '.ctf':
        scan = load_ctf_file(fname)
    else:
        raise Exception('File extension "{}" not supported'.format(ext))

    return scan
