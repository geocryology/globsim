import numpy as np
import re

try:
    import ESMF

    # Check ESMF version.  7.0.1 behaves differently than 7.1.0r
    ESMFv = int(re.sub("[^0-9]", "", ESMF.__version__))
    ESMFnew = ESMFv > 701

except ModuleNotFoundError:
        print("*** ESMF not imported, trying esmpy. ***")
        try:
            import esmpy as ESMF
        except ImportError:
            print('Could not import ESMF or esmpy')
            pass

from globsim.boundingbox import BoundingBox


def grid_create_from_coordinates(xcoords, ycoords,
                                 xcorners=False, ycorners=False, corners=False,
                                 domask=False, doarea=False, ctk=ESMF.TypeKind.R8):
    """
    Create a 2 dimensional Grid using the bounds of the x and y coordiantes.
    :param xcoords: The 1st dimension or 'x' coordinates at cell centers, as
        a Python list or numpy Array
    :param ycoords: The 2nd dimension or 'y' coordinates at cell centers, as
        a Python list or numpy Array
    :param xcorners: The 1st dimension or 'x' coordinates at cell corners,
        as a Python list or numpy Array
    :param ycorners: The 2nd dimension or 'y' coordinates at cell corners,
        as a Python list or numpy Array
    :param domask: boolean to determine whether to set an arbitrary mask or
        not
    :param doarea: boolean to determine whether to set an arbitrary area
        values or not
    :param ctk: the coordinate typekind
    :return: grid
    """
    [x, y] = [0, 1]
    
    # create a grid given the number of grid cells in each dimension, the center stagger location is allocated, the
    # Cartesian coordinate system and type of the coordinates are specified
    max_index = np.array([len(xcoords), len(ycoords)])
    grid = ESMF.Grid(max_index, staggerloc=[ESMF.StaggerLoc.CENTER], coord_sys=ESMF.CoordSys.CART, coord_typekind=ctk)
    
    # set the grid coordinates using numpy arrays, parallel case is handled using grid bounds
    gridXCenter = grid.get_coords(x)
    x_par = xcoords[grid.lower_bounds[ESMF.StaggerLoc.CENTER][x]:grid.upper_bounds[ESMF.StaggerLoc.CENTER][x]]
    gridXCenter[...] = x_par.reshape((x_par.size, 1))
    gridYCenter = grid.get_coords(y)
    y_par = ycoords[grid.lower_bounds[ESMF.StaggerLoc.CENTER][y]:grid.upper_bounds[ESMF.StaggerLoc.CENTER][y]]
    gridYCenter[...] = y_par.reshape((1, y_par.size))
    
    # create grid corners in a slightly different manner to account for the bounds format common in CF-like files
    if corners:
        grid.add_coords([ESMF.StaggerLoc.CORNER])
        lbx = grid.lower_bounds[ESMF.StaggerLoc.CORNER][x]
        ubx = grid.upper_bounds[ESMF.StaggerLoc.CORNER][x]
        lby = grid.lower_bounds[ESMF.StaggerLoc.CORNER][y]
        uby = grid.upper_bounds[ESMF.StaggerLoc.CORNER][y]
        gridXCorner = grid.get_coords(x, staggerloc=ESMF.StaggerLoc.CORNER)
        
        for i0 in range(ubx - lbx - 1):
            gridXCorner[i0, :] = xcorners[i0+lbx, 0]
            gridXCorner[i0 + 1, :] = xcorners[i0+lbx, 1]
            gridYCorner = grid.get_coords(y, staggerloc=ESMF.StaggerLoc.CORNER)
        
        for i1 in range(uby - lby - 1):
            gridYCorner[:, i1] = ycorners[i1+lby, 0]
            gridYCorner[:, i1 + 1] = ycorners[i1+lby, 1]
    
    # add an arbitrary mask
    if domask:
        mask = grid.add_item(ESMF.GridItem.MASK)
        mask[:] = 1
        mask[np.where((1.75 <= gridXCenter.any() < 2.25) &
        (1.75 <= gridYCenter.any() < 2.25))] = 0
    
    # add arbitrary areas values
    if doarea:
        area = grid.add_item(ESMF.GridItem.AREA)
        area[:] = 5.0
    
    return grid


def grid_create_from_coordinates_periodic(longitudes, latitudes, 
                                          lon_corners=False, lat_corners=False, corners=False,
                                          domask=False):
    """
    Create a 2 dimensional periodic Grid using the 'longitudes' and
        'latitudes'.
    :param longitudes: longitude coordinate values at cell centers
    :param latitudes: latitude coordinate values at cell centers
    :param lon_corners: longitude coordinate values at cell corners
    :param lat_corners: latitude coordinate values at cell corners
    :param corners: boolean to determine whether or not to add corner
        coordinates to this grid
    :param domask: boolean to determine whether to set an arbitrary mask or
        not
    :return: grid
    """
    [lon, lat] = [0, 1]
    
    # create a grid given the number of grid cells in each dimension the center stagger location is allocated
    max_index = np.array([len(longitudes), len(latitudes)])
    grid = ESMF.Grid(max_index, num_peri_dims=1, staggerloc=[ESMF.StaggerLoc.CENTER])
    
    # set the grid coordinates using numpy arrays, parallel case is handled using grid bounds
    gridXCenter = grid.get_coords(lon)
    lon_par = longitudes[grid.lower_bounds[ESMF.StaggerLoc.CENTER][lon]:grid.upper_bounds[ESMF.StaggerLoc.CENTER][lon]]
    gridXCenter[...] = lon_par.reshape((lon_par.size, 1))
    gridYCenter = grid.get_coords(lat)
    lat_par = latitudes[grid.lower_bounds[ESMF.StaggerLoc.CENTER][lat]:grid.upper_bounds[ESMF.StaggerLoc.CENTER][lat]]
    gridYCenter[...] = lat_par.reshape((1, lat_par.size))
    
    # create grid corners in a slightly different manner to account for the bounds format common in CF-like files
    if corners:
        grid.add_coords([ESMF.StaggerLoc.CORNER])
        lbx = grid.lower_bounds[ESMF.StaggerLoc.CORNER][lon]
        ubx = grid.upper_bounds[ESMF.StaggerLoc.CORNER][lon]
        lby = grid.lower_bounds[ESMF.StaggerLoc.CORNER][lat]
        uby = grid.upper_bounds[ESMF.StaggerLoc.CORNER][lat]
        gridXCorner = grid.get_coords(lon, staggerloc=ESMF.StaggerLoc.CORNER)
        
        for i0 in range(ubx - lbx - 1):
            gridXCorner[i0, :] = lon_corners[i0+lbx, 0]
            gridXCorner[i0 + 1, :] = lon_corners[i0+lbx, 1]
            gridYCorner = grid.get_coords(lat, staggerloc=ESMF.StaggerLoc.CORNER)
        
        for i1 in range(uby - lby - 1):
            gridYCorner[:, i1] = lat_corners[i1+lby, 0]
            gridYCorner[:, i1 + 1] = lat_corners[i1+lby, 1]
    
    # add an arbitrary mask
    if domask:
        mask = grid.add_item(ESMF.GridItem.MASK)
        mask[:] = 1
        mask[np.where((1.75 <= gridXCenter.any() < 2.25) &
        (1.75 <= gridYCenter.any() < 2.25))] = 0
    
    return grid


def get_buffered_slices(grid: "ESMF.Grid",
                        lon_indices:'np.ndarray',
                        lat_indices:'np.ndarray',
                        b=3) -> "tuple[slice,slice]":
    """ get  longitude, latitude slices for 'grid' that cover '*_indices' plus a buffer of 'b' grid cells """
    # Buffer by 'b' grid cells to ensure enough room for interpolation
    # But ensure data boundaries are not exceeded
    x_start = max(0, lon_indices[0] - b)
    x_end = min(len(grid.coords[0][0]), lon_indices[1] + b + 1)
    y_start = max(0, lat_indices[0] - b)
    y_end = min(len(grid.coords[0][1]), lat_indices[1] + b + 1)

    lon_slice = slice(x_start, x_end)
    lat_slice = slice(y_start, y_end)

    return lon_slice, lat_slice


def clip_grid_to_indices(grid: "ESMF.Grid",
                         lon_indices:'np.ndarray',
                         lat_indices:'np.ndarray') -> "ESMF.Grid":

    lon_slice, lat_slice = get_buffered_slices(grid, lon_indices, lat_indices)

    sliced_lon = grid.coords[0][0][lon_slice, lat_slice]
    sliced_lat = grid.coords[0][1][lon_slice, lat_slice]

    new_grid = grid_create_from_coordinates_periodic(sliced_lon[:,0],sliced_lat[0,:])

    return new_grid


def clipped_grid_indices(grid: "ESMF.Grid", bbox: "BoundingBox") -> tuple:
    """ Get indices of grid that are needed to cover bbox"""

    latitudes = grid.coords[0][1]
    longitudes = grid.coords[0][0]
    
    # Consider cases where all points within a single grid cell
    new_bbox = BoundingBox(max(longitudes[longitudes <= bbox.xmin]),
                           min(longitudes[longitudes >= bbox.xmax]),
                           max(latitudes[latitudes <= bbox.ymin]),
                           min(latitudes[latitudes >= bbox.ymax]))

    valid_lat = np.where((latitudes[0,:] >= new_bbox.ymin) & (latitudes[0,:] <= new_bbox.ymax))[0]
    valid_lon = np.where((longitudes[:,0] >= new_bbox.xmin) & (longitudes[:,0] <= new_bbox.xmax))[0]

    return valid_lon, valid_lat
