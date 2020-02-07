import numpy as np

from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform



def make_othro_array(nx, ny, dx, dy, sp=[0, 0, 0]):
    """This function makes an nx by ny orthogonal array of points,
        spaced dx and dy, and starting from sp. The z coordinate
        will be assumed to be the same as the starting point sp. 

        Parameters
        nx(int): The number of points in the x direction
        ny(int): The number of points in the y direction
        dx(float): The distance between points in the x direction
        dy(float): The distance between points in the y direction
        sp(list): The x, y, z coordinates of the starting point

        Returns
        points(ndarray): The point array
    """

    points = [[i * dx + sp[0], j * dy + sp[1], sp[2]] for j in range(ny) for i in range(nx)]
    return np.array(points)

def point_array_distance(points):
    """This function returns an m x m array of distances between
    m points. The diagonal will be zero and the off diagonal terms
    will be the distances between each pair of points. 

    Parameters
    points(ndarray): The m x 1 array of points to calculate distances.
    """
    d = pdist(points, 'euclidean')
    return squareform(d)