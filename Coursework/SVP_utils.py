# %%
"""
# Imports and set up
"""

import numpy as np
import random
from copy import deepcopy

B = np.array([
 [37, 20, 96, 20, 34, 64, 82, 56, 47, 21, 50, 49],
 [39, 24, 19, 49, 82, 97, 88, 84, 41, 51, 36, 74],
 [19, 56, 37, 73,  4, 12, 72, 18, 46,  8, 54, 94],
 [13, 46, 26,  8, 83, 71, 45, 84, 21, 32, 53, 80],
 [65, 39, 25, 56, 52, 44, 84, 30, 69, 33, 13,  5],
 [59, 56, 90,  1, 42, 58, 90, 92,  2,  6,  7, 80],
 [18, 14, 26, 31, 91, 93, 77, 64, 95, 36, 23,  5],
 [11, 58, 22, 51, 90, 13, 93, 43, 21, 81, 12, 77],
 [42, 65, 99,  6, 23, 43, 94, 30, 37, 66, 34, 66],
 [99, 31, 24, 44, 18, 58, 17, 27, 70, 88, 59, 11],
 [30, 43, 21, 70, 48, 47, 13, 93, 94, 48, 69, 58],
 [ 7, 12, 94, 88, 59, 95, 43, 62, 71, 36, 91, 70]
])

# %%
"""
# Gram Schmidt function for Basis Orthogonalisation
"""

def gram_schmidt(B):
    """
    Function to orthogonalise a basis B
    
    Parameters:
    B (np.array): Input lattice basis vectors
    
    Returns:
    np.array: Orthogonalised lattice basis vectors
    """

    def proj(u, v): 
        return u * np.dot(u, v)
        
    V = B.T
    us = []
    for i in range(1,len(B[0])+1):
        us.append( V[i-1] - sum([proj(us[j], V[i-1]) for j in range(i-1)]) )
    
    B_ = np.array(us)
    return B_.T

# B = gram_schmidt(B)
# %%
"""
# SVPoint class
"""

class SVPoint():
    def __init__(self, x):
        """
        Create a point p on the lattice p = Bx

        Parameters:
        x (np.array): The vector x s.t. Bx = p
        """
        self.x = np.array(x)
        self.p = np.dot(B, x)
        self.norm = self.find_norm()

    def __eq__(self, other):
        c = self.p == other.p
        return c.all()
        
    def __neq__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return f"p={self.p},\nx={self.x},\nnorm={self.norm}\n"

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(tuple(self.x))

    def find_norm(self):
        """
        Function to find the length of a vector/point
        
        Parameters:
        p (np.array): vector
        
        Returns:
        int: length of the vector
        """
        return np.linalg.norm(self.p)

# %%
"""
# File input output functions
"""

import numpy as np

def read_basis(filename):
    with open(filename) as f:
        lines = [l[2:-3] for l in f.readlines()[5:] if l[0] not in '#\n']

    B = [[int(ele) for ele in l.split(', ')] for l in lines]
    B = np.array(B)
    return B

B = read_basis('task/latticeBasis.txt')

def write_output(B, u, norm, x, filename):
    out_str = f"""
# DISCLAIMER
# This is an example of short vector output
# PLease include:
# -- the basis B
# -- the point on the lattice u
# -- the norm of u
# -- the vector x such that u = Bx
# Again, you needn't write the disclaimer

B =
{B}

u =
{u}

norm =
{norm}

x =
{x}
    """
    with open(filename,'w') as f:
        f.write(out_str)

write_output(B, np.array([14,5,5,23,24,423,4,23,5,]), 100, np.array([1,3,4,56,4]), 'test.txt')

# %%
"""
# Sieving functions
"""

def generate_basis_vectors(B, n, SV):
    """
    Function to generate points which are the basis vectors from B
    
    Parameters:
    B (np.array): Basis for lattice
    n (int): number of points to generate
    SV (SVPoint): short vector on the lattice
 
    Returns:
    set: set of distinct points
    SVPoint: short vector on the lattice
    """
    points = set()

    for i in range(B.shape[0]):
        tmp = np.zeros_like(B[1])
        tmp[i] = 1
        new_p = SVPoint(tmp)
        points.add(new_p)
        if 0 < new_p.norm < SV.norm:
            SV = new_p

    return points, SV


def generate_random_points(B, n, l=-2, h=3):
    """
    Function to generate random points from the lattice defined by B
    
    Parameters:
    B (np.array): Basis for lattice
    n (int): number of points to generate
    l (int): low number for range of basis vector multiples
    h (int): high number for range of basis vector multiples
    
    Returns:
    set: points
    SVPoint: short vector on the lattice
    """
    points = set()
    SV = SVPoint(np.random.randint(l, h, size=(len(B[0]),)))

    while len(points) < n:
        new_p = SVPoint(np.random.randint(l, h, size=(len(B[0]),)))
        points.add(deepcopy(new_p))
        if 0 < new_p.norm < SV.norm:
            # SV = SVPoint(new_p.x)
            SV = deepcopy(new_p)

    return points, SV


def find_difference(p, q):
    """
    Function to find the point which is the difference of p and another q
    
    Parameters:
    p (SVPoint): Point p
    q (SVPoint): Point q
    
    Returns:
    SVPoint: average of self and other on B
    """
    diff_x = p.x - q.x 
    return SVPoint(diff_x)


def find_average(p, q):
    """
    Function to find the point which is the average of p and q
    
    Parameters:
    p (SVPoint): Point p
    q (SVPoint): Point q
    
    Returns:
    SVPoint: average of p and q on B
    None: if the average isn't a point, it returns None
    """
    for i, j in zip(p.x, q.x):
        if (i - j) % 2 != 0:
            return None    
    avg_x = np.array((p.x + q.x)/2)
    return SVPoint(avg_x)


def find_modified_average(p, q):
    """
    Function to find the point which is the average of self and another point
    
    Parameters:
    p (SVPoint): Point p
    q (SVPoint): Point q
    
    Returns:
    SVPoint: average of p and q
    """
    for i, (ui, vi) in enumerate(zip(p.x, q.x)):
        if (ui - vi) % 2 != 0:
            q.x[i] += 1
    
    avg_x = np.array((p.x + q.x)/2)
    return SVPoint(avg_x)
    
def augment(points, n, SV, f=find_modified_average):
    """
    Function to find n more points
    
    Parameters:
    points (set): Points to consider
    n (int): number of points to find
    SV (SVPoint): current shortest vector on the lattice
    f (function): function to combine two vectors

    Returns:
    set: more points
    SVPoint: current shortest vector on the lattice
    """
    new_SV = deepcopy(SV)
    points_list = list(points)

    # Generate a new set. Options: empty set, some from the old set, include shortest vector
    new_points = set(random.sample(points, int(0.5*n)))
    # new_points = set()
    new_points.add(new_SV)

    # Fill up the new set with averages
    while len(new_points) < n:
        p = random.choice(points_list)
        q = random.choice(points_list)

        if p != q:
            new_p = f(p, q)
            new_points.add(new_p)
            if 0 < new_p.norm < new_SV.norm:
                new_SV = SVPoint(new_p.x)

    return new_points, new_SV

