# %%
"""
# Imports and set up
"""

import numpy as np
import random
from copy import deepcopy

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
    def __init__(self, x, B):
        """
        Create a point p on the lattice p = Bx

        Parameters:
        x (np.array): The vector x s.t. Bx = p
        B (np.array): The basis vectors
        """
        self.x = x
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

# B = read_basis('task/latticeBasis.txt')

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

# write_output(B, np.array([14,5,5,23,24,423,4,23,5,]), 100, np.array([1,3,4,56,4]), 'test.txt')

# %%
"""
# Sieving functions
"""

def generate_basis_vectors(B, n):
    """
    Function to generate points which are the basis vectors from B
    
    Parameters:
    B (np.array): Basis for lattice
    n (int): number of points to generate
 
    Returns:
    set: set of distinct points
    SVPoint: short vector on the lattice
    """
    points = set()

    for i in range(B.shape[0]):
        tmp = np.zeros_like(B[1])
        tmp[i] = 1
        new_p = SVPoint(tmp, B)
        points.add(new_p)

    SV = min(points, key=lambda x: x.norm)
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

    while len(points) < n:
        x = np.random.randint(l, h, size=(len(B[0]),))
        new_p = SVPoint(x, B)
        assert (new_p.p == np.dot(B, new_p.x)).all()
        if new_p.norm != 0:
            points.add(new_p)

    SV = min(points, key=lambda x: x.norm)

    return points, SV


def find_difference(p, q, B):
    """
    Function to find the point which is the difference of p and another q
    
    Parameters:
    p (SVPoint): Point p
    q (SVPoint): Point q
    B (np.array): Basis for lattice
    
    Returns:
    SVPoint: average of self and other on B
    """
    diff_x = p.x - q.x 
    return SVPoint(diff_x, B)


def find_average(p, q, B):
    """
    Function to find the point which is the average of p and q
    
    Parameters:
    p (SVPoint): Point p
    q (SVPoint): Point q
    B (np.array): Basis for lattice
    
    Returns:
    SVPoint: average of p and q on B
    None: if the average isn't a point, it returns None
    """
    for i, j in zip(p.x, q.x):
        if (i - j) % 2 != 0:
            return None    
    avg_x = np.array((p.x + q.x)/2)
    return SVPoint(avg_x, B)


def find_modified_average(p, q, B):
    """
    Function to find the point which is the average of self and another point
    
    Parameters:
    p (SVPoint): Point p
    q (SVPoint): Point q
    B (np.array): Basis for lattice
    
    Returns:
    SVPoint: average of p and q
    """
    for i, (ui, vi) in enumerate(zip(p.x, q.x)):
        if (ui - vi) % 2 != 0:
            tmp = q.x
            tmp[i] += random.choice([-1,1])
            q.x = tmp
    
    avg_x = np.array((p.x + q.x)/2)
    return SVPoint(avg_x, B)
    
import time

def augment(points, SV, B, n, p, f=find_modified_average):
    """
    Function to find n more points
    
    Parameters:
    points (set): Points to consider
    SV (SVPoint): current shortest vector on the lattice
    B (np.array): Basis for lattice
    n (int): number of points to find
    p (float): proportion of the old points to keep for the new set
    f (function): function to combine two vectors

    Returns:
    set: more points
    SVPoint: current shortest vector on the lattice
    """
    points_list = list(points)
    points_list = sorted(points_list, key=lambda x: x.norm)

    # Generate a new set. Options: empty set, some from the old set, include shortest vector
    # new_points = set(random.sample(points, int(p*len(points))))
    # new_points = set()
    # new_points.add(new_SV)

    # Keep the first p shortest vectors
    points_list = points_list[:-int((1-p) * len(points_list))]

    # Fill up the new set with new vectord
    s = time.time()
    print('starting')
    
    while len(new_points) < n:
        print('loop')
        p = random.choice(points_list)
        q = random.choice(points_list)

        if time.time() - s > 20:
            print(p.p)
            print(q.p)

        if p != q:
            new_p = f(p, q, B)
            if new_p.norm < p.norm and new_p.norm < q.norm:
                new_points.add(new_p)
    print('end')
    points_list = list(points)
    points_list = sorted(points_list, key=lambda x: x.norm)[:-1]
    points = set(points_list)
    
    SV = min(points, key=lambda x: x.norm)

    return new_points, SV

