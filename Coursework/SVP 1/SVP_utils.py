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
    list: points ordered by norm
    """
    points = []

    while len(points) < n:
        x = np.random.randint(l, h, size=(len(B[0]),))
        new_p = SVPoint(x, B)

        if new_p.norm != 0 and new_p not in points:
            points.append(new_p)

    points.sort(key=lambda x: x.norm)
    return points


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
    tmp = deepcopy(q.x)
    for i, (ui, vi) in enumerate(zip(p.x, q.x)):    
        if (ui - vi) % 2 != 0:
            tmp[i] += random.choice([-1,1])
    
    avg_x = np.array((p.x + tmp)/2)
    new_p = SVPoint(avg_x, B)
    return new_p


def find_modified_average_random(p, q, B):
    """
    Function to find a random point on the lattice from the modified average

    Parameters:
    p (SVPoint): Point p
    q (SVPoint): Point q
    B (np.array): Basis for lattice
    
    Returns:
    SVPoint: random point with x values between -3 and 3
    """
    tmp = deepcopy(q.x)
    for i, (ui, vi) in enumerate(zip(p.x, q.x)):    
        if (ui - vi) % 2 != 0:
            tmp[i] += random.choice([-1,1])

    # Add some random noise
    for _ in range(3):
        tmp[random.randint(0,11)] += random.choice([-2, 2])

    avg_x = np.array((p.x + tmp)/2)

    new_p = SVPoint(avg_x, B)
    return new_p
    
# import time

def augment(points, B, n, p, generate_point=find_modified_average):
    """
    Function to find n more points
    
    Parameters:
    points (list): Ordered list of points to consider
    B (np.array): Basis for lattice
    n (int): number of points to put into the list to return
    p (float): proportion of the old points to keep for the new list
    generate_point (function): function to combine two vectors and generate a new one

    Returns:
    list: more points
    """

    # Generate a new set. Options: empty set, some from the old set, include shortest vector
    # new_points = set(random.sample(points, int(p*len(points))))
    # new_points = set()
    # new_points.add(new_SV)

    # Keep the first p shortest vectors
    new_points = points[:int(p * len(points))]

    attempts = 0

    # Fill up the new set with new vectors
    # while len(new_points) < n and attempts < 10000:
    while len(new_points) < n:
        attempts += 1

        p1 = random.choice(points[:int(p * len(points))])
        p2 = random.choice(points[:int(p * len(points))])

        # We require that: new_p is not the 0 vector, is not in the list, is shorter than p1 and p2
        new_p = generate_point(p1, p2, B)
        if (new_p.norm == 0):
            # print("0 norm")
            continue
        if (p1 == p2):
            # print('parents are the same')
            continue
        if (new_p.norm > p1.norm and new_p.norm > p2.norm): # TODO experiment with this being an and or an or MAKES NO DIFFERENCE?
            # print('new_p longer than parents')
            continue
        if (new_p in new_points):
            # print('new_p in new_points')
            continue
        new_points.append(new_p)

    new_points.sort(key=lambda x: x.norm)
    return new_points