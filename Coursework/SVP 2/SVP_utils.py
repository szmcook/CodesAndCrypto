# %%
"""
# Imports and set up
"""

import numpy as np
import random

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

def find_norm(x, B):
    """
    Function to find the length of a vector/point
    
    Parameters:
    x (np.array): vector
    B (np.array): Basis for lattice

    Returns:
    int: length of the vector Bx
    """
    p = np.dot(B, x)
    return np.linalg.norm(p)


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
    tuple: short vector on the lattice
    """
    points = set()

    while len(points) < n:
        x = np.random.randint(l, h, size=(len(B[0]),))
        norm = find_norm(x, B)

        if norm != 0:
            points.add( (x, norm) )

    SV = min(points, key=lambda x: x[1])

    return points, SV


def find_difference(p, q, B):
    """
    Function to find the point which is the difference of p and another q
    
    Parameters:
    p (tuple): Point p
    q (tuple): Point q
    B (np.array): Basis for lattice

    Returns:
    tuple: difference of p and q on B
    """
    diff_x = p[1] - q[1]
    norm = find_norm(diff_x, B)
    return (diff_x, norm)


def find_average(p, q, B):
    """
    Function to find the point which is the average of p and q
    
    Parameters:
    p (tuple): Point p
    q (tuple): Point q
    B (np.array): Basis for lattice
    
    Returns:
    tuple: average of p and q on B
    None: if the average isn't a point, it returns None
    """
    for i, j in zip(p[0], q[0]):
        if (i - j) % 2 != 0:
            return None    
    avg_x = np.array((p[0] + q[0])/2)
    norm = find_norm(avg_x, B)
    return (avg_x, norm)


def find_modified_average(p, q):
    """
    Function to find the point which is the average of self and another point
    
    Parameters:
    p (tuple): Point p
    q (tuple): Point q
    
    Returns:
    tuple: average of p and q
    """
    for i, (ui, vi) in enumerate(zip(p[0], q[0])):
        if (ui - vi) % 2 != 0:
            q[0][i] += 1
    
    avg_x = np.array((q[0] + q[0])/2)
    return SVPoint(avg_x, B)
    
def augment(points, n, SV, f=find_modified_average):
    """
    Function to find n more points
    
    Parameters:
    points (set): Points to consider
    n (int): number of points to find
    SV (tuple): current shortest vector on the lattice
    f (function): function to combine two vectors

    Returns:
    set: more points
    tuple: current shortest vector on the lattice
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
                new_SV = SVPoint(new_p.x, B)

    return new_points, new_SV

