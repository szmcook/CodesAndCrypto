
"""
# Imports and set up
"""

import numpy as np
import random
from copy import deepcopy
import re


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

"""
# SVPoint class
"""

class SVPoint():
    def __init__(self, x, B):
        """
        Create a vector p on the lattice p = Bx

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
        Function to find the length of a vector
        
        Returns:
        int: length of the vector
        """
        return np.linalg.norm(self.p)


"""
# File input output functions
"""


def read_basis(filename):
    with open(filename) as f:
        lines = [re.findall('[0-9]+', l) for l in f.readlines() if l[0] not in '#B\n']
    return np.array(lines, dtype=np.int64)


def write_output(B, u, norm, x, filename):
    out_str = f"""# DISCLAIMER
# This is an example of short vector output
# PLease include:
# -- the basis B
# -- the vector on the lattice u
# -- the norm of u
# -- the vector x such that u = Bx
# Again, you needn't write the disclaimer

B =
{np.array2string(B, separator=", ")}

u =
{np.array2string(np.atleast_2d(u).T, separator=", ")}

norm =
{norm}

x =
{np.array2string(np.atleast_2d(x).T, separator=", ")}
    """
    with open(filename,'w') as f:
        f.write(out_str)


"""
# Sieving functions
"""

def generate_basis_vectors(B, n):
    """
    Function to generate vectors which are the basis vectors from B
    
    Parameters:
    B (np.array): Basis for lattice
    n (int): number of vectors to generate
    
    Returns:
    list: vectors ordered by norm
    """
    vectors = []

    while len(vectors) < n:
        vs = random.randint(0, 4)
        for _ in range(vs):
            tmp = np.zeros_like(B[1])
            i = random.randint(0, 11)
            tmp[i] = random.choice([-1, 1])
            new_p = SVPoint(tmp, B)
            if new_p.norm != 0 and new_p not in vectors:
                vectors.append(new_p)

    vectors.sort(key=lambda x: x.norm)
    return vectors


def generate_random_vectors(B, n, l=-2, h=3):
    """
    Function to generate random vectors from the lattice defined by B
    
    Parameters:
    B (np.array): Basis for lattice
    n (int): number of vectors to generate
    l (int): low number for range of basis vector multiples
    h (int): high number for range of basis vector multiples
    
    Returns:
    list: vectors ordered by norm
    """
    vectors = []

    while len(vectors) < n:
        x = np.random.randint(l, h, size=(len(B[0]),))
        new_p = SVPoint(x, B)

        if new_p.norm != 0 and new_p not in vectors:
            vectors.append(new_p)

    vectors.sort(key=lambda x: x.norm)
    return vectors


def find_random(p, q, B):
    """
    Function to find a random vector on the lattice
    
    Parameters:
    p (SVPoint): vector p
    q (SVPoint): vector q
    B (np.array): Basis for lattice
    
    Returns:
    SVPoint: random vector on lattice defined by B
    """
    vectors = []    
    while len(vectors) < 1:
        x = np.random.randint(-4, 4, size=(len(B[0]),))
        new_p = SVPoint(x, B)
        if new_p.norm != 0:
            vectors.append(new_p)

    return vectors[0]


def find_difference(p, q, B):
    """
    Function to find the vector which is the difference of p and another q
    
    Parameters:
    p (SVPoint): vector p
    q (SVPoint): vector q
    B (np.array): Basis for lattice
    
    Returns:
    SVPoint: average of self and other on B
    """
    diff_x = p.x - q.x 
    return SVPoint(diff_x, B)


def find_average(p, q, B):
    """
    Function to find the vector which is the average of p and q
    
    Parameters:
    p (SVPoint): vector p
    q (SVPoint): vector q
    B (np.array): Basis for lattice
    
    Returns:
    SVPoint: average of p and q on B
    None: if the average isn't a vector, it returns None
    """
    for i, j in zip(p.x, q.x):
        if (i - j) % 2 != 0:
            return None    
    avg_x = np.array((p.x + q.x)/2)
    return SVPoint(avg_x, B)


def find_modified_average(p, q, B):
    """
    Function to find the vector which is the average of self and another vector
    
    Parameters:
    p (SVPoint): vector p
    q (SVPoint): vector q
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
    Function to find a random vector on the lattice from the modified average

    Parameters:
    p (SVPoint): vector p
    q (SVPoint): vector q
    B (np.array): Basis for lattice
    
    Returns:
    SVPoint: random vector with x values between -3 and 3
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
    

def augment(vectors, B, n, p, generate_vector=find_modified_average, timeout=None):
    """
    Function to find n more vectors
    
    Parameters:
    vectors (list): Ordered list of vectors to consider
    B (np.array): Basis for lattice
    n (int): number of vectors to put into the list to return
    p (float): proportion of the old vectors to keep for the new list
    generate_vector (function): function to combine two vectors and generate a new one
    timeout (int): seconds to timeout after

    Returns:
    list: more vectors
    """

    # Generate a new set. Options: empty set, random sample from the old set, shortest from the old set
    # new_vectors = set()
    # new_vectors = set(random.sample(vectors, int(p*len(vectors))))
    # new_vectors.add(new_SV)

    # Keep the first p shortest vectors
    new_vectors = vectors[:int(p * len(vectors))]

    # Fill up the new set with new vectors
    while len(new_vectors) < n:
        p1 = random.choice(vectors[:int(p * len(vectors))])
        p2 = random.choice(vectors[:int(p * len(vectors))])

        # We require that: new_p is not the 0 vector, is not in the list, is shorter than p1 and p2
        new_p = generate_vector(p1, p2, B)
        if generate_vector==find_random:
            new_vectors.append(new_p)
        elif (new_p.norm == 0):
            # print("0 norm")
            continue
        elif (p1 == p2):
            # print('parents are the same')
            continue
        elif (new_p.norm > p1.norm) and (new_p.norm > p2.norm): # Using or in place of and makes little difference here.
            # print('new_p longer than parents')
            continue
        elif (new_p in new_vectors):
            # print('new_p in new_vectors')
            continue
        else:
            new_vectors.append(new_p)

    new_vectors.sort(key=lambda x: x.norm)
    return new_vectors