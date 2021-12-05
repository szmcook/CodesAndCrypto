# %%
"""
# Lattice based cryptography
"""

# %%
import numpy as np
from SVP_utils import generate_random_points, augment_mod_avg, SVPoint

# %%
B = [
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
]

B = np.array(B)

# %%
# Upper bound on the shortest vector by minkowski's theorem
np.sqrt([12])[0] * (np.linalg.det(B) ** (1/12))

# %%
def sieve(B):
    """
    Function to sieve for short vectors on a lattice with basis B
    
    Parameters:
    B (np.array): Input lattice basis vectors
    SV (np.array): short vector on the lattice
    L (int): length of the short vector
    
    Returns:
    np.array: short vector on the lattice
    int: length of the short vector
    """
    # 0 Orthogonalise the basis with gram-schmidt

    # 1 Generate n points on the lattice
    n = 10000
    points, SV = generate_random_points(B, n)

    # 2 Replace the points with new ones until done TODO when is this done?
    for _ in range(1):
        
        if _ % 300 == 0:
            print(_, SV)
            assert (SV.p == np.dot(B, SV.x)).all()

        points, SV = augment_mod_avg(points, n, SV)
    
    return SV

SV = sieve(B)
assert (SV.p == np.dot(B, SV.x)).all()

with open("SVP_out.txt", "w") as f:
    f.write(f"{SV}")