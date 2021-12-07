# %%
"""
# ECC Functions and Classes
"""
# %%
def EEA(a,b):
    rs = [a, b]
    qs = [0, 0]
    ss = [1, 0]
    ts = [0, 1]

    i = 1 

    while rs[i] > 0:
        i = i + 1
        qs.append( rs[i-2] // rs[i-1] )
        rs.append( rs[i-2] - qs[i] * rs[i-1] )
        ss.append( ss[i-2] - qs[i] * ss[i-1] )
        ts.append( ts[i-2] - qs[i] * ts[i-1] )

    #      gcd      x       y
    # s.t. ax + by = gcd
    return ts[i-1]


def inverseModp(x, p):
    y = x % p
    t = EEA(p, y)
    return t % p
# %% 
"""
# Point and Elliptic Curve Classes
"""
# %%
class Point():
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        if self.x == other.x and self.y == other.y:
            return True
        else:
            return False
        
    def __neq__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return f"({self.x}, {self.y})"

    def __repr__(self):
        return str(self)


class ECC():
    def __init__(self, p, a, b):
        """
        Create an Elliptic Curve with y**3 = x**3 + ax + b

        Parameters:
        p (int): The modulo prime
        a (int): a
        b (int): b

        Returns:
        int: The value to return
        """
        self.p = p
        self.a = a
        self.b = b
        assert 4*a**3 + 27*b**2 != 0

    def ECPointAddition(self, P, Q):
        """
        Add P to Q

        Parameters:
        P (Point): Point from the curve
        Q (Point): Point from the curve
        
        Returns:
        Point: The value of P + Q
        """
        # If P is infinity, return Q
        if (P.x, P.y) == (self.p, self.p):
            return Q
        # If Q is infinity, return P
        if (Q.x, Q.y) == (self.p, self.p):
            return P
        # If P = -Q, return infinity
        if P.x == Q.x and P.y == (-Q.y) % self.p:
            return Point(self.p, self.p)
        # If P = -Q, return infinity
        if (P.x, P.y) == (Q.x, Q.y) and P.y == 0:
            return Point(self.p, self.p)
        # 2P, for P != -P # POINT DOUBLING
        if (P.x, P.y) == (Q.x, Q.y) and P.y != 0:
            z = (3 * P.x**2 + self.a) * inverseModp(2*P.y, self.p)
            x3 = (z**2 - 2*P.x) % self.p
            y3 = (z*(P.x - x3) - P.y) % self.p
            return Point(x3, y3)
        else:
            z = ( (Q.y - P.y) * inverseModp( Q.x-P.x, self.p ) ) % self.p
            x3 = ( z**2 - P.x - Q.x ) % self.p
            y3 = ( z * (P.x - x3) - P.y ) % self.p
            return Point(x3, y3)

    def ECPointDoubling(self, P):
        """
        Double P

        Parameters:
        P (Point): Point from the curve
        
        Returns:
        Point: The value of 2P
        """
        if P.y == -P.y:
            return Point(self.p, self.p)
        else:
            z = ( (3*(P.x**2)%self.p + self.a) * inverseModp( 2*P.y, self.p ) ) % self.p
            x3 = ( (z**2) - 2*P.x ) % self.p
            y3 = ( z * (P.x - x3) - P.y ) % self.p
            return Point(x3, y3)


    def ECPointMult(self, k, P):
        """
        Multiply P by k (quickly)

        Parameters:
        P (Point): Point from the curve
        k (int): Scalar to multiply p by
        
        Returns:
        Point: The value of kP
        """
        sum = Point(self.p, self.p)
        X = P
        # Quick and dirty way to get the binary expansion of x
        f = format(k, 'b')[::-1]
        w = len(f)

        for i in range(w):
            i = int( f[i] )
            if i == 1:
                sum = self.ECPointAddition(sum, X)
            X = self.ECPointDoubling(X)
            
        return sum
# %% 
"""
Binary digit manipulation functions
"""
# %% 
def bits_to_bytes(bits):
    """
    Takes a string of 8 bits and returns a single byte
    """
    x = int(bits, 2)
    b = x.to_bytes(8, byteorder='big')
    print(b)
    return b
# %% 
"""
# File input and output functions
"""
# %% 
def read_ECC_instance(filename):
    with open(filename) as f:
        lines = f.readlines()
    lines = [l for l in lines if l[0]!="#"]
    for l in lines:
        if " = " in l:
            name, val = l.split(" = ")
            if name == "p":
                p = int(val)
            elif name == "a":
                a = int(val)
            if name == "b":
                b = int(val)
            elif name == "P":
                x, y = val.replace("(","").replace(")","").split(",")
                P = Point(int(x), int(y))
            if name == "n":
                n = int(val)
            elif name == "Q":
                x, y = val.replace("(","").replace(")","").split(",")
                Q = Point(int(x), int(y))

    return p, a, b, P, n, Q


def write_output(p, a, b, P, n, Q, c, d, c_, d_, filename):
    out_str = f"""# DISCLAIMER
# Example of output for Basic Pollard Rho
# You needn't write the disclaimer but please write the inputs so it's easier for me to check your results
# IMPORTANT: 
# 1. You may not necessarily get the same collision than I did for this example. 
# 2. A "collision" where c=c' and d=d' is not valid.

Input:
p = {p}
a = {a}
b = {b}
P = ({P.x}, {P.y})
n = {n}
Q = ({Q.x}, {Q.y})

Collision:
c = {c}
d = {d} 
c'= {c_}
d'= {d_}"""
    with open(filename, 'w') as f:
        f.write(out_str)


def write_output_full(p, a, b, P, n, Q, c, d, c_, d_, l, filename):
    out_str = f"""# DISCLAIMER
# Example of output for Basic Pollard Rho
# You needn't write the disclaimer but please write the inputs so it's easier for me to check your results
# IMPORTANT: 
# 1. You may not necessarily get the same collision than I did for this example. 
# 2. A "collision" where c=c' and d=d' is not valid.

Input:
p = {p}
a = {a}
b = {b}
P = ({P.x}, {P.y})
n = {n}
Q = ({Q.x}, {Q.y})

Collision:
c = {c}
d = {d} 
c'= {c_}
d'= {d_}   

Discrete logarithm:
l = {l}"""
    with open(filename, 'w') as f:
        f.write(out_str)