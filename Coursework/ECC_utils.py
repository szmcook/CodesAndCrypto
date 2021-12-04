# %%
"""
# ECC Functions and Classes
"""
# %%
def EEA(a,b):
    r = [a, b]
    q = [0, 0]
    s = [1, 0]
    t = [0, 1]

    i = 1 

    while r[i] > 0:
        i = i + 1
        q.append( r[i-2] // r[i-1] )
        r.append( r[i-2] - q[i] * r[i-1] )
        s.append( s[i-2] - q[i] * s[i-1] )
        t.append( t[i-2] - q[i] * t[i-1] )

    #      gcd      x       y
    # s.t. ax + by = gcd
    return(r[i-1], s[i-1], t[i-1])


def inverseModp(x, p):
    y = x % p
    (r, s, t) = EEA(p, y)
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
        # self.Points = self.points_from_P(Point(0, 7))
        # self.size = len(self.Points)

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
        # If P = -Q
        if (P.x, P.y) == (Q.x, Q.y) and P.y == 0:
            return Point(self.p, self.p)
        # 2P, for P != P
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
        if P.y == -P.y:
            return Point(self.p, self.p)
        else:
            z = ( (3*(P.x**2)%self.p + self.a) * inverseModp( 2*P.y, self.p ) ) % self.p
            x3 = ( (z**2) - 2*P.x ) % self.p
            y3 = ( z * (P.x - x3) - P.y ) % self.p
            return Point(x3, y3)

    def points_from_P(self, P):
        Points = [ (0, Point(self.p, self.p)), (1, P) ]

        k = 1
        Q = Point(P.x, P.y)
        while Q.x != self.p and Q.y != self.p: # and k <= 4000:
            k = k + 1
            Q = self.ECPointAddition(P, Q)
            Points.append( (k, Q) )

        Points.pop()
        return Points

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

    def scalarMultP_by_addition(self, k, P):
        """
        Multiply P by k

        Parameters:
        P (Point): Point from the curve
        k (int): Scalar to multiply p by
        
        Returns:
        Point: The value of kP
        """
        if k == 0:
            return Point(self.p, self.p)
        if k == 1:
            return P
        else:
            i = 1
            Q = P
            while i < k:
                i += 1
                Q = self.ECPointAddition(Q, P)
            return Q
            
# %% 
"""
# File input and output functions
"""

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

    print(p, a, b, n)
    print(P, Q)
    print(type(p), type(a), type(b), type(n))
    print(type(P), type(Q))

    return p, a, b, P, n, Q

# p, a, b, P, n, Q = read_ECC_instance('task/exampleInputRho.txt')


def write_output(p, a, b, P, Q, c, d, c_, d_, filename):
    out_str = f"""
# DISCLAIMER
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
Q = ({Q.x}, {Q.x})

Collision:
c = {c}
d = {d} 
c'= {c_}
d'= {d_}   
    """
    with open(filename, 'w') as f:
        f.write(out_str)

def write_output_full(p, a, b, P, Q, c, d, c_, d_, l, filename):
    out_str = f"""
# DISCLAIMER
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
Q = ({Q.x}, {Q.x})

Collision:
c = {c}
d = {d} 
c'= {c_}
d'= {d_}   

Discrete logarithm:
l = {l}
    """
    with open(filename, 'w') as f:
        f.write(out_str)

write_output_full(16001, 10, 1, Point(1654, 7208), Point(5000, 1283), 159, 161, 271, 349, 1024, 'test.txt')