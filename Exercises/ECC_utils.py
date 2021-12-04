# %%
"""
# ECC I
"""

# %%
"""
## 9.1

K is a finite field made froma generator p

It's the set (Zp, +%p, x%p, 0, 1) for p > 3

E defined over K=GF(p) is the collection of points (x,y) in K2 such that

```y^2 = x^3 + ax + b```

```4a^3 + 27b^2 != 0``` for an EC
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
        self.p = p
        self.a = a
        self.b = b
        assert 4*a**3 + 27*b**2 != 0
        # self.Points = self.points_from_P(Point(0, 7))
        # self.size = len(self.Points)

    def ECPointAddition(self, P, Q):
        # If P is infinity, return Q
        if (P.x, P.y) == (self.p, self.p):
            return Q
        # If Q is infinity, return infinity
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

    def scalarMultP(self, k, P):
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

    def curve(self):
        discriminant = ( 4 * (self.a ** 3) + 27 * (self.b ** 2) ) % self.p
        if discriminant == 0:
            raise Exception("discriminant = 0, invalid curve")
        
        curve = [ Point(self.p, self.p) ]
        
        xcurve = []
        for x in range(self.p):
            xonthecurve = (x**3 + self.a*x + self.b) % self.p
            xcurve.append( (x, xonthecurve) )
        
        ycurve = []
        for y in range(self.p):
            yonthecurve = ( y**2 ) % self.p
            ycurve.append( (y, yonthecurve) )

        xcurve.sort(key=lambda x: x[1])

        # TODO FINISH THIS METHOD