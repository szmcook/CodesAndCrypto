# The unmasked repeated squaring algorithm
def repeatedSquaring(y,x,n):
    s = 1
    r = 0

    # Quick and dirty way to get the binary expansion of x
    f = format(x, 'b')
    w = len(f)
    
    for k in range(w):
        xk = int( f[k] )
        if xk == 1:
            r = (s * y) % n
        else:
            r = s
        s = (r ** 2) % n

    return r

# The EEA, useful to find inverses mod n
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


def modular_inverse(b, m):
    return EEA(b, m)[1] % m


def mod_exponent_neg(b, e, m):
    '''
    returns b**e % m for a negative e
    
    according to wikipedia
    c = b**e % m = d**(-e) % m
    when e < 0 and b * d = 1

    we know that the modular inverse of 2 % 101 is 51

    '''
    assert e < 0
    d = modular_inverse(b, m)
    assert b*d%m == 1

    c = d**(-e) % m
    return c

    
import random

# I presume that the input is correct, i.e. ex = 1 mod phi(n)
def maskedExponentiation(y,x,n,e,vi = 0,vf = 0):

    # if vi or vf is not given, we generate the pair vf, vi at random
    if vi == 0 or vf == 0:
        r = 0

        while r != 1:
            # vf is random coprime to n
            vf = random.randrange(n)
            (r,s,t) = EEA(n, vf)
            
        vi = repeatedSquaring(t, e, n)
    # if vi and vf are given, we update them
    else:
        vi = (vi * vi) % n
        vf = (vf * vf) % n

    # we now use maked exponentiation
    z = (y * vi) % n
    r = repeatedSquaring(z,x,n)
    r = (vf * r) % n

    return r
