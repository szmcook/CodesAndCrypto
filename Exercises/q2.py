from copy import copy

def plus(a, b):
    return (a+b)%2

def mult(a, b):
    return (a*b)%2

def concat(a,b,c,d):
    return str(a)+str(b)+str(c)+str(d)

def split(a):
    '''takes an int or a string and returns 4 ints'''
    s = str(a)
    a,b,c,d = int(s[0]), int(s[1]), int(s[2]), int(s[3])
    return a,b,c,d

def halve(a):
    '''takes a string and returns L, R'''
    L = a[:4]
    R = a[4:]
    return L, R

def f(R, K):
    '''each input is a 0 or a 1'''
    r1,r2,r3,r4 = split(R)
    k1,k2,k3,k4 = split(K)
    a = mult(plus(r1, k1), k2)
    b = mult(plus(r2, k2), k3)
    c = mult(plus(r3, k3), k4)
    d = mult(plus(r4, k4), k1)
    res = concat(a,b,c,d)
    return res

def xor(a, b):
    '''a, b are 4 bit binary numbers represented as strings'''
    a1,b1,c1,d1 = split(a)
    a2,b2,c2,d2 = split(b)

    a = plus(a1,a2)
    b = plus(b1,b2)
    c = plus(c1,c2)
    d = plus(d1,d2)

    res = concat(a,b,c,d)
    return res


def encrypt(M, K):
    L0, R0 = halve(M)
    K1, K2 = halve(K)

    L1 = copy(R0)
    R1 = xor(L0, f(R0, K1))

    # print(L1, R1)

    L2 = copy(R1)
    R2 = xor(L1, f(R1, K2))

    print(L2, R2)
    return L2 + R2


def decrypt(C, K):
    L2, R2 = halve(C)
    K1, K2 = halve(K)

    L1, R1 = xor(R2, f(L2, K2)), L2
    L0, R0 = xor(R1, f(L1, K1)), L1

    print(L0, R0)
    return L0 + R0


# M = "01110001"
# K = "01011101"
# C = encrypt(M, K)
# M = decrypt(C, K)

