# Python 3 program to find a prime factor of composite using
# Pollard's Rho algorithm
import math

# Function to calculate (base^exponent)%modulus
def f2(x):
    # initialize result
    exponent = 2
    result = 1
    while (exponent > 0):
        # if y is odd, multiply base with result
        if (exponent & 1):
            result = (result * x) % n
        # exponent = exponent/2
        exponent = exponent >> 1
        # base = base * base
        x = (x * x) % n
    return result + 1

# method to return prime divisor for n
def pollard_rho(n):
	x = 1
	y = x
	p = 1

	while p == 1:
		x = f2(x) % n
		y = f2( f2(y) % n ) % n
		p = math.gcd(abs(x - y), n)

	return p

n = 7171
print(f"{n} -> {pollard_rho(n)}")