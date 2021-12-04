# meet in the middle attack

# Method:
# 1. Compute all possible encryptions of the sample M (with all possible keys)
# 2. Compute all possible decryptions of the sample C
# 3. Look for matches

from q2p1 import encrypt, decrypt

def generate_keyspace(n):
    '''generates strings representing the keys with n bits length'''

    def binary(i):
        b = bin(i)[2:]
        s = "0"*(n-len(b)) + b
        return s

    keys = []

    for i in range(2**n):
        keys.append(binary(i))

    return keys

def meet_in_the_middle(M, C, keyspace):
    '''takes a sample of plaintext, sampletext and returns a key from the keyspace'''

    encryptions = {}
    decryptions = {}
    
    for key1 in keyspace:
        e = encrypt(M, key1)
        encryptions[key1] = e

    for key2 in keyspace:
        d = decrypt(C, key2)
        decryptions[key2] = d

    candidates = []
    for key1 in keyspace:
        if encryptions[key1] in decryptions.values():
            key2 = list(decryptions.keys())[list(decryptions.values()).index(encryptions[key1])]
            candidates.append((key1, key2))

    print(candidates)
    print(len(candidates))

    return candidates[0]

key_size = 8
M = "01110001"
C = "10000111"
(key1, key2) = meet_in_the_middle(M, C, generate_keyspace(key_size))

print(encrypt(M, key1))
print(decrypt(C, key2))

