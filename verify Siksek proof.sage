cremona_labels = ["128A1", "128B1", "128C1", "128D1"] # Cremona labels of the four elliptic curves corresponding to newforms of level 128

def a(l, E): # a_l(E) function for elliptic curve E over finite field of order l
    return l+1-E.count_points()

def refine(l, F):
    """
    Inputs:
        integer l congruent to 1 mod p and +/-1 mod 8.
        string F: one of the four Cremona labels in cremona_labels
    Behavior:
        If l is a prime satisfying conditions 1-4 in Proposition 15.7.1, computes the set R_l(F),
        and replaces current_set with its intersection with R_l(F) and appends l to used_primes
        if it makes current_set smaller.
    """
    global current_set, used_primes
    if not l.is_prime(): return
    alF = a(l, EllipticCurve(GF(l), F))
    if alF%p == (l+1)%p or alF%p == (-l-1)%p: return # condition 3
    n = (l-1)//p # value of n so that l = np+1
    theta = sqrt(mod(2, l))
    if (1+theta)^n == 1: # condition 4
        if (1-theta)^n == 1: return # If theta doesn't work, try -theta. If -theta doesn't work either, return.
        else: theta = -theta
    X1 = [d for d in Integers(l) if (d^2-2)^n == 1] # The set X_l', d is for delta
    R.<x,y> = GF(l)[]
    X = [d for d in X1 if a(l, EllipticCurve(y^2-(x^3+2*d*x^2+2*x)))%p == alF%p] # The set X_l
    g = mod(primitive_root(l), l)
    phi = lambda x: mod(discrete_log(x, g), p) # The map phi: the composition of the discrete log with reduction mod p
    den = phi(1+theta)^(-1) # denominator used in R_l(F)
    R = set(phi(d+theta)*den for d in X) # The set R_l(F)
    if not current_set.issubset(R): # If current_set is not a subset of R_l(F), replace with intersection and append l to used_primes
        used_primes.append(l)
        current_set = current_set.intersection(R)
    
dict_of_primes_used = {} # dictionary mapping tuple (p, F) for p a prime and F in cremona_labels to the list of primes used in the proof
p = 11
next_thousand = 1000 # Every time p exceeds the next multiple of 1000, we write to "output.txt"
while p < 20000:
    # We want l = 1 (mod p), l = 1 or 7 (mod 8). We can find the congruence class mod 8p guaranteed by the Chinese Remainder Theorem
    # If l = np+1 = 1 mod 8 then n = 0 (mod 8), so l = 1 (mod 8p)
    # If l = np+1 = -1 mod 8 then if p^{-1} is the inverse of p mod 8, n = -2p^{-1} (mod 8), then l = -2p^{-1}p + 1 (mod 8p)
    # Fix s1 and s2 to be the first integers > 1 in these classes mod 8p, and then we will repeatedly add 8p to them.
    s1 = Integer(mod(-2*mod(p, 8)^(-1),8*p)*p+1)
    s2 = 1+8*p
    target = {mod(1, p), mod(-1, p)} # target set for intersection of R_l(F)
    for F in cremona_labels:
        current_set = set(Integers(p))
        used_primes = []
        r1 = s1
        r2 = s2
        while not current_set.issubset(target): # loop to repeatedly refine current_set until it is a subset of target
            refine(r1, F)
            if current_set.issubset(target): break
            refine(r2, F)
            r1 += 8*p
            r2 += 8*p
        dict_of_primes_used[(p, F)] = used_primes
    p = p.next_prime()
    if p > next_thousand:
        print(f"Done up to {p}")
        with open("output.txt", 'a') as outfile:
            for key, val in dict_of_primes_used.items():
                outfile.write(f"p={key[0]}, F={key[1]}: {', '.join(map(str, val))}\n")
        dict_of_primes_used = {}
        next_thousand += 1000
