def verify_y_greater_than(p, y_bound): # verifies that y can be taken to be greater than y_bound for the given value of p
    theta = sqrt(2) * (1 - 2 / ((1 + sqrt(2))**(-2 / p) + 1))
    cf = continued_fraction(theta)
    lb = p * 2**((3 * p - 7) / 2) - 2 # lower bound on q_{k+1} from (6.4)
    Cp = 2**((3 - p) / 2) / p
    k = 1
    while True:
        if cf.quotient(k + 1) > lb:
            print(f"Verification for {p} failed")
            return False
        k += 1
        L = cf.denominator(k)
        if 1.9 * L**2 >= y_bound:
            print(f"Verified {p}")
            return True
        
verify_y_greater_than(1847, 10^400)

verify_y_greater_than(1861, 10^200)

verify_y_greater_than(1867, 10^200)

for p in [1871, 1873, 1877, 1879]:
    verify_y_greater_than(p, 10^150)

for p in range(1889, 2300):
    if Integer(p).is_prime():
        verify_y_greater_than(p, 10^100)

for p in range(2300, 3559):
    if Integer(p).is_prime():
        verify_y_greater_than(p, 10^50)

for p in range(17, 1832):
    if Integer(p).is_prime():
        verify_y_greater_than(p, 10^1000)
