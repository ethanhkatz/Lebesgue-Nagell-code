def verify_parameters(p, y, Kp, L, R1, R2, mu, rho): # verifies Propositions 5.26 and 5.27 with a given choice of parameters
    p, y, Kp, L, R1, R2, mu, rho = RBF(p), RBF(y), RBF(Kp), RBF(L), RBF(R1), RBF(R2), RBF(mu), RBF(rho)
    a2 = (rho + 1) * log(sqrt(RBF(2)) + 1)
    S1 = ceil(L / R1)
    R = R1 + R2 - 1
    C1 = S1 + 2 + (1 - 2 * R2) / min(R2, p)
    C2 = Kp * L / min(R2, p)
    C3 = RBF(0.5) * (2 * (R - 1) + (S1 + (1 - L - 2 * R2) / min(R2, p) + 1) * p)
    C4 = (L * p) / (2 * min(R2, p))
    C5 = min(R2, p) / (12 * R)
    C6 = ((S1 + 2) * min(R2, p) + 1 - L - 2 * R2) / L
    C7 = (rho + 1) * (2 * log(sqrt(RBF(2)) + 1) + RBF(10)^(-50)) / p
    sigma = 1 - (1 - mu)**2 / 2
    K = ceil(Kp * log(y)) # lower bound for K
    if not (2 * R2 <= (K - 1) * L + 1):
        print(f"Verification failed for {Integer(p)}")
        return False
    log_bu = log(C3 / K + C4) + RBF(1.5) # upper bound for log(b)
    gu = RBF(0.25) - C5 + C5 * C6 / (C6 + K) # upper bound for g
    N = K * L # lower bound for N
    epsilonu = (3 * log(N) + log(2 * RBF.pi()) + 1 / (6 * N) + 2 * log(1 + (1 - RBF(-1).exp())**N)) / N # upper bound for epsilon
    C8 = Kp * (sigma * L - 1) * log(rho) - 2 * Kp * log_bu - gu * L * (2 * R + C2 * a2)
    C9 = 3 * log(L) + gu * L * (R * C7 + C1 * a2) + epsilonu
    if not (C3 > 0 and C6 > 0 and C8 > 0 and log(y) > 3 / C8 - 1 / Kp and N >= 3):
        print(f"Verification failed for {Integer(p)}")
        return False

    # proposition 5.26
    if not (C8 * log(y) - 3 * log(Kp * log(y) + 1) > C9):
        print(f"Verification failed for {Integer(p)}")
        return False

    # proposition 5.27
    if not (C1 > 0):
        print(f"Verification failed for {Integer(p)}")
        return False
    if not (log(rho) * mu * L * Kp - p / 2 + C2 / (C1 + C2 * log(y)) + RBF(7165)/10000 * L * C2 * y**(-p/2) < 0):
        print(f"Verification failed for {Integer(p)}")
        return False
    if (log(rho) * mu * L * Kp - p / 2) * log(y) + RBF(1053)/1000 + log(L * (C1 + C2 * log(y)) / 4) + RBF(7165)/10000 * L * (C1 + C2 * log(y)) * y**(-p/2) + log(rho) * mu * L <= 0:
        print(f"Verified {Integer(p)}")
        return True
    else:
        print(f"Verification failed for {Integer(p)}")
        return False

verify_parameters(919, RBF(10)^800, 26.64, 9, 1, 64, 0.58, 27.22)

verify_parameters(929, RBF(10)^300, 26.71, 9, 1, 64, 0.58, 27.8)

for p in [937, 941]:
    verify_parameters(p, RBF(10)^150, 26.76, 9, 1, 64, 0.58, 28.55)

verify_parameters(947, RBF(10)^100, 26.83, 9, 1, 64, 0.58, 29.3)

verify_parameters(953, RBF(10)^100, 26.42, 9, 1, 64, 0.59, 29.8)

for p in range(967, 998):
    if Integer(p).is_prime():
        verify_parameters(p, RBF(10)^100, 26.67, 9, 1, 64, 0.59, 29.8)

for p in range(1000, 1200):
    if Integer(p).is_prime():
        verify_parameters(p, RBF(10)^100, 27.04, 10, 1, 64, 0.57, 26.3)

for p in range(1200, 1952):
    if Integer(p).is_prime():
        verify_parameters(p, RBF(10)^50, 28.69, 10, 1, 69, 0.59, 33)
