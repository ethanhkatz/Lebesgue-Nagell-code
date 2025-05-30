#Newf is the list of four newforms of weight 2, level 128, trivial character
Newf = Newforms(128)

#The newforms are all quadratic twists of each other, so just choose one of them and call it F
F = Newf[0]


#For part (1) of Proposition 8.2, we need to compute a_1093(F)
print(F.coefficient(1093))
# 30


#Check for part (2) of Proposition 8.2
#For each prime 5 <= p < 5000, check that |a_p(F)|^2 is not 1 mod p
P = prime_range(5,100)
for p in P:
    coeff = F.coefficient(p)
    if Integer(mod(coeff^2,p)) == 1:
        print(p)

#No primes are printed