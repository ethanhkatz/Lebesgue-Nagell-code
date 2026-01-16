for p in [3, 5, 7, 11, 13]:
    for r in range(-(p-1)//2, (p-1)//2+1):
        f = (1/(2*sqrt(2)) * ((1+sqrt(2))^r*(x+sqrt(2))^p-(1-sqrt(2))^r*(x-sqrt(2))^p)).simplify_rational()
        S = gp.thueinit(f, flag=1) # flag=1 means we get unconditional results. See the GP/PARI documentation
        print(f"p={p}, r={r}: {gp.thue(S, 1)}")