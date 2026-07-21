from math import comb
def Creal(x, a):           # generalized binomial C(x,a) for real x
    if a == 0: return 1.0
    if x < a - 1e-12: return 0.0
    p = 1.0
    for i in range(a): p *= (x - i)
    return p / __import__('math').factorial(a)
def invC(m, a):            # unique real x*>=a with C(x*,a)=m  (bisection)
    if m <= 0: return float(a) - 1.0
    lo, hi = float(a) - 1.0, float(a)
    while Creal(hi, a) < m: hi *= 2
    for _ in range(300):
        mid = (lo + hi) / 2
        if Creal(mid, a) < m: lo = mid
        else: hi = mid
    return hi

print("(a) monotonicity of C(x,a) on [a,inf):")
for a in [2,4,5]:
    xs = [a + 0.1*i for i in range(30)]
    mono = all(Creal(xs[i],a) < Creal(xs[i+1],a) for i in range(len(xs)-1))
    print(f"    a={a}: strictly increasing = {mono}")

print("\n(e) WORKED EXAMPLE r=3, c=9, kappa_(3,7)=16:")
r, c, s, m = 3, 9, 7, 16
a = s - r                       # 4
floor_s   = comb(c - r, a)      # C(6,4)
floor_s1  = comb(c - r, a + 1)  # C(6,5)
print(f"    a={a}  floor_s=C({c-r},{a})={floor_s}  -> T3 fires? {m == floor_s}  (m={m})")
xs = invC(m, a)
print(f"    x* with C(x*,{a})={m}: x*={xs:.6f}   check C(x*,{a})={Creal(xs,a):.6f}")
ub = Creal(xs, a + 1)
print(f"    C(x*,{a+1})={ub:.6f}  floor={int(ub)}   floor_s1=C({c-r},{a+1})={floor_s1}")
print(f"    T3+ fires? {int(ub) == floor_s1}  => kappa_(3,8) certified = {floor_s1}")

print("\n(d) does T3+ SUBSUME T3?  (m == floor_s should give x* == c-r exactly)")
for (rr,cc,ss) in [(3,9,7),(2,8,5),(4,12,9),(1,6,3)]:
    aa = ss - rr; fl = comb(cc-rr, aa)
    x = invC(fl, aa)
    ok = abs(x - (cc-rr)) < 1e-6 and int(Creal(x, aa+1)+1e-9) == comb(cc-rr, aa+1)
    print(f"    r={rr} c={cc} s={ss}: m=floor={fl} -> x*={x:.6f} (c-r={cc-rr})  T3+ reproduces T3: {ok}")

print("\n(f) SOUNDNESS SWEEP: over many (r,c,s,m), does T3+ ever fire with an upper bound")
print("    BELOW the T1 floor (which would be a contradiction => unsound)?")
bad = 0; fired = 0; tested = 0
for r in range(1,6):
  for c in range(r+2, 16):
    for s in range(r+1, min(c, 12)):
      a = s - r
      if c < s + 1: continue
      fl_s, fl_s1 = comb(c-r, a), comb(c-r, a+1)
      for m in range(fl_s, fl_s + 40):     # kappa >= T1 floor
        tested += 1
        x = invC(m, a); ub = int(Creal(x, a+1) + 1e-9)
        if ub == fl_s1:
            fired += 1
            if ub < fl_s1: bad += 1        # upper < lower => impossible
        if ub < fl_s1: bad += 1            # ANY case where derived upper < T1 lower
print(f"    tested={tested}  T3+ fired={fired}  UNSOUND cases (upper < T1 lower)={bad}")
