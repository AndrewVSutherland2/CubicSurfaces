\\ pentareduce.gp -- (K5, lambda)-model reduction for cubic surfaces.
\\
\\ A smooth cubic surface with pentahedral quintic pol5 (separable, sigma_5 != 0)
\\ is Q-isomorphic to  { Tr_{A/Q}(lambda * y^3) = 0 } on the hyperplane
\\ { Tr_{A/Q}(y) = 0 } of A = Q[T]/(pol5), with lambda the class of T.  The gauge
\\ y -> c*y (c in A^*) replaces (lambda, hyperplane) by (lambda*c^3, Tr(c*y)=0)
\\ and is a GL_4(Q) change of coordinates: same surface over Q.  This program
\\ minimizes over that gauge: polredbest each factor of pol5 (using the KNOWN
\\ prime support -- the huge discriminant of pol5 is index junk), strips the
\\ cube part of the ideal (lambda) at the known primes, and LLL-reduces the
\\ trace-zero lattice under a |lambda|^(1/3)-weighted metric.  No integer larger
\\ than the invariants E, Delta_Cl is ever factored (and that happens upstream).
\\
\\ Run with cwd = a scratch dir containing pentared_in.gp, which must define
\\   PC : 6 integer coefficients of the monic pentahedral quintic, constant first
\\   SP : vector of (probable) primes -- known support of E, Delta_Cl and scaling
\\ Writes pentared_out.txt:
\\   line 1: "OK <maxdigits>"  or  "FAIL <message>"
\\   line 2: 20 integer coefficients of the reduced cubic, space-separated, in
\\           the monomial order x_i x_j x_k for 1 <= i <= j <= k <= 4
\\   line 3: diagnostics
default(parisize, "1G");
default(parisizemax, "8G");
OUT = "pentared_out.txt";
read("pentared_in.gp");

\\ real embedding vector of an element given by coordinates on nf.zk
embvec(nf, x) =
{
  my(cv = nfeltembed(nf, x), r1 = nf.r1, r2 = nf.r2, out = vector(r1 + 2*r2));
  for (i = 1, r1, out[i] = real(cv[i]));
  for (i = 1, r2,
    out[r1+2*i-1] = sqrt(2)*real(cv[r1+i]);
    out[r1+2*i]   = sqrt(2)*imag(cv[r1+i]));
  out;
}

\\ weight block |sigma(lam)|^(1/3) matching embvec's row layout
wtvec(nf, lam) =
{
  my(cv = nfeltembed(nf, lam), r1 = nf.r1, r2 = nf.r2, out = vector(r1 + 2*r2));
  for (i = 1, r1, out[i] = abs(real(cv[i]))^(1/3));
  for (i = 1, r2,
    out[r1+2*i-1] = abs(cv[r1+i])^(1/3);
    out[r1+2*i]   = abs(cv[r1+i])^(1/3));
  out;
}

main() =
{
  my(pol5 = Polrev(PC, 'T));
  if (poldegree(pol5) != 5, error("pol5 has degree != 5"));
  my(hd = vecmax(apply(x -> #digits(abs(x)+1), PC)));
  default(realprecision, 2*hd + 200);
  my(SPv = Set(SP));
  my(big = select(x -> x > 10^15, Vec(SPv)));
  if (#big, addprimes(big));
  my(fa = factor(pol5));
  if (vecmax(fa[,2]) > 1, error("pol5 not separable"));
  my(nfac = matsize(fa)[1]);
  my(FD = vector(nfac), FNF = vector(nfac), FLAM = vector(nfac), FCAND = vector(nfac));

  for (j = 1, nfac,
    my(P = fa[j,1], d = poldegree(P));
    FD[j] = d;
    if (d == 1,
      my(r0 = -polcoef(P, 0), c = 1);
      FNF[j] = 0; FLAM[j] = r0;
      foreach (SPv, p, my(v = valuation(r0, p)); if (v >= 3, c /= p^(v\3)));
      FCAND[j] = if (c == 1, [1], [c, 1]),
    \\ else: a genuine field factor
      my(vv = polredbest([P, SPv], 1));
      my(nf = nfinit([vv[1], SPv]));
      my(lamb = nfalgtobasis(nf, lift(vv[2])));
      FNF[j] = nf; FLAM[j] = lamb;
      my(J = matid(d));
      foreach (SPv, p,
        my(dec = idealprimedec(nf, p));
        foreach (dec, pr,
          my(e = idealval(nf, lamb, pr));
          if (e >= 3, J = idealmul(nf, J, idealpow(nf, pr, e\3)))));
      my(cands = [vectorv(d, i, i==1)]);          \\ c = 1
      if (J != matid(d),
        my(Bi = idealhnf(nf, idealinv(nf, J)));
        my(Rm = Mat(vector(d, t, Col(embvec(nf, Bi[,t])))));
        my(Bu = Bi * qflll(Rm));
        cands = concat([[Bu[,1]], [Bu[,2]], cands]));
      FCAND[j] = cands));

  \\ block offsets into the concatenated 5-dim coordinate space
  my(offs = vector(nfac+1)); offs[1] = 0;
  for (j = 1, nfac, offs[j+1] = offs[j] + FD[j]);
  if (offs[nfac+1] != 5, error("factor degrees do not sum to 5"));

  \\ fixed embedding block matrix (5 x 5 real, block diagonal per factor)
  my(E5 = matrix(5, 5));
  for (j = 1, nfac,
    if (FD[j] == 1,
      E5[offs[j]+1, offs[j]+1] = 1.0,
      for (t = 1, FD[j],
        my(ev = embvec(FNF[j], vectorv(FD[j], i, i==t)));
        for (r = 1, FD[j], E5[offs[j]+r, offs[j]+t] = ev[r]))));

  \\ enumerate gauge combos (mixed radix over per-factor candidate lists)
  my(ncomb = 1); for (j = 1, nfac, ncomb *= #FCAND[j]);
  my(tried = 0, bestdig = oo, bestcf = 0);
  for (ci = 0, min(ncomb, 48) - 1,
    my(t = ci, cs = vector(nfac));
    for (j = 1, nfac, cs[j] = FCAND[j][1 + (t % #FCAND[j])]; t \= #FCAND[j]);
    \\ lambda' and the coupled trace row / weights
    my(lamP = vector(nfac), trrow = vector(5), wts = vector(5), ok = 1);
    for (j = 1, nfac,
      if (FD[j] == 1,
        lamP[j] = FLAM[j] * cs[j]^3;
        if (lamP[j] == 0, ok = 0; break);
        trrow[offs[j]+1] = cs[j];
        wts[offs[j]+1] = abs(lamP[j])^(1/3),
      \\ else
        my(nf = FNF[j]);
        lamP[j] = nfeltmul(nf, FLAM[j], nfeltpow(nf, cs[j], 3));
        for (t2 = 1, FD[j],
          trrow[offs[j]+t2] = nfelttrace(nf,
            nfeltmul(nf, cs[j], vectorv(FD[j], i, i==t2))));
        my(wb = wtvec(nf, lamP[j]));
        for (r = 1, FD[j], wts[offs[j]+r] = wb[r])));
    if (!ok, next);
    my(den = lcm(apply(denominator, trrow)));
    my(K = matkerint(Mat(vector(5, t2, trrow[t2]*den))));
    if (matsize(K)[2] != 4, next);
    tried++;
    \\ two lattice metrics: lambda-weighted and plain
    for (wv = 0, 1,
      my(W = if (wv, matdiagonal(wts), matid(5)));
      my(Bred = K * qflll(W * (E5 * K)));
      \\ coefficients of Tr(lambda' * (sum x_i b_i)^3)
      my(cf = vector(20), idx = 0, bad = 0);
      for (i2 = 1, 4, for (j2 = i2, 4, for (k2 = j2, 4,
        idx++;
        my(val = 0);
        for (j = 1, nfac,
          my(o = offs[j]);
          if (FD[j] == 1,
            val += lamP[j] * Bred[o+1, i2] * Bred[o+1, j2] * Bred[o+1, k2],
          \\ else
            my(nf = FNF[j]);
            my(b1 = Bred[o+1..o+FD[j], i2], b2 = Bred[o+1..o+FD[j], j2],
               b3 = Bred[o+1..o+FD[j], k2]);
            my(pr = nfeltmul(nf, nfeltmul(nf, b1, b2), b3));
            val += nfelttrace(nf, nfeltmul(nf, lamP[j], pr))));
        my(mm = if (i2==j2 && j2==k2, 1, if (i2==j2 || j2==k2, 3, 6)));
        cf[idx] = mm * val)));
      my(d2 = lcm(apply(denominator, cf)));
      my(cfz = cf * d2);
      my(g = gcd(cfz));
      if (g == 0, next);
      cfz /= g;
      my(dig = vecmax(apply(x -> #digits(abs(x)+1), cfz)));
      if (dig < bestdig, bestdig = dig; bestcf = cfz)));
  if (bestdig == oo, error("no gauge combo produced a nonzero model"));
  [bestcf, bestdig,
   Str("factors=", FD, " combos_tried=", tried, " prec=", default(realprecision))];
}

{
  my(r = iferr(main(), ERR, [0, 0, Str(ERR)]));
  if (r[1] != 0,
    write(OUT, Str("OK ", r[2]));
    write(OUT, concat(vector(#r[1], i, Str(r[1][i], " "))));
    write(OUT, r[3]),
  \\ else
    write(OUT, Str("FAIL ", r[3])));
}
quit
