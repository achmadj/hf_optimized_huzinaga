import math
import numpy as np
from scipy import special 

def buildG(basis, G, K):
  Ntei = 0 # create a two-electron integral counter, as TEIs can take a long time
  for A, bA in enumerate(basis): # retrieve atomic orbital A from basis
    for B, bB in enumerate(basis):
      for C, bC in enumerate(basis):
        for D, bD in enumerate(basis):
  
          Ntei +=1
          if Ntei % 250 == 0:
            print ('Computed '+ str(Ntei) + ' of ' + str(K**4) + ' integrals.')
  
          for a, dA in zip(bA['a'],bA['d']): # retrieve alpha and contract coefficients for atomic orbital A
            for b, dB in zip(bB['a'],bB['d']):
              for c, dC in zip(bC['a'],bC['d']):
                for d, dD in zip(bD['a'],bD['d']):
   
                  RA, RB, RC, RD = bA['R'], bB['R'], bC['R'], bD['R'] # hold coordinates of AOs A, B, C, and D
  
                  ## variables for angular momenta, in terms of x,y,z, for each orbital
                  lA, mA, nA = bA['l'], bA['m'], bA['n']
                  lB, mB, nB = bB['l'], bB['m'], bB['n']
                  lC, mC, nC = bC['l'], bC['m'], bC['n']
                  lD, mD, nD = bD['l'], bD['m'], bD['n']
  
                  tei  = dA * dB * dC * dD # multiply together the contraction coefficients
                  tei *= Gxyz(lA,mA,nA,lB,mB,nB,lC,mC,nC,lD,mD,nD,a,b,c,d,RA,RB,RC,RD) # multiply by integral over primitives
   
                  G[A,B,C,D] += tei

  return G

def Gxyz(lA,mA,nA,lB,mB,nB,lC,mC,nC,lD,mD,nD,a,b,c,d,RA,RB,RC,RD):
  gP = a + b
  gQ = c + d

  delta = 1/(4*gP) + 1/(4*gQ)

  RP = gaussianProduct(a,RA,b,RB,gP)
  RQ = gaussianProduct(c,RC,d,RD,gQ)

  ABsq = IJsq(RA,RB)
  CDsq = IJsq(RC,RD)
  PQsq = IJsq(RP,RQ)

  Gxyz = 0
  for l in range(0,lA+lB+1):
    for r in range(0,int(l/2)+1):
      for lp in range(0,lC+lD+1):
        for rp in range(0,int(lp/2)+1):
          for i in range(0,int((l+lp-2*r-2*rp)/2)+1):
            gx = gi(l,lp,r,rp,i, lA,lB,RA[0],RB[0],RP[0],gP, lC,lD,RC[0],RD[0],RQ[0],gQ )

            for m in range(0,mA+mB+1):
              for s in range(0,int(m/2)+1):
                for mp in range(0,mC+mD+1):
                  for sp in range(0,int(mp/2)+1):
                    for j in range(0,int((m+mp-2*s-2*sp)/2)+1):
                      gy = gi(m,mp,s,sp,j, mA,mB,RA[1],RB[1],RP[1],gP, mC,mD,RC[1],RD[1],RQ[1],gQ)

                      for n in range(0,nA+nB+1):
                        for t in range(0,int(n/2)+1):
                          for np in range(0,nC+nD+1):
                            for tp in range(0,int(np/2)+1):
                              for k in range(0,int((n+np-2*t-2*tp)/2)+1):
                                gz = gi(n,np,t,tp,k, nA,nB,RA[2],RB[2],RP[2],gP, nC,nD,RC[2],RD[2],RQ[2],gQ)

                                nu    = l+lp+m+mp+n+np-2*(r+rp+s+sp+t+tp)-(i+j+k)
                                F     = BoysFunction(nu, PQsq/(4*delta))
                                Gxyz += gx * gy * gz * F

  Gxyz *= ( 2 * math.pi**2 ) / ( gP * gQ ) 
  Gxyz *= math.sqrt( math.pi / ( gP + gQ ) )
  Gxyz *= math.exp( -(a*b*ABsq)/gP ) 
  Gxyz *= math.exp( -(c*d*CDsq)/gQ )

  Na = N(a,lA,mA,nA)
  Nb = N(b,lB,mB,nB)
  Nc = N(c,lC,mC,nC)
  Nd = N(d,lD,mD,nD)

  Gxyz *= Na * Nb * Nc * Nd

  return Gxyz

def gi(l,lp,r,rp,i, lA,lB,Ai,Bi,Pi,gP, lC,lD,Ci,Di,Qi,gQ):
  delta = 1/(4*gP) + 1/(4*gQ)

  gi  = (-1)**l 
  gi *= theta(l,lA,lB,Pi-Ai,Pi-Bi,r,gP) * theta(lp,lC,lD,Qi-Ci,Qi-Di,rp,gQ)
  gi *= (-1)**i * (2*delta)**(2*(r+rp))
  gi *= special.factorial(l+lp-2*r-2*rp,exact=True) * delta**i
  gi *= (Pi-Qi)**(l+lp-2*(r+rp+i))
  gi /= (4*delta)**(l+lp) * special.factorial(i,exact=True)
  gi /= special.factorial(l+lp-2*(r+rp+i),exact=True)

  return gi

def theta(l,lA,lB,PA,PB,r,g):
  theta  = ck(l,lA,lB,PA,PB) 
  theta *= special.factorial(l,exact=True) * g**(r-l) 
  theta /= special.factorial(r,exact=True) * special.factorial(l-2*r,exact=True) 

  return theta

def ck(j,l,m,a,b):
  coefficient = 0.0
 
  for k in range(0,l+1):
    for i in range(0,m+1):
      if i + k == j:
        coefficient += special.binom(l,k) * special.binom(m,i) * a**(l-k) * b**(m-i)

  return coefficient

def N(a,l,m,n):
  N = (2*a/math.pi)**(3/4) * (((8*a)**(l+m+n) * special.factorial(l,exact=True) * special.factorial(m,exact=True) * special.factorial(n,exact=True))/(special.factorial(2*l,exact=True) * special.factorial(2*m,exact=True) * special.factorial(2*n,exact=True)))**(1/2)
  return N

def BoysFunction(nu, x):
  """
  The analytical function coded herein was suggested by CCL forum; similar to
  result when evaluating the Boys Function integral in Mathematica. 
  Depends on gamma functions, which are easily computed with SciPy's 
  library of special functions.
  """
  if x < 1e-7:
    return (2*nu+1)**(-1) - x*(2*nu+3)**(-1) # (Handout 4, Eq. 17)
  else:
    return (1/2) * x**(-(nu+0.5)) * special.gamma(nu+0.5) * special.gammainc(nu+0.5,x) # (Handout 4, Eq. 16)
  
def gaussianProduct(a,RA,b,RB,g):
  P = []
  for i in range(3):
    P.append( (a*RA[i]+b*RB[i])/g )

  return P

def IJsq(RI,RJ):
  return sum( (RI[i]-RJ[i])**2 for i in (0,1,2) )
