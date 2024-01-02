import numpy as np
import scipy as scp
from icecream import ic

# specify the basis set
alpha1s = np.array([0.075139, 0.231031, 0.994203])
coeff1s = np.array([0.700115, 0.399513, -0.09997])
alpha2s = np.array([0.109818, 0.405771, 2.227660])
coeff2s = np.array([0.444635, 0.535328, 0.154329])

def optim(a1, c1, a2, c2, s1, s2):
  a1 = a1 * s1**2
  a2 = a2 * s2**2
  a = np.concatenate((a1, a2))
  c = np.concatenate((c1, c2))
  S = np.zeros((len(a), len(a)))
  for i in range(len(a)):
    for j in range(len(a)):
      Norm = (2*a[i]/np.pi)**(3/4) * (2*a[j]/np.pi)**(3/4)
      S[i][j] = Norm * (np.pi/(a[i]+a[j]))**(3/2) * c[i] * c[j]
      
       #* np.exp(-a[i]*a[j]/(a[i]+a[j])**2)
  T = np.zeros((len(a), len(a)))
  for i in range(len(a)):
    for j in range(len(a)):
      Norm = (2*a[i]/np.pi)**(3/4) * (2*a[j]/np.pi)**(3/4)
      T[i][j] += 3 * a[j] * S[i,j]
      T[i][j] -= 2 * a[j] * a[j] * S[i,j] * (0.5/(a[i]+a[j])) 
      T[i][j] -= 2 * a[j] * a[j] * S[i,j] * (0.5/(a[i]+a[j])) 
      T[i][j] -= 2 * a[j] * a[j] * S[i,j] * (0.5/(a[i]+a[j])) 

  V = np.zeros((len(a), len(a)))
  for i in range(len(a)):
    for j in range(len(a)):
      Norm = (2*a[i]/np.pi)**(3/4) * (2*a[j]/np.pi)**(3/4)
      V[i,j] = Norm * -2 * np.pi / (a[i]+a[j]) * c[i] * c[j]
  
  H = T + V
  E, wf = scp.linalg.eigh(H, S)
  kinetic = wf[:, 1] @ T @ wf[:, 1].T
  potential = wf[:, 1] @ V @ wf[:, 1].T
  e = kinetic + potential
  return e, kinetic, potential

result = optim(alpha1s, coeff1s, alpha2s, coeff2s, 0.5, 1)
ic(result)