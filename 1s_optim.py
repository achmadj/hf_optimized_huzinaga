import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as la

# specify the basis set
a = np.array([0.0822, 0.2247, 0.6733, 2.3465, 10.2465, 68.1600])
c1 = np.array([0.2426, 0.4922, 0.2943, 0.0928, 0.0194, 0.0026])
c = np.ones(6)

r = np.linspace(-3, 3, 200)
sto = np.sqrt(1/np.pi) * np.exp(-np.abs(r))
g = np.zeros((200, 6))
hmatrix = np.zeros((6, 6))
smatrix = np.zeros((6, 6))

for i in range(6):
    g[:,i] = (2*a[i]/np.pi)**(3/4) * np.exp(-a[i]*(r**2))
    for j in range(i+1):
        hmatrix[i,j] = (4*a[i]**3*a[j]**3)**(1/4)/(a[i]+a[j])*(6*a[j]/np.sqrt(a[i]+a[j])*(1-a[j]/(a[i]+a[j]))-4/np.sqrt(np.pi))
        hmatrix[j,i] = hmatrix[i,j]
        smatrix[i,j] = 2*np.sqrt(2)*(a[i]*a[j]/(a[i]+a[j])**2)**(3/4)
        smatrix[j,i] = smatrix[i,j]

E, wf = la.eigh(hmatrix, smatrix)
sto6g = g @ wf[:, 0].reshape(-1, 1)
# breakpoint()  
plt.plot(r, sto, r, sto6g, linewidth=2)
# plt.figure()

sto6g_huzinaga = g @ c1.reshape(-1, 1)
plt.plot(r, sto, r, sto6g_huzinaga, linewidth=2)
plt.show()
# print(np.column_stack((E[:,0], c1)))