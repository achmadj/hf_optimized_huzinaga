import numpy as np

zeta = 1
a = np.array([0.0247, 0.0798, 0.3371]) * zeta**2
c = np.array([1, 1, 1]) * 128**(1/4) * a**(5/4) / np.pi**(3/4)

dsquare = np.zeros((3, 3))
aplus = np.zeros((3, 3))
smatrix = np.zeros((3, 3))
ionic = np.zeros((3, 3))
ll1pr2 = np.zeros((3, 3))
tij = np.zeros((3, 3))
hmatrix = np.zeros((3, 3))

r = np.linspace(0, 20, 3000)
weight = np.zeros(3000)
weight[::2] = 2
weight[1::2] = 4
weight[0] = 1
weight[-1] = 1

def kinetic(alphai, alphaj):
    ij = r**2 * np.exp(-alphai * r**2)
    grad2 = 4 * alphaj**2 * r**4 - 10 * alphaj * r**2 + 2
    return np.sum(weight * ij * grad2) * (r[1] - r[0]) * 4 * np.pi / 9

def angular(alphai):
    ij = r**2 * np.exp(-alphai * r**2)
    return np.sum(weight * ij) * (r[1] - r[0]) * 4 * np.pi / 9

def potential(alphai):
    ij = r**2 * np.exp(-alphai * r**2)
    return np.sum(r * weight * ij) * (r[1] - r[0]) * 4 * np.pi / 9

def overlap(alphai):
    ij = r**2 * np.exp(-alphai * r**2)
    return np.sum(r**2 * weight * ij) * (r[1] - r[0]) * 4 * np.pi / 9

for i in range(3):
    for j in range(i+1):
        dsquare[i,j] = c[i] * c[j]
        dsquare[j,i] = dsquare[i,j]
        aplus[i,j] = a[i] + a[j]
        aplus[j,i] = aplus[i,j]
        smatrix[i,j] = overlap(a[i] + a[j]) * c[i] * c[j]
        ionic[i,j] = -potential(a[i] + a[j]) * c[i] * c[j]
        ll1pr2[i,j] = angular(a[i] + a[j]) * c[i] * c[j]
        tij[i,j] = kinetic(a[i] + a[j], a[j]) * (-0.5) * c[i] * c[j]
        smatrix[j,i] = smatrix[i,j]
        ionic[j,i] = ionic[i,j]
        ll1pr2[j,i] = ll1pr2[i,j]
        tij[j,i] = tij[i,j]
        hmatrix[i,j] = kinetic(a[i] + a[j], a[j]) * (-0.5) * c[i] * c[j] - potential(a[i] + a[j]) * c[i] * c[j] + angular(a[i] + a[j]) * c[i] * c[j]
        hmatrix[j,i] = hmatrix[i,j]

ij_tot = np.sum(smatrix) - np.trace(smatrix)
v_tot = np.sum(ionic) - np.trace(ionic)
t_tot = np.sum(tij) - np.trace(tij)
ll1pr2_tot = np.sum(ll1pr2) - np.trace(ll1pr2)

Etot = v_tot + t_tot + ll1pr2_tot
print(Etot)

E, wf = np.linalg.eig(np.linalg.inv(smatrix) @ hmatrix)
print(E, wf)