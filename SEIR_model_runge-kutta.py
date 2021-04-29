import numpy as np
from numpy import zeros, linspace
import matplotlib.pyplot as plt

beta = 0.4
gamma = 0.1
mu = 0.2
b = 200
n = 10000
h = b/n
t = np.linspace(0, n*h, n)

S = zeros(n)
E = zeros(n)
I = zeros(n)
R = zeros(n)

#Initial condition
N = 1000
S[0] = 999
I[0] = 0
E[0] = 1
R[0] = 0
def dSdt(S ,E, I, R):
    return -(beta * S * I)/N

def dEdt(S, E, I, R):
    return beta*S*I/N - mu*E

def dIdt(S, E, I, R):
    return mu*E- gamma * I

def dRdt(S, E, I, R):
    return gamma * I

for j in range(0, n-1, 1):
    k0 = h * (dSdt(S[j], E[j], I[j], R[j]))
    m0 = h * (dEdt(S[j], E[j], I[j], R[j]))
    l0 = h * (dIdt(S[j], E[j], I[j], R[j]))
    o0 = h * (dRdt(S[j], E[j], I[j], R[j]))

    k1 = h * dSdt(S[j]+(1/2)*k0, E[j]+(1/2)*m0, I[j] + (1/2)*l0, R[j] + (1/2)*o0)
    m1 = h * dEdt(S[j]+(1/2)*k0, E[j]+(1/2)*m0, I[j] + (1/2)*l0, R[j] + (1/2)*o0)
    l1 = h * dIdt(S[j]+(1/2)*k0, E[j]+(1/2)*m0, I[j] + (1/2)*l0, R[j] + (1/2)*o0)
    o1 = h * dRdt(S[j]+(1/2)*k0, E[j]+(1/2)*m0, I[j] + (1/2)*l0, R[j] + (1/2)*o0)

    k2 = h * dSdt(S[j] + (1 / 2) * k1, E[j] + (1/2) * m1, I[j] + (1 / 2) * l1, R[j] + (1 / 2) * o1)
    m2 = h * dEdt(S[j] + (1 / 2) * k1, E[j] + (1/2) * m1, I[j] + (1 / 2) * l1, R[j] + (1 / 2) * o1)
    l2 = h * dIdt(S[j] + (1 / 2) * k1, E[j] + (1/2) * m1, I[j] + (1 / 2) * l1, R[j] + (1 / 2) * o1)
    o2 = h * dRdt(S[j] + (1 / 2) * k1, E[j] + (1/2) * m1, I[j] + (1 / 2) * l1, R[j] + (1 / 2) * o1)

    k3 = h * dSdt(S[j] + k2, E[j] + m2, I[j] + l2, R[j] + o2)
    m3 = h * dSdt(S[j] + k2, E[j] + m2, I[j] + l2, R[j] + o2)
    l3 = h * dIdt(S[j] + k2, E[j] + m2, I[j] + l2, R[j] + o2)
    o3 = h * dRdt(S[j] + k2, E[j] + m2, I[j] + l2, R[j] + o2)

    k = (k0 + 2*k1 + 2*k2 + k3) * (1/6)
    m = (m0 + 2*m1 + 2*m2 + m3) * (1/6)
    l = (l0 + 2*l1 + 2*l2 + l3) * (1/6)
    o = (o0 + 2*o1 + 2*o2 + o3) * (1/6)

    S[j+1] = S[j] + k
    E[j+1] = E[j] + m
    I[j+1] = I[j] + l
    R[j+1] = R[j] + o

fig = plt.figure(figsize= [6,4])
plt.plot(t, S, label = "S(t)")
plt.plot(t, E, label = "E(t)", color = "red")
plt.plot(t, I, label = "I(t)")
plt.plot(t, R, label = "R(t)")
plt.legend(loc = [0.7,0.5])
plt.grid()
plt.show()
