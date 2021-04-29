import numpy as np
import matplotlib.pyplot as plt

beta = 0.4
gamma = 1/10
b = 150
n = 2000
h = b/n
t = np.linspace( 0, n*h, n )

S = np.zeros(n)
I = np.zeros(n)
R = np.zeros(n)
S[0] = 999
I[0] = 1
R[0] = 0
N = S[0] + I[0] + R[0]

def dSdt(S , I, R):
    return -(beta * S * I)/N

def dIdt(S, I, R):
    return (beta * S * I)/N - gamma * I

def dRdt(S, I, R):
    return gamma * I

for j in range(0, n-1, 1):
    k0 = h * (dSdt(S[j], I[j], R[j]))
    l0 = h * (dIdt(S[j], I[j], R[j]))
    o0 = h * (dRdt(S[j], I[j], R[j]))

    k1 = h * dSdt(S[j]+(1/2)*k0, I[j] + (1/2)*l0, R[j] + (1/2)*o0)
    l1 = h * dIdt(S[j]+(1/2)*k0, I[j] + (1/2)*l0, R[j] + (1/2)*o0)
    o1 = h * dRdt(S[j]+(1/2)*k0, I[j] + (1/2)*l0, R[j] + (1/2)*o0)

    k2 = h * dSdt(S[j] + (1 / 2) * k1, I[j] + (1 / 2) * l1, R[j] + (1 / 2) * o1)
    l2 = h * dIdt(S[j] + (1 / 2) * k1, I[j] + (1 / 2) * l1, R[j] + (1 / 2) * o1)
    o2 = h * dRdt(S[j] + (1 / 2) * k1, I[j] + (1 / 2) * l1, R[j] + (1 / 2) * o1)

    k3 = h * dSdt(S[j] + k2, I[j] + l2, R[j] + o2)
    l3 = h * dIdt(S[j] + k2, I[j] + l2, R[j] + o2)
    o3 = h * dRdt(S[j] + k2, I[j] + l2, R[j] + o2)

    k = (k0 + 2*k1 + 2*k2 + k3) * (1/6)
    l = (l0 + 2*l1 + 2*l2 + l3) * (1/6)
    o = (o0 + 2*o1 + 2*o2 + o3) * (1/6)

    S[j+1] = S[j] + k
    I[j+1] = I[j] + l
    R[j+1] = R[j] + o

fig = plt.figure(figsize= [6,4])
plt.plot(t, S, label = "S(t)")
plt.plot(t, I, label = "I(t)")
plt.plot(t, R, label = "R(t)")
plt.legend(loc = [0.7,0.5])
plt.grid()
plt.show()