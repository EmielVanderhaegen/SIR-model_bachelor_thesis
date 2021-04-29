import numpy as np
import matplotlib.pyplot as plt

beta = 0.4
gamma = 1/10
b = 150
n = 10000
h = b/n
t = np.linspace( 0, n*h, n )

S = np.zeros(n)
I = np.zeros(n)

S[0] = 999
I[0] = 1
N = S[0] + I[0]

def dSdt(S , I):
    return -(beta * S * I)/N + gamma*I

def dIdt(S, I):
    return (beta * S * I)/N - gamma * I


for j in range(0, n-1, 1):
    k0 = h * (dSdt(S[j], I[j]))
    l0 = h * (dIdt(S[j], I[j]))


    k1 = h * dSdt(S[j]+(1/2)*k0, I[j] + (1/2)*l0)
    l1 = h * dIdt(S[j]+(1/2)*k0, I[j] + (1/2)*l0)


    k2 = h * dSdt(S[j] + (1 / 2) * k1, I[j] + (1 / 2) * l1)
    l2 = h * dIdt(S[j] + (1 / 2) * k1, I[j] + (1 / 2) * l1)


    k3 = h * dSdt(S[j] + k2, I[j] + l2)
    l3 = h * dIdt(S[j] + k2, I[j] + l2)

    k = (k0 + 2*k1 + 2*k2 + k3) * (1/6)
    l = (l0 + 2*l1 + 2*l2 + l3) * (1/6)

    S[j+1] = S[j] + k
    I[j+1] = I[j] + l


fig = plt.figure(figsize= [6,4])
plt.plot(t, S, label = "S(t)")
plt.plot(t, I, label = "I(t)")
plt.legend(loc = [0.7,0.5])
plt.grid()
plt.show()
