import numpy as np
from numpy import zeros, linspace
import matplotlib.pyplot as plt

beta = 0.4
gamma = 0.1
b= 150
n= 8000
h= b/n
t = np.linspace(0, n*h, n)

S = zeros(n)
I = zeros(n)
R = zeros(n)

#Initial condition
N = 1000
S[0] = 999
I[0] = 1
R[0] = 0

for k in range(n-1):
    S[k + 1] = S[k] - h * (beta * S[k] * I[k])/N
    I[k + 1] = I[k] + h * (beta * S[k] * I[k])/N - h * gamma * I[k]
    R[k + 1] = R[k] + h * gamma * I[k]


fig = plt.figure(figsize= [6,4])
plt.plot(t, S, label = "S(t)")
plt.plot(t, I, label = "I(t)")
plt.plot(t, R, label = "R(t)")
plt.legend(loc = [0.7,0.5])
plt.grid()
plt.show()


