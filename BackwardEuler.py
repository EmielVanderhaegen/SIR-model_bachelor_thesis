import numpy as np
from numpy import zeros, linspace
import matplotlib.pyplot as plt

beta = 0.4
gamma = 0.1

t = np.linspace(0, 600, 10)
h=t[1]-t[0]

S=np.zeros_like(t)
I=np.zeros_like(t)
R=np.zeros_like(t)

#Initial condition
S[0] = 0.9
I[0] = 0.1
R[0] = 0

for k in range(len(t)-1):
    S[k+1]=S[k]/(1+beta*I[k]*h)
    I[k+1]=I[k]/(1-(beta*S[k+1]-gamma)*h)
    R[k+1]=R[k]+(gamma*I[k+1])*h


fig, axs = plt.subplots(1)
axs.plot(t[0:5], S[0:5],'-', label='S(t)')
axs.plot(t[0:5], I[0:5],'-', label='I(t)')
axs.plot(t[0:5], R[0:5],'-', label='R(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
