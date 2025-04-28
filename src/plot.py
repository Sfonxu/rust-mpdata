import numpy as np
import matplotlib.pyplot as plt

data = np.load("output.npy")
data2 = np.load("initial.npy")
data3 = np.load("output_mpdata.npy")
xs, dx = np.linspace(-100, 300, 45, retstep=True)
print(dx)
plt.step(xs, data, where='mid', label="upwind")
plt.step(xs, data2, where='mid', label="initial")
plt.step(xs, data3, where='mid', label="mpdata")
plt.legend()
plt.show()
