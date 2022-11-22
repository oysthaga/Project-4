import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt

epsilon1 = pa.mat()
epsilon1.load("epsilonin6T1.0L20ord1.bin")

hist, bins = np.histogram(epsilon1, bins=100, density=True)
width = (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.figure()
plt.subplot(211)
plt.title("$T = 1.0 J/k_B$")
plt.bar(center, hist, align='center', width=width)
plt.xlabel("$\epsilon$")
plt.show()

epsilon2 = pa.mat()
epsilon2.load("epsilonin6T2.4L20ord1.bin")
hist, bins = np.histogram(epsilon2, bins=200, density=True)
width = (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.subplot(212)
plt.title("$T = 2.4 J/k_B$")
plt.bar(center, hist, align='center', width=width)
plt.xlabel("$\epsilon$")
plt.show()
