import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt

M3 = pa.mat()
M4 = pa.mat()
M5 = pa.mat()
M6 = pa.mat()
M7 = pa.mat()

M3.load("FileMatin3T1.0L2ord0.bin")
M4.load("FileMatin4T1.0L2ord0.bin")
M5.load("FileMatin5T1.0L2ord0.bin")
M6.load("FileMatin6T1.0L2ord0.bin")
#M7.load("FileMat7.bin")

M = pa.join_rows(M3, M4, M5, M6) # Can only join 4 at a time
#M = pa.join_rows(M, M7)

MeanEpsilon = M[0,:]
MeanEpsilonSquared = M[1,:]
MeanAbs_m = M[2,:]
Mean_mSqaured = M[3,:]
C_V = M[4,:]
chi = M[5,:]

num = np.array([3,4,5,6])#,7])
plt.figure()
plt.subplot(611)
plt.plot(num, MeanEpsilon)
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<\epsilon>$')
plt.subplot(612)
plt.plot(num, MeanEpsilonSquared)
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<\epsilon^2>$')
plt.subplot(613)
plt.plot(num, MeanAbs_m)
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<|m|>$')
plt.subplot(614)
plt.plot(num, Mean_mSqaured)
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<m^2>$')
plt.subplot(615)
plt.plot(num, C_V)
plt.xlabel('log10(NumCycles)'); plt.ylabel('$C_V$')
plt.subplot(616)
plt.plot(num, chi)
plt.xlabel('log10(NumCycles)'); plt.ylabel('$\chi$')
