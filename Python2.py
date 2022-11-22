import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt


num = np.array([2,3,4,5,6])#, 7])

M = pa.mat()

M2 = pa.mat()
M3 = pa.mat()
M4 = pa.mat()
M5 = pa.mat()
M6 = pa.mat()
#M7 = pa.mat()

M2.load("FileMatin2T1.0L20ord0.bin")
M3.load("FileMatin3T1.0L20ord0.bin")
M4.load("FileMatin4T1.0L20ord0.bin")
M5.load("FileMatin5T1.0L20ord0.bin")
M6.load("FileMatin6T1.0L20ord0.bin")
#M7.load("FileMatin7T1.0L20ord0.bin")

M = pa.join_rows(M2, M3, M4, M5)#, M6)
M = pa.join_rows(M, M6)


MeanEpsilon = M[0,:]
MeanAbs_m = M[2,:]

num = np.array([2,3,4,5,6])#, 7])

plt.figure()
plt.subplot(241)
plt.title("$T=1.0 J/k_B$, Unordered initial state, $<\epsilon>$")
plt.plot(num, MeanEpsilon, 'o')
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<\epsilon>$')
plt.subplot(242)
plt.plot(num, MeanAbs_m, 'o')
plt.title("$T=1.0 J/k_B$, Unordered initial state, $<|m|>$")
plt.xlabel("log10(NumCycles)"); plt.ylabel("$<|m|>$")





M2b = pa.mat()
M3b = pa.mat()
M4b = pa.mat()
M5b = pa.mat()
M6b = pa.mat()
#M7b = pa.mat()

M2b.load("FileMatin2T2.4L20ord0.bin")
M3b.load("FileMatin3T2.4L20ord0.bin")
M4b.load("FileMatin4T2.4L20ord0.bin")
M5b.load("FileMatin5T2.4L20ord0.bin")
M6b.load("FileMatin6T2.4L20ord0.bin")
#M7b.load("FileMatin7T2.4L20ord0.bin")

Mb = pa.join_rows(M2b,M3b, M4b, M5b)#, M6b)
Mb = pa.join_rows(Mb, M6b)

MeanEpsilonb = Mb[0,:]
MeanAbs_mb = Mb[2,:]



plt.subplot(243)
plt.title("$T=2.4 J/k_B$, Unordered initial state, $<\epsilon>$")
plt.plot(num, MeanEpsilonb, 'o')
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<\epsilon>$')
plt.subplot(244)
plt.title("$T=2.4 J/k_B$, Unordered initial state, $<|m|>$")
plt.plot(num, MeanAbs_mb, 'o')
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<|m|>$')



M2c = pa.mat()
M3c = pa.mat()
M4c = pa.mat()
M5c = pa.mat()
M6c = pa.mat()
#M7c = pa.mat()

M2c.load("FileMatin2T1.0L20ord1.bin")
M3c.load("FileMatin3T1.0L20ord1.bin")
M4c.load("FileMatin4T1.0L20ord1.bin")
M5c.load("FileMatin5T1.0L20ord1.bin")
M6c.load("FileMatin6T1.0L20ord1.bin")
#M7c.load("FileMatin7T1.0L20ord1.bin")

Mc = pa.join_rows(M2c, M3c, M4c, M5c)
Mc = pa.join_rows(Mc, M6c)

MeanEpsilonc = Mc[0,:]
MeanAbs_mc = Mc[2,:]


plt.subplot(245)
plt.title("$T=1.0 J/k_B$, Ordered initial state, $<\epsilon>$")
plt.plot(num, MeanEpsilonc, 'o')
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<\epsilon>$')
plt.subplot(246)
plt.title("$T=1.0 J/k_B$, Ordered initial state, $<|m|>$")
plt.plot(num, MeanAbs_mc, 'o')
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<|m|>$')

M2d = pa.mat()
M3d = pa.mat()
M4d = pa.mat()
M5d = pa.mat()
M6d = pa.mat()
#M7d = pa.mat()

M2d.load("FileMatin2T2.4L20ord1.bin")
M3d.load("FileMatin3T2.4L20ord1.bin")
M4d.load("FileMatin4T2.4L20ord1.bin")
M5d.load("FileMatin5T2.4L20ord1.bin")
M6d.load("FileMatin6T2.4L20ord1.bin")
#M7d.load("FileMatin7T2.4L20ord1.bin")

Md = pa.join_rows(M2d, M3d, M4d, M5d)
Md = pa.join_rows(Md, M6d)

MeanEpsilond = Md[0,:]
MeanAbs_md = Md[2,:]


plt.subplot(247)
plt.title("$T=2.4 J/k_B$, Ordered initial state, $<\epsilon>$")
plt.plot(num, MeanEpsilond, 'o')
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<\epsilon>$')
plt.subplot(248)
plt.title("$T=2.4 J/k_B$, Ordered initial state, $<|m|>$")
plt.plot(num, MeanAbs_md, 'o')
plt.xlabel('log10(NumCycles)'); plt.ylabel('$<|m|>$')


CV = M[4,:]
chi = M[5,:]
CVb = Mb[4,:]
chib = Mb[5,:]
CVc = Mc[4,:]
chic = Mc[5,:]
CVd = Md[4,:]
chid = Md[5,:]

plt.figure()
plt.subplot(241)
plt.title("$T=1.0 J/k_B$, Unordered initial state, $C_V$")
plt.plot(num, CV, 'o')
plt.xlabel("log10(NumCycles)"); plt.ylabel("$C_V$")
plt.subplot(242)
plt.title("$T=1.0 J/k_B$, Unordered initial state, $\chi$")
plt.plot(num, chi, 'o')
plt.xlabel("log10(NumCycles)"); plt.ylabel("$\chi$")
plt.subplot(243)
plt.title("$T=2.4 J/k_B$, Unordered initial state, $C_V$")
plt.plot(num, CVb, 'o')
plt.xlabel("log10(NumCycles)"); plt.ylabel("$C_V$")
plt.subplot(244)
plt.title("$T=2.4 J/k_B$, Unordered initial state, $\chi$")
plt.plot(num, chib, 'o')
plt.xlabel("log10(NumCycles)"); plt.ylabel("$\chi$")
plt.subplot(245)
plt.title("$T=1.0 J/k_B$, Ordered initial state, $C_V$")
plt.plot(num, CVc, 'o')
plt.xlabel("log10(NumCycles)"); plt.ylabel("$C_V$")
plt.subplot(246)
plt.title("$T=1.0 J/k_B$, Ordered initial state, $\chi$")
plt.plot(num, chic, 'o')
plt.xlabel("log10(NumCycles)"); plt.ylabel("$\chi$")
plt.subplot(247)
plt.title("$T=2.4 J/k_B$, Ordered initial state, $C_V$")
plt.plot(num, CVd, 'o')
plt.xlabel("log10(NumCycles)"); plt.ylabel("$C_V$")

plt.title("$T=2.4 J/k_B$, Unordered initial state, $C_V$")
plt.subplot(248)
plt.title("$T=2.4 J/k_B$, Ordered initial state, $\chi$")
plt.plot(num, chid, 'o')
plt.xlabel("log10(NumCycles)"); plt.ylabel("$\chi$")