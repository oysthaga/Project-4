import pyarma as pa
import matplotlib.pyplot as plt

T100 = pa.mat()
T100.load('Tk100.bin')
T10 = pa.mat()
T10.load('Tk10.bin')

eps40 = pa.mat()
m40 = pa.mat()
C_V40 = pa.mat()
chi40 = pa.mat()


eps40.load('epsilonin4L40ord1.bin')
m40.load('min4L40ord1.bin')
C_V40.load('C_Vin4L40ord1.bin')
chi40.load('chiin4L40ord1.bin')



eps60 = pa.mat()
m60 = pa.mat()
C_V60 = pa.mat()
chi60 = pa.mat()

eps60.load('epsilonin4L60ord1.bin')
m60.load('min4L60ord1.bin')
C_V60.load('C_Vin4L60ord1.bin')
chi60.load('chiin4L60ord1.bin')


eps80 = pa.mat()
m80 = pa.mat()
C_V80 = pa.mat()
chi80 = pa.mat()

eps80.load('epsilonin4L80ord1.bin')
m80.load('min4L80ord1.bin')
C_V80.load('C_Vin4L80ord1.bin')
chi80.load('chiin4L80ord1.bin')

eps100 = pa.mat()
m100 = pa.mat()
C_V100 = pa.mat()
chi100 = pa.mat()

eps100.load('epsilonin4L100ord1.bin')
m100.load('min4L100ord1.bin')
C_V100.load('C_Vin4L100ord1.bin')
chi100.load('chiin4L100ord1.bin')




plt.figure()
plt.subplot(411)
plt.plot(T100, eps40, label="$L=40$")
plt.plot(T100, eps60, label="$L=60$")
plt.plot(T10, eps80, label="$L=80$")
plt.plot(T10, eps100, label="$L=100$")
plt.xlabel("$T$ [$J/k_b$]"); plt.ylabel("$<\epsilon>$ [$J$]")
plt.legend()
plt.subplot(412)
plt.plot(T100, m40, label="$L=40$")
plt.plot(T100, m60, label="$L=60$")
plt.plot(T10, m80, label="$L=80$")
plt.plot(T10, m100, label="$L=100$")
plt.xlabel("$T$ [$J/k_b$]"); plt.ylabel("$<|m|>$")
plt.legend()

plt.subplot(413)
plt.plot(T100, C_V40, label="$L=40$")
plt.plot(T100, C_V60, label="$L=60$")
plt.plot(T10, C_V80, label="$L=80$")
plt.plot(T10, C_V100, label="$L=100$")
plt.xlabel("$T$ [$J/k_b$]"); plt.ylabel("$C_V$ [$J^2/k_b T^2$]")
plt.legend()
plt.subplot(414)
plt.plot(T100, chi40, label="$L=40$")
plt.plot(T100, chi60, label="$L=60$")
plt.plot(T10, chi80, label="$L=80$")
plt.plot(T10, chi100, label="$L=100$")
plt.xlabel("$T$ [$J/k_b$]"); plt.ylabel("$\chi$ [$N/k_b T$]")
plt.legend()



T100 = pa.mat()
T100.load('Tk100Tmin2.24Tmax2.32.bin')
T10 = pa.mat()
T10.load('Tk10Tmin2.24Tmax2.32.bin')

eps40 = pa.mat()
m40 = pa.mat()
C_V40 = pa.mat()
chi40 = pa.mat()


eps40.load('epsilonin4L40ord1Tmin2.24Tmax2.32.bin')
m40.load('min4L40ord1Tmin2.24Tmax2.32.bin')
C_V40.load('C_Vin4L40ord1Tmin2.24Tmax2.32.bin')
chi40.load('chiin4L40ord1Tmin2.24Tmax2.32.bin')



eps60 = pa.mat()
m60 = pa.mat()
C_V60 = pa.mat()
chi60 = pa.mat()

eps60.load('epsilonin4L60ord1Tmin2.24Tmax2.32.bin')
m60.load('min4L60ord1Tmin2.24Tmax2.32.bin')
C_V60.load('C_Vin4L60ord1Tmin2.24Tmax2.32.bin')
chi60.load('chiin4L60ord1Tmin2.24Tmax2.32.bin')


eps80 = pa.mat()
m80 = pa.mat()
C_V80 = pa.mat()
chi80 = pa.mat()

eps80.load('epsilonin4L80ord1Tmin2.24Tmax2.32.bin')
m80.load('min4L80ord1Tmin2.24Tmax2.32.bin')
C_V80.load('C_Vin4L80ord1Tmin2.24Tmax2.32.bin')
chi80.load('chiin4L80ord1Tmin2.24Tmax2.32.bin')

eps100 = pa.mat()
m100 = pa.mat()
C_V100 = pa.mat()
chi100 = pa.mat()

eps100.load('epsilonin4L100ord1Tmin2.24Tmax2.32.bin')
m100.load('min4L100ord1Tmin2.24Tmax2.32.bin')
C_V100.load('C_Vin4L100ord1Tmin2.24Tmax2.32.bin')
chi100.load('chiin4L100ord1Tmin2.24Tmax2.32.bin')




plt.figure()
plt.subplot(411)
plt.plot(T100, eps40, label="$L=40$")
plt.plot(T10, eps60, label="$L=60$")
plt.plot(T10, eps80, label="$L=80$")
plt.plot(T10, eps100, label="$L=100$")
plt.xlabel("$T$ [$J/k_b$]"); plt.ylabel("$<\epsilon>$ [$J$]")
plt.legend()

plt.subplot(412)
plt.plot(T100, m40, label="$L=40$")
plt.plot(T10, m60, label="$L=60$")
plt.plot(T10, m80, label="$L=80$")
plt.plot(T10, m100, label="$L=100$")
plt.xlabel("$T$ [$J/k_b$]"); plt.ylabel("$<|m|>$")
plt.legend()

plt.subplot(413)
plt.plot(T100, C_V40, label="$L=40$")
plt.plot(T10, C_V60, label="$L=60$")
plt.plot(T10, C_V80, label="$L=80$")
plt.plot(T10, C_V100, label="$L=100$")
plt.xlabel("$T$ [$J/k_b$]"); plt.ylabel("$C_V$ [$J^2/k_b T^2$]")
plt.legend()

plt.subplot(414)
plt.plot(T100, chi40, label="$L=40$")
plt.plot(T10, chi60, label="$L=60$")
plt.plot(T10, chi80, label="$L=80$")
plt.plot(T10, chi100, label="$L=100$")
plt.xlabel("$T$ [$J/k_b$]"); plt.ylabel("$\chi$ [$N/k_b T$]")
plt.legend()