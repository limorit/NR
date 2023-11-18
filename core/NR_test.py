import matplotlib.pyplot as plt
import numpy as np
from NR_model import NR, NR1, NR2, NR_rios, NRbyK


# tau = 167.5 * 1e3
# Kth = 1.8 * 1e3
# tau_FL = 80 * 1e3
# Scomp = 460 * 1e3
# tau_0 = tau / 12
# A = 260
# m = 1.78
# D = 6e-6
# a0 = 162e-6

# Smax = tau
# Sfl = tau_FL
# S2 = Scomp

# tau_Scomp = tau / Scomp
# tau_tau0 = tau / tau_0
# taufl_Scomp = tau_FL / Scomp
# tau_Scomp1 = 0.6
# tau_tau01 = 12
# taufl_Scomp1 = 0.5

# tau_Scomps = [0.1, 0.2, 0.3, 0.4, 0.5]
# tau_tau0s = [5, 6, 7, 8, 9, 10]
# taufl_Scomps = [ 0.5,0.6,0.7,0.8,0.9,1]

# for i in range(len(tau_Scomps)):
#     x, y = NR(tau_Scomps[i], tau_tau01, taufl_Scomp1, A=1, m=1)
#     plt.loglog(x, y, label="tau_Scomp = " + str(tau_Scomps[i]))
# plt.legend()
# plt.show()

# for i in range(len(tau_tau0s)):
#     x, y = NR(tau_Scomp1, tau_tau0s[i], taufl_Scomp1, A=1, m=1)
#     plt.loglog(x, y, label="tau_tau0 = " + str(tau_tau0s[i]))
# plt.legend()
# plt.show()

# for i in range(len(taufl_Scomps)):
#     x, y = NR(tau_Scomp1, tau_tau01, taufl_Scomps[i], A=1, m=1)
#     plt.loglog(x, y, label="taufl_Scomp = " + str(taufl_Scomps[i]))
# plt.legend()
# plt.show()

# x,y = NR(tau_Scomp1, tau_tau01, taufl_Scomp1, 1, 1)
# plt.loglog(x,y)
# plt.show()

# x, y = NR1(Smax, Sfl, S2, a0, D, A, m)

# G = 73.4e9 / (2 * (1 + 0.33))
# E = 73.4e9

# y_real = 4 / (np.pi * E) * Smax * 0.5 * D * y

# plt.loglog(x, y_real)
# plt.show()

# D = 8.71e-6
# E = 220e9
# G = 85e9
# miu = 0.294
# Smax = 340e6
# Kth = 0.79e6
# Sfl = Kth / (np.pi * D * 0.5) ** 0.5
# Scomp = 1.5e9
# S2 = Scomp
# a0 = 139e-6

# p1 = 1.6e-13
# p2 = 0.92
# dS = Smax * 0.9
# rou = p1 * (1.55 * dS) ** p2


# x, y, y_real = NR(Smax, Sfl, S2, D, miu, G)
# x = x * 1e6
# y_real *= 1e6*rou

# plt.loglog(x, y_real)
# plt.show()


# tau, tauFL, Kth, delt, miu, D, G = (
#     167.5 * 1e3,
#     80 * 1e3,
#     2.57e6,
#     0.5,
#     0.33,
#     6e-6,
#     73.4e9 / (2 * 1.33),
# )

# x, y, y_real = NRbyK(tau, tauFL, Kth, delt, miu, D, G)
# x /= Kth
# plt.loglog(x, y)
# plt.show()

G, miu, D, s1i, s2i, tauFL, r0, S, A2, m2 = (
    73.4e9 / (2 * 1.33),
    0.33,
    6e-6,
    0,
    460 * 1e3,
    80 * 1e3,
    0.5e-6,
    167.5 * 1e3,
    1,
    1,
)
a, CTOD, dadN = NR_rios(G, miu, D, s1i, s2i, tauFL, r0, S, A2, m2)
plt.loglog(a, CTOD)
plt.xlabel("a")
plt.ylabel("CTOD")
plt.grid(which="minor", axis="both", ls="-.")
plt.show()

plt.loglog(CTOD, dadN)
plt.xlabel("CTOD")
plt.ylabel("dadN")
plt.grid(which="minor", axis="both", ls="-.")
plt.show()

plt.loglog(a, dadN)
plt.xlabel("a")
plt.ylabel("dadN")
plt.grid(which="minor", axis="both", ls="-.")
plt.show()
