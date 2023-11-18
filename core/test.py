from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from CT import get_dadN
import material as mat
from NR_model import NR_rios


a, CTOD, dadN = NR_rios(
    G=mat.G,
    miu=mat.miu,
    D=mat.D,
    s1i=0,
    s2i=mat.s2i,
    tauFL=mat.tauFL,
    r0=mat.r0,
    S=mat.S,
    A2=1,
    m2=1,
)

a *= 1e3
dadN *= 1e3


plt.plot(a, dadN, marker="x", linewidth=1, markersize=3)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("a (mm)")
plt.ylabel("CTOD (mm/cycle)")
plt.grid(which="minor", axis="both", ls="-.")
plt.legend()
plt.savefig("./NR.svg", format="svg", dpi=600)
plt.show()


# n1 = np.arange(0, 1, 1e-3)
# D = mat.D
# r0 = mat.r0
# s2 = 1245.7428
# S = mat.S
# n2 = D * 0.5 / (D * 0.5 + r0)
# alfa1 = s2 / np.arccos(n2)
# alfa2 = -s2
# alfa3 = np.pi * 0.5 * S / np.arccos(n2)
# s3 = alfa1 * np.arccos(n1) + alfa2 + alfa3
# CTOD = (D * 0.5 + r0) * (
#     s2 * 2 * n1 * np.arccosh(abs((1 + n1**2) / (2 * n1)))
#     - (s3 - s2)
#     * (
#         (n1 - n2) * np.arccosh(abs((1 - n1 * n2) / (n2 - n1)))
#         - (n1 + n2) * np.arccosh(abs((1 + n1 * n2) / (n2 + n1)))
#     )
# )
# plt.plot(n1, CTOD,label='i = 1')
# plt.xlabel("n1")
# plt.ylabel("CTOD_bar")
# plt.legend()
# plt.show()

# n1 = np.arange(0, 1, 1e-3)
# D = mat.D
# r0 = mat.r0
# s2 = 1245.7428
# S = mat.S
# n2 = D * 0.5 / (D * 0.5 + r0)
# alfa1 = s2 / np.arccos(n2)
# alfa2 = -s2
# alfa3 = np.pi * 0.5 * S / np.arccos(n2)
# plt.plot(
#     n1,
#     alfa1 * np.arccos(n1) + alfa2 + alfa3,
#     label="alfa1 = "
#     + str(alfa1)
#     + "\nalfa2 = "
#     + str(alfa2)
#     + "\nalfa3 = "
#     + str(alfa3)
#     + "\ni = 1",
# )
# plt.xlabel("n1")
# plt.ylabel("s3 (MPa)")
# plt.legend()
# plt.show()

# i = np.arange(1, 100, 2)

# mi = 1 + 0.5 * np.log(i)

# tauLi_0p1 = mi / mi[0] * 601 / i**0.5
# tauLi_minus1 = mi / mi[0] * 409 / i**0.5

# x = []
# y1 = []
# y2 = []
# for _ in range(len(i) - 1):
#     x.append(i[_])
#     x.append(i[_])
#     y1.append(tauLi_0p1[_])
#     y1.append(tauLi_0p1[_ + 1])
#     y2.append(tauLi_minus1[_])
#     y2.append(tauLi_minus1[_ + 1])

# plt.loglog(x, y1, label="R = 0.1", marker="s", linewidth=1, markersize=1)
# plt.loglog(x, y2, label="R = -1", marker="o", linewidth=1, markersize=1)
# plt.xlabel("i")
# plt.ylabel("tauLi (MPa)")
# plt.grid(which="both", axis="both", ls="-.")
# plt.legend()
# plt.show()


# s2i, r0, S, A2, m2 = (
#     1528518311.9880035,
#     2.788711023863334e-09,
#     431051452.48999935,
#     0.7833846323472584,
#     1.0163315674781184,
# )
# NR_a, NR_CTOD, NR_dadN = NR_rios(
#     mat.G, mat.miu, mat.D, mat.s1i, s2i, mat.tauFL, r0, S, A2, m2
# )


# NR_a *= 1e3
# NR_dadN *= 1e3

# plt.plot(NR_a, NR_dadN, label="Short crack", marker="x", linewidth=1, markersize="2")
# plt.xscale("log")
# plt.yscale("log")
# plt.xlabel("a (mm)")
# plt.ylabel("dadN (mm/cycle)")
# plt.grid(which="minor", axis="both", ls="-.")
# plt.legend()
# plt.show()


# D = 6e-6
# r0 = 0.5e-6
# s2i = 500 * 1e3
# s1i = 0
# G = 73.4e9 / (2 * 1.33)
# miu = 0.33
# S = 167.5 * 1e3

# ci = D / 2
# n1i = np.linspace(0, 1, 1000)
# n2i = 1 * D / 2 / (1 * D / 2 + r0)

# s3i = (
#     (s2i - s1i) * np.arcsin(n1i) - s2i * np.arcsin(n2i) + np.pi * 0.5 * S
# ) / np.arccos(n2i)

# ai = n1i * ci

# CTOD = 2 * ai * (1 - miu) / (G * np.pi * n1i)
# tmp = 0
# tmp -= (s2i - s1i) * (-(n1i + n1i) * np.arccosh(abs((1 + n1i**2) / (n1i + n1i))))
# # tt = abs((1 - n2i * n1i) / (n2i - n1i))
# # print(tt.max(), tt.min())
# # tt2 = abs((1 + n2i * n1i) / (n2i + n1i))
# # print(tt2.max(), tt2.min())
# tmp -= (s3i - s2i) * (
#     (n1i - n2i) * np.arccosh(abs((1 - n2i * n1i) / (n2i - n1i)))
#     - (n1i + n2i) * np.arccosh(abs((1 + n2i * n1i) / (n2i + n1i)))
# )
# CTOD *= tmp

# plt.plot(n1i, CTOD, marker="o")
# plt.show()
