import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

G = 86e3
miu = 0.303
d = 27.3e-3
Sy = 1070
Smax = 800
Smin = 800 * 0.1
A = G / (2 * np.pi * (1 - miu))


m = 1.53778
C = 0.64525


s = Smax - Smin  # 外载


# A:平衡滑移带阶段, a < c < d
s1A = 2 * s / 0.5  # 滑移带中晶粒A向上移动时的摩擦阻力
s2B = 2 * 2 * s1A  # 晶粒A和晶粒B之间的摩擦力
a_limA = d * np.cos(np.pi * s / (2 * s1A))
aA = np.arange(10e-3, a_limA, 1e-5)
cA = aA / np.cos(np.pi * s / (2 * s1A))
sec = 1 / (np.cos(np.pi * s / (2 * s1A)))
fai_tA = (2 * s1A * aA / (np.pi**2 * A)) * np.log(sec)

# B:塞积滑移带阶段, a < c = d
s1B =2 * s / 0.9  # 滑移带中晶粒A向上移动时的摩擦阻力
s2B =2 * 2 * s1B  # 晶粒A和晶粒B之间的摩擦力
a_limB = d * np.cos(np.pi * s1B / (2 * s1A))
cB = d
aB = np.arange(a_limA, a_limB, 1e-5)
beta = 1 - (2 * s1B / (np.pi * s)) * np.arccos(aB / cB)
fai_tB = (beta * s / (np.pi * A)) * (cB**2 - aB**2) ** 0.5 + (
    2 * s1B * aB / (np.pi**2 * A)
) * np.log(cB / aB)

# C:扩展滑移带阶段, a < d < c
s1C =2 * s / 0.9  # 滑移带中晶粒A向上移动时的摩擦阻力
s2C =2 * 2 * s1C  # 晶粒A和晶粒B之间的摩擦力
aC = np.arange(a_limB, d, 1e-5)


def func14(_x):
    return (
        np.arccos(aC / _x) + (s2C / s1C - 1) * np.arccos(d / _x) - np.pi * s / (2 * s1C)
    )


cC = np.array(fsolve(func14, aC))
tmpa = (cC**2 - aC**2) ** 0.5
tmpd = (cC**2 - d**2) ** 0.5
tmp1 = d * np.log(abs((tmpd + tmpa) / (tmpd - tmpa)))
tmp2 = aC * np.log(abs((aC * tmpd + d * tmpa) / (aC * tmpd - d * tmpa)))
gacd = tmp1 - tmp2
fai_tC = (2 * s1C * aC / (np.pi**2 * A)) * np.log(cC / aC) + (
    (s2C - s1C) / (np.pi**2 * A)
) * gacd


a = list(aA * 1e3) + list(aB * 1e3) + list(aC * 1e3)
fai = list(fai_tA) + list(fai_tB) + list(fai_tC)
plt.scatter(a, fai, marker="x", linewidths=0.2)
plt.xscale("log")
plt.yscale("log")
plt.show()
