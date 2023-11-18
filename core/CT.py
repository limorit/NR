import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from NR_model import NR_rios

# 单位均为mm


def get_K(a, P, B, W):
    F = (
        29.6 * (a / W) ** 0.5
        - 185.5 * (a / W) ** (3 / 2)
        + 655.7 * (a / W) ** (5 / 2)
        - 1017 * (a / W) ** (7 / 2)
        + 638.9 * (a / W) ** (9 / 2)
    )
    K0 = P / (B * W**0.5)

    return F * K0


def get_dadN(a, P, B, W, C, n, R):
    dK = get_K(a, P, B, W) - get_K(a, P * R, B, W)
    return C * dK**n, dK


# P = 100
# B = 10
# W = 40
# Sy = 920

# an = 0.2 * W
# a = np.arange(an, W / 2, 1e-3)

# ans = get_K(a, P, B, W)


# df = pd.DataFrame({"a": a, "ans": ans})
# df.to_excel("./CT.xlsx", header=None)

# plt.plot(a, ans)
# plt.show()

# P = np.arange(0, 100000, 1)


# P = 4e3
# B = 10
# W = 40

# R = 0.1
# n = 3.33
# C = (9.37e-9) / 10 ** (3 / 2 * n)

# an = 0.2 * W
# a0 = 1 / np.pi * (10 * 0.9 * 1e6 / (601 * 1e6)) ** 2 * 1e3
# a_start = a0
# a_end = 1
# a = np.arange(a_start, a_end, 1e-3)


# # print(a0)


# dadN, dK = get_dadN(a, P, B, W, C, n, R)
# dK *= 10 ** (-3 / 2)

# min1 = dadN.min()
# max1 = dadN.max()

# plt.loglog(P, dK)
# plt.xlabel("P")
# plt.ylabel("dK")
# plt.grid(which="minor", axis="both", ls="-.")
# plt.show()

# a = a - an

# plt.loglog(dK, dadN)
# plt.xlabel("dK (MPa*m^0.5)")
# plt.ylabel("dadN mm/cycle")
# plt.grid(which="minor", axis="both", ls="-.")
# plt.show()

# plt.loglog(a, dK)
# plt.xlabel("a (mm)")
# plt.ylabel("dK (MPa*m^0.5)")
# plt.grid(which="minor", axis="both", ls="-.")
# plt.show()

# plt.loglog(a, dadN, label="Long crack")
# plt.xlabel("a (mm)")
# plt.ylabel("dadN (mm/cycle)")
# plt.grid(which="minor", axis="both", ls="-.")
# # plt.show()

# # 此处为标准单位
# HK = 371
# HV = 3033 * HK / (3066 / (6.255555 / (1088 - HK) + 0.9077) + abs(HK - 355))
# s2i = 3.85 * HV - 136.9
# s2i *= 1e6

# G, miu, D, s1i, s2i, tauFL, r0, S, A2, m2 = (
#     44e9,
#     0.29,
#     25e-6,
#     0,
#     s2i,
#     601e3,
#     0.5e-6,
#     500e6,
#     0.1,
#     1,
# )


# a, CTOD, dadN = NR_rios(G, miu, D, s1i, s2i, tauFL, r0, S, A2, m2)
# a *= 1e3
# dadN *= 1e3

# max2 = dadN.max()
# min2 = dadN.min()

# ytmp = np.linspace(min([min1, min2]), max([max1, max2]), 2)
# plt.plot([a0] * len(ytmp), ytmp, color="r", linestyle="-.")

# plt.loglog(a, dadN, label="Short crack")
# plt.xlabel("a (mm)")
# plt.ylabel("dadN (mm/cycle)")
# plt.grid(which="minor", axis="both", ls="-.")
# plt.legend()
# plt.show()
