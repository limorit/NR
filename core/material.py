import numpy as np

# CT 单位mm
P = 4e3
B = 10
W = 40

R = 0.1
n = 3.33
C = (9.37e-9) / 10 ** (3 / 2 * n)

an = 0.2 * W
a0 = 1 / np.pi * (10 * 0.9 * 1e6 / (601 * 1e6)) ** 2 * 1e3
# print(a0)


# NR-TC4    标准单位

s2i = 1245.7428  # MPa
s2i *= 1e6

G, miu, D, s1i, s2i, tauFL, r0, S, A2, m2 = (
    44e9,  # G
    0.29,  # miu
    25e-6,  # D
    0,  # s1i
    s2i,  # s2i
    601e3,  # tauFL
    0.5e-6,  # r0, unknow
    1000e6,  # S
    5000,  # A2, unknow
    1.8,  # m2, unknow
)
