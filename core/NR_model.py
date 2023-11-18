import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# if __name__ == 'main':
#     ln = np.log


#     # eq3 : 作用在塑性区前方r0处位错源的应力
#     #   tau_tao0 = tau / tau0, tau >> tao0
#     tau_tao0s = np.array([5, 1.1])
#     x = np.arange(1, 1.1, 1e-5)
#     ns = np.array([1.0, 0.75, 0.5, 0.25])

#     for tau_tau0 in tau_tao0s:
#         for n in ns:
#             y = (
#                 1
#                 / (2**0.5)
#                 / ((x - 1) ** 0.5)
#                 * (1 - 2 / np.pi / tau_tau0 * np.arccos(n))
#                 + 1 / tau_tau0
#             )
#         #     plt.plot(x, y, label="n = " + str(n))
#         # plt.title("tau/tau0 = " + str(tau_tau0))
#         # plt.ylim(1, 25)
#         # plt.legend()
#         # plt.show()

#     # eq4 : tau_c为激活位错所需要的应力.当n=nc, 即裂尖接近晶界时，就会在r0处达到激活位错源的临界应力, 代表了滑移障碍的强度
#     # 重要的是要考虑它的变化范围
#     m_star = 3.1  # taylor type的均值, fcc金属为3.1
#     r0 = 1
#     tau_c = 1
#     S_zeta0 = 0.5 * m_star * tau_c
#     #   许多因素(位错锁定、变形过程解锁源的产生，位错密度的增加、应变硬化、晶粒内障碍的发展和滑移带宽度的变化)
#     # 都会影响tau_c和r0的值。在NR模型中，这两项为常数


#     # 考虑裂纹在后续晶粒中的扩展，注意到塑性区总是与晶界重合，即
#     i = np.arange(1, 100, 2)
#     D = 1
#     c = i * D / 2

#     # 由于r0假设为常数，则zeta0随着塑性区极限的位移的接近一个值。源zeta0在每个新晶界中的相对位置可以写成
#     zeta0 = 1
#     zeta0_i = 1 + 1 / i * (zeta0 - 1)

#     # 对于相同的n, 可以得到每个新源的应力
#     S_zeta0_i = i**0.5 * S_zeta0 * (zeta0_i - 1) ** 0.5 / ((zeta0 - 1) ** 0.5)


#     # 裂纹扩展的边界条件：疲劳极限、长裂纹扩展阈值
#     # 当n=1时, zeta0为
#     ai = i * D / 2
#     zeta0 = 1 + r0 / ai

#     # 此时得到使滑移带扩展所需要的应力，eq12
#     tau_Li = 0.5 * m_star * tau_c * (2 * r0) ** 0.5 / (i * D / 2) ** 0.5
#     # 特别地，当i=1，有
#     tau_FL = 0.5 * m_star * tau_c * (2 * r0) ** 0.5 / (D / 2) ** 0.5

#     y = []
#     x = []

#     for _ in range(len(i) - 1):
#         x.append(i[_])
#         x.append(i[_])
#         y.append(1 / i[_] ** 0.5)
#         y.append(1 / i[_ + 1] ** 0.5)
#     # plt.loglog(x, y)
#     # plt.show()

#     # 建立疲劳极限与应力强度因子阈值的关系
#     Kth = tau_FL * (np.pi * D * 0.5) ** 0.5
#     tau_FL=Kth/((np.pi * D * 0.5) ** 0.5)
#     #!!!裂纹扩展速率!!!
#     # 裂纹扩展速率认为与CTOD成正比，而CTOD等于进入塑性区的位错数×伯格斯矢量. 其无量纲形式为
#     tau_tau0 = 12
#     fai_bar = i * n * ln(1 / n) + (1 - n**2) ** 0.5 * (
#         np.arcsin(n) + np.pi * 0.5 * (tau_tau0 - 1)
#     )

#     # 当进入新晶粒时，n的临界值nc下降到一个新值ns
#     # n当前的临界值nc为 eq29
#     taufl_Scomp = 0.1
#     tau_Scomp = 0.3
#     nc = np.cos(np.pi * 0.5 * (tau_Scomp - taufl_Scomp / i**0.5))
#     ns = nc * i / (i + 2)

#     # plot
#     x = []
#     y = []
#     for _ in range(len(nc)):
#         if _ == 0:
#             n = np.arange(0.4, nc[_] + 1e-5, 1e-5)
#         else:
#             n = np.arange(ns[_ - 1], nc[_] + 1e-5, 1e-5)
#         fai_bar = i[_] * n * ln(1 / n) + (1 - n**2) ** 0.5 * (
#             np.arcsin(n) + np.pi * 0.5 * (tau_tau0 - 1)
#         )

#         y += list(fai_bar)
#         if _ == 0:
#             x += list(np.linspace(0.4, i[_], len(fai_bar)))
#         else:
#             x += list(np.linspace(i[_ - 1], i[_], len(fai_bar)))


#     plt.plot(x, y)
#     plt.xscale("log")
#     plt.yscale("log")
#     plt.show()

#     #!几个主要假设：
#     # nc的计算方法仅为假设，并不具有排他性
#     # 对于Scomp, Lardner(1968)和Tomkins(1968)采用了Su（极限强度）而不是Sy（屈服强度）
#     # Scomp的唯一要求：Scomp > tau

# 综上


def NR_old(tau_Scomp, tau_tau0, taufl_Scomp, A, m):
    # tau_FL = Kth / ((np.pi * D * 0.5) ** 0.5)
    # D = Kth / tau_FL**2 * 2 / np.pi
    # n0 = 162e-6 / D / 0.5
    ln = np.log
    i = np.arange(1, 1000, 2)
    nc = np.cos(np.pi * 0.5 * (tau_Scomp - taufl_Scomp / i**0.5))
    ns = np.append(0, nc * i / (i + 2))
    x = []
    y = []
    for _ in range(len(nc)):
        n = np.arange(ns[_], nc[_] + 1e-5, 1e-5)

        fai_bar = i[_] * (
            n * ln(1 / n)
            + (1 - n**2) ** 0.5 * (np.arcsin(n) + np.pi * 0.5 * (tau_tau0 - 1))
        )

        y += list(fai_bar)
        if _ == 0:
            x += list(np.linspace(0.1, i[_], len(fai_bar)))
        else:
            x += list(np.linspace(i[_ - 1], i[_], len(fai_bar)))

    y = np.array(y)
    x = np.array(x)
    y = A * y**m
    return x, y
    # plt.plot(x, y)
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.show()


def NR(tau, tauFL, S2i, D, miu, G):
    # 新论文
    ln = np.log
    i = np.arange(1, 100, 2)
    tauLi = tauFL / i**0.5
    nc = np.cos(np.pi * 0.5 * (tau - tauLi) / S2i)
    ns = np.append(0, nc * i / (i + 2))
    ai = 0.5 * D * i
    x = []
    y = []
    for _ in range(len(nc)):
        n = np.arange(ns[_], nc[_] + 1e-5, 1e-5)

        fai_bar = i[_] * (
            n * ln(1 / n)
            + (1 - n**2) ** 0.5 * (np.arcsin(n) - np.pi * 0.5 * (tau / S2i - 1))
        )

        y += list(fai_bar)
        if _ == 0:
            x += list(np.linspace(1e-6, ai[_], len(fai_bar)))
        else:
            x += list(np.linspace(ai[_ - 1], ai[_], len(fai_bar)))

    if len(x) == 0:
        print("stop!")
    y = np.array(y)
    x = np.array(x)
    k = 2 * (1 - miu) * D * S2i / np.pi / G
    y_real = k * y

    return x, y, y_real


def NRbyK(tau, tauFL, Kth, delt, miu, D, G):
    # 新论文，基于K
    ln = np.log
    i = np.arange(1, 1000, 2)

    K = np.linspace(Kth, 100 * Kth, len(i))

    def func(x):
        return np.cos(0.5 * np.pi * (1 - Kth * x**2 / K) ** (1 - delt)) - x

    nc = fsolve(func, [0.1] * len(K))

    ns = nc / (1 + 2 * nc * (tau / tauFL) ** 2 * (Kth / K) ** 2)
    ns = np.append([0], ns)

    x = np.array([])
    y = np.array([])

    for _ in range(len(i)):
        # n = np.arange(ns[_], nc[_] + 1e-5, 1e-5)
        n = ns[_]

        fai_bar = (
            (tauFL / tau)
            * (K[_] / Kth) ** 2
            * (ln(1 / n) + (1 - n**2) ** 0.5 / n * np.arcsin(n))
        )

        if np.isnan(fai_bar):
            continue

        if isinstance(fai_bar, float):
            fai_bar = np.array([fai_bar])
        else:
            fai_bar = np.array(fai_bar)

        y = np.append(y, fai_bar)

        # if _ == 0:
        #     x = np.append(x, np.linspace(0, K[_], len(fai_bar)))
        # else:
        #     x = np.append(x, np.linspace(K[_ - 1], K[_], len(fai_bar)))
        x = np.append(x, K[_])

    k = 2 * (1 - miu) * D * tauFL / np.pi / G
    y_real = k * y

    return x, y, y_real


def NR1(Smax, Sfl, S2, a0, D, A, m):
    # 计算nc的方法被修改，见
    # [1] DE PANNEMAECKER A, FOUVRY S, BUFFIERE J Y, 等.
    # Modelling the fretting fatigue crack growth: From short crack
    # correction strategies to microstructural approaches[J/OL].
    # International Journal of Fatigue, 2018, 117: 75-89.
    # DOI:10.1016/j.ijfatigue.2018.07.034.
    ln = np.log
    i = np.arange(1, 1000, 2)
    x = []
    y = []

    # 连续a的上界
    ai = 0.5 * D * i

    # 晶界障碍
    f = 2.5
    d = 0.5 * D
    Sth_mult = a0**0.5 / (ai**f + a0**f - d**f) ** (1 / (2 * f))
    nc = np.cos(np.pi * 0.5 * (Smax / S2 - Sfl / S2 * Sth_mult))
    ns = np.append([0], nc * i / (i + 2))

    for _ in range(len(nc)):
        # if _ == 0:
        #     n = np.linspace(0, nc[0], 1000)
        # else:
        #     n = np.linspace(nc[_ - 1], nc[_], 1000)
        n = np.linspace(ns[_], nc[_], 1000)

        # fai_bar = i[_] * (
        #     S2 / Smax * n * ln(1 / n)
        #     + np.pi
        #     * 0.5
        #     * (1 - n**2) ** 0.5
        #     * (1 - 2 / np.pi * S2 / Smax * np.arccos(n))
        # )

        fai_bar = i[_] * (
            n * ln(1 / n)
            + (1 - n**2) ** 0.5 * (np.arcsin(n) + np.pi * 0.5 * (-Smax / S2 + 1))
        )

        y += list(fai_bar)
        if _ == 0:
            x += list(np.linspace(1e-6, ai[_], len(fai_bar)))
        else:
            x += list(np.linspace(ai[_ - 1], ai[_], len(fai_bar)))

    # plt.scatter(i, ai)
    # plt.show()

    y = np.array(y)
    x = np.array(x)
    y = A * y**m
    return x, y


def NR2(Smax, Sfl, S2, a0, D):
    # 速率计算不同，见：
    # MAO J, HU D, MENG F, 等.
    # Multiscale modeling of transgranular short crack growth
    # during fatigue in polycrystalline metals[J/OL].
    # International Journal of Fatigue, 2018, 116: 648-658.
    # DOI:10.1016/j.ijfatigue.2018.07.017.
    ln = np.log
    i = np.arange(1, 1000, 2)
    x = []
    y = []

    # 连续a的上界
    ai = 0.5 * D * i

    # 晶界障碍
    f = 2.5
    d = 0.5 * D
    Sth_mult = a0**0.5 / (ai**f + a0**f - d**f) ** (1 / (2 * f))
    nc = np.cos(np.pi * 0.5 * (Smax / S2 - Sfl / S2 * Sth_mult))
    ns = np.append([0], nc * i / (i + 2))

    for _ in range(len(nc)):
        # if _ == 0:
        #     n = np.linspace(0, nc[0], 1000)
        # else:
        #     n = np.linspace(nc[_ - 1], nc[_], 1000)
        n = np.linspace(ns[_], nc[_], 1000)

        # fai_bar = i[_] * (
        #     S2 / Smax * n * ln(1 / n)
        #     + np.pi
        #     * 0.5
        #     * (1 - n**2) ** 0.5
        #     * (1 - 2 / np.pi * S2 / Smax * np.arccos(n))
        # )

        fai_bar = i[_] * (
            n * ln(1 / n)
            + (1 - n**2) ** 0.5 * (np.arcsin(n) + np.pi * 0.5 * (-Smax / S2 + 1))
        )

        y += list(fai_bar)
        if _ == 0:
            x += list(np.linspace(1e-6, ai[_], len(fai_bar)))
        else:
            x += list(np.linspace(ai[_ - 1], ai[_], len(fai_bar)))

    # plt.scatter(i, ai)
    # plt.show()

    y = np.array(y)
    x = np.array(x)
    return x, y


def NR_rios(G, miu, D, s1i, s2i, tauFL, r0, S, A2, m2):
    # 见文献：
    # [1] DE LOS RIOS E R, TRULL M, LEVERS A.
    # Modelling fatigue crack growth in shot‐peened components of Al 2024‐T351[J/OL].
    # Fatigue & Fracture of Engineering Materials & Structures, 2000, 23(8): 709-716.
    # DOI:10.1046/j.1460-2695.2000.00287.x.

    i = np.arange(1, 2000, 2)
    # tauLi = tauFL / i**0.5
    m = 1 + 0.5 * np.log(i)
    tauLi = s1i + m / m[0] * (tauFL - s1i) / i**0.5
    c = i * D / 2 + r0
    nc = np.cos(np.pi * 0.5 * (S - tauLi) / s2i)
    ns = np.append([0.2], nc * i / (i + 2))

    x = np.array([])
    y = np.array([])

    for _ in range(len(i)):
        nsi = ns[_]
        nci = nc[_]
        ci = c[_]
        n1i = np.linspace(nsi, nci, 50)
        n2i = i[_] * D / 2 / (i[_] * D / 2 + r0)
        alfa1 = s2i / np.arccos(n2i)
        alfa2 = -s2i
        alfa3 = np.pi * 0.5 * S / np.arccos(n2i)
        s3 = alfa1 * np.arccos(n1i) + alfa2 + alfa3

        ai = n1i * ci

        CTOD = (D * 0.5 + r0) * (
            s2i * 2 * n1i * np.arccosh(abs((1 + n1i**2) / (2 * n1i)))
            - (s3 - s2i)
            * (
                (n1i - n2i) * np.arccosh(abs((1 - n1i * n2i) / (n2i - n1i)))
                - (n1i + n2i) * np.arccosh(abs((1 + n1i * n2i) / (n2i + n1i)))
            )
        )
        CTOD *= 2 * (1 - miu) / G

        y = np.append(y, CTOD)
        x = np.append(x, ai)

    dadN = A2 * y**m2

    return x, y, dadN


# NR1(0, 0, 0, 0, 6e-6, 0, 0)

# a = np.arange(1e-6, 1e-2, 1e-8)
# Sfl = 80 * 1e6
# a0 = 162 * 1e-6
# D = 6e-6
# d = 0.5 * D
# fs = [1, 1.5, 2, 2.5, 3, 3.5, 4]
# # f = 2.5
# for f in fs:
#     Sth = a0**0.5 / (a**f + a0**f - d**f) ** (1 / (2 * f))
#     plt.loglog(a, Sth, label=str(f))
#     plt.plot([a0] * 1000, np.linspace(Sth.min(), Sth.max(), 1000))
# plt.legend()
# plt.show()

# S2 = 460 * 1e3
# Smax = 167.5 * 1e3

# n = np.arange(0, 1, 1e-5)
# fai_bar = 1 * (
#     n * np.log(1 / n)
#     + (1 - n**2) ** 0.5 * (np.arcsin(n) + np.pi * 0.5 * (-Smax / S2 + 1))
# )

# plt.plot(n, fai_bar)
# plt.show()

# tauFL = 1
# tau = 0.3
# Kth = 1
# K = np.linspace(1, 20, 1000)
# ln = np.log
# delt = 0.5


# def func(x):
#     return np.cos(0.5 * np.pi * (1 - Kth * x**2 / K) ** (1 - delt)) - x


# nc = fsolve(func, [1] * len(K))

# n = nc[0]

# fai_bar = (
#     (tauFL / tau)
#     * (K / Kth) ** 2
#     * (ln(1 / n) + (1 - n**2) ** 0.5 / n * np.arcsin(n))
# )

# plt.plot(K, fai_bar)
# plt.show()
