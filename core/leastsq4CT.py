import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from CT import get_dadN
from NR_model import NR, NR_rios

import material as mat

time_start = time.time()

# CT长裂纹试验数据 标准单位
a_start = mat.a0
a_end = 3
a = np.linspace(a_start, a_end, 100)
dadN, dK = get_dadN(a, mat.P, mat.B, mat.W, mat.C, mat.n, mat.R)
a *= 1e-3
dadN *= 1e-3

try_times = 0


# 需要拟合的函数
def func(p, x):
    _r0, _S, _A2, _m2 = p
    _a, _CTOD, _dadN = NR_rios(
        mat.G, mat.miu, mat.D, mat.s1i, mat.s2i, mat.tauFL, _r0, _S, _A2, _m2
    )

    _ans = []
    cant_find = False
    for xi in x:
        for i in range(len(_a) - 1):
            if xi <= _a[0]:
                _ans.append(_dadN[0])
                cant_find = False
                break
            elif xi >= _a[-1]:
                _ans.append(_dadN[-1])
                cant_find = False
                break
            elif _a[i] <= xi <= _a[i + 1]:
                if abs((_dadN[i] - _dadN[i + 1]) / max([_dadN[i], _dadN[i + 1]])) > 0.8:
                    _ans.append(max([_dadN[i], _dadN[i + 1]]))
                else:
                    _y0 = _dadN[i]
                    _y1 = _dadN[i + 1]
                    _x0 = _a[i]
                    _x1 = _a[i + 1]
                    tmp = _y0 + (_y1 - _y0) * (xi - _x0) / (_x1 - _x0)
                    _ans.append(tmp)
                cant_find = False
                break
            else:
                cant_find = True
        if cant_find:
            raise ValueError(
                "Cant find! _a[0] = "
                + str(_a[0])
                + ", _a[-1] = "
                + str(_a[-1])
                + ", xi = "
                + str(xi)
            )
    if len(_ans) == 0:
        raise ValueError(
            "Wrong in finding _ans!"
            + "\n"
            + "_a[0] = "
            + str(_a[0])
            + ", _a[-1] = "
            + str(_a[-1])
            + ", xi = "
            + str(xi)
        )
    return np.array(_ans)


# 定义偏差函数
def error(p, x, y):
    global try_times
    try_times += 1
    print("try times = " + str(try_times))
    _r0, _S, _A2, _m2 = p
    if (
        _r0 <= 0
        or _S <= 0
        or _A2 <= 0
        or _m2 <= 0
        # or (_S - mat.tauFL) < 0
        # or (_S - mat.tauFL) / mat.s2i < -1
        # or (_S - mat.tauFL) / mat.s2i > 1
        or np.cos(np.pi * 0.5 * (_S - mat.tauFL) / mat.s2i) < 0
    ):
        return [1e6] * len(y)
    else:
        return func(p, x) - y


# 设初始值
p0 = [mat.r0, mat.S, mat.A2, mat.m2]

# 最小二乘拟合
params = leastsq(error, p0, args=(a, dadN))

r0, S, A2, m2 = params[0]
print("r0, S, A2, m2 = " + str(r0) + ", " + str(S) + ", " + str(A2) + ", " + str(m2))
print("fit code = " + str(params[1]))

a *= 1e3
dadN *= 1e3
plt.scatter(
    a, dadN, c="none", label="Long crack", marker="o", edgecolors="r", linewidths=1, s=2
)

NR_a, NR_CTOD, NR_dadN = NR_rios(
    mat.G, mat.miu, mat.D, mat.s1i, mat.s2i, mat.tauFL, r0, S, A2, m2
)


time_end = time.time()

time_sum = time_end - time_start
print("running time = " + str(time_sum))

NR_a *= 1e3
NR_dadN *= 1e3

plt.plot(NR_a, NR_dadN, label="Short crack", marker="x", linewidth=1, markersize="2")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("a (mm)")
plt.ylabel("dadN (mm/cycle)")
plt.grid(which="minor", axis="both", ls="-.")
plt.legend()
plt.savefig("./NR.svg", format="svg", dpi=600)
plt.show()
