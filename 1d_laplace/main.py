import matplotlib.pyplot as plt
import numpy as np

def u(x):
    return x ** 6

def func(x):
    return -30 * x ** 4

def TDA(a, b, c, f):
    n = len(f)
    alpha = np.zeros(n)
    beta = np.zeros(n)
    x = np.zeros(n)

    alpha[0] = -c / b
    beta[0] = f[0] / b
    for i in range(1, n - 1):
        alpha[i] = -c / (b + a * alpha[i - 1])
        beta[i] = (f[i] - a * beta[i - 1]) / (b + a * alpha[i - 1])
    x[-1] = (f[-1] - a * beta[-2]) / (b + a * alpha[-2])

    for i in range(n - 2, -1, -1):
        x[i] = alpha[i] * x[i + 1] + beta[i]
    return x

l2_norms = []
c_norms = []
grid = np.linspace(10, 10000, 10, dtype='int')
for n in grid:
    x = np.linspace(0, 1, n + 1)
    f = np.array([func(x[i]) for i in range(1, n)]) / n ** 2
    f[0] += u(0)
    f[-1] += u(1)
    
    ans = np.array([u(x[i]) for i in range(1, n)])
    sol = TDA(-1, 2, -1, f)
    l2_norms.append(np.linalg.norm(abs(sol - ans)))
    c_norms.append(max(abs(sol - ans))) 
    
plt.plot(grid, l2_norms, label='l2 norm')
plt.plot(grid, c_norms, label='C norm')
plt.xlabel('len')
plt.ylabel('norms')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('res.png')
plt.show()
