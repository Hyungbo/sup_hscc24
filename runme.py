import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.signal import place_poles

def ip(t, x, mode):
    M = 0.5
    m = 0.2
    b = 0.1
    g = 9.8
    L = 0.1 + np.mod(t, 0.5)
    I = m * L**2
    p = (M + m) * I + M * m * L**2

    A = np.array([[0, 1, 0, 0],
                  [0, -(I + m * L**2) * b / p, (m**2 * g * L**2) / p, 0],
                  [0, 0, 0, 1],
                  [0, -(m * L * b) / p, m * g * L * (M + m) / p, 0]])

    B = np.array([[0],
                  [(I + m * L**2) / p],
                  [0],
                  [m * L / p]])

    poles = [-8.5 + 7.9j, -8.5 - 7.9j, -4.8 + 0.8j, -4.8 - 0.8j]
    K = place_poles(A, B, poles).gain_matrix

    Kprime = K + 20 * np.sin(2 * np.pi * 20 * t) * np.ones((1, 4))

    if mode == 0:
        dxdt = A @ x - B @ (K @ x)
    elif mode == 1:
        dxdt = A @ x - B @ (Kprime @ x)
    elif mode == 2:
        dxdt = np.array([L, I, *K.flatten(), *Kprime.flatten()])

    return dxdt

# Simulation time and initial conditions
tw = [0, 2]
ic = [-0.2, 0, 0, 0]

# Solving the ODEs
t0 = np.linspace(tw[0], tw[1], 100)
x0 = odeint(lambda x, t: ip(t, x, 0), ic, t0)

t1 = np.linspace(tw[0], tw[1], 100)
x1 = odeint(lambda x, t: ip(t, x, 1), ic, t1)

# Get values for L, I, K, and K'
rec = np.array([ip(t, 0, 2) for t in t1])

# Plotting
plt.figure(figsize=(10, 8))

plt.subplot(4, 1, 1)
plt.plot(t1, rec[:, 0], '-', label='Moment of Inertia')
plt.plot(t1, rec[:, 1], '--', label='Ball position')
plt.legend(fontsize=12)
plt.title('Time varying parameters', fontsize=15, fontweight='bold')

plt.subplot(4, 1, 2)
plt.plot(t0, x0[:, 2], '-', label='Angle of Tube')
plt.plot(t0, x0[:, 0], '--', label='Position of Cart')
plt.axis([0, 2, -0.25, 0.1])
plt.grid()
plt.legend(loc='lower right', fontsize=12)
plt.title('Response of Average System', fontsize=15, fontweight='bold')

plt.subplot(4, 1, 3)
plt.plot(t1, -rec[:, 6], '-', label='Used Gain K_1(t)')
plt.plot(t1, -rec[:, 2], '--', label='Designed Gain K_1(t)')
plt.legend(fontsize=12)
plt.title('Comparison between Designed Gain and Used Gain', fontsize=15, fontweight='bold')

plt.subplot(4, 1, 4)
plt.plot(t1, x1[:, 2], '-', label='Angle of Tube')
plt.plot(t1, x1[:, 0], '--', label='Position of Cart')
plt.axis([0, 2, -0.25, 0.1])
plt.grid()
plt.legend(loc='lower right', fontsize=12)
plt.title('Response of Closed-loop System', fontsize=15, fontweight='bold')
plt.xlabel('time', fontsize=15)

plt.tight_layout()
plt.show()
