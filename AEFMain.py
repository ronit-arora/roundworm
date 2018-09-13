import numpy as np
import matplotlib.pyplot as plt
from ModelAEIF import *

tstep = 0.05E-3

# Check Syntax
# Types = {'RS','RS','RS','RS','RS','RS','RS','RS','RS','RS'}
Types = np.array(['RS', 'RS', 'RS', 'RS', 'RS', 'RS', 'RS', 'RS', 'RS', 'RS'])

tau = 15E-3 * np.ones((10))
# print(tau)
tau_s = tau / 4
tau_v = 15E-3
Tmax = 40

# From TempData1.mat file - directly placed into code here
Mean = np.array([[-0.0269, 0.0247],
                 [ 0.0013, 0.0020],
                 [ 0.0113, 0.0199],
                 [ 0.0162, -0.0139],
                 [ 0.0259, 0.0144],
                 [-0.0260, -0.0371],
                 [ 0.0316, -0.0304]])

Weight = np.array([[-0.4795],
                   [ 0.4326],
                   [ 0.4237],
                   [-0.3365],
                   [-0.3193],
                   [-0.2446],
                   [ 0.1537]])

Sigma = 0.0100

X = np.zeros((800001, 2))
Y = np.zeros((800001, 10, 2))

# Weird Function Handle Stuff PLEASE FIX
# Temp = @(x,y) (5*Temperature(x,y,Mean,Weight,Sigma)-80);
# I = @(x,y) Iinput(Temp, x, y);
def Temp(x, y):
    return 5 * Temperature(x, y, Mean, Weight, Sigma) - 80

def I(x, y):
    return Iinput(Temp(x,y), x, y)

N = 10
Syn_const_mat_2 = np.zeros((2, N, N))
Syn_const_mat_2[0, :, :] = 1E-12

# I0 for each synapse, W65 is -300, W98 is -400 (other is -180)
Syn_const_mat_2[1, :, :] = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                     [-205, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                     [207, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                     [0, 120, 0, 0, 0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, -50, 0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 200, -180, -200, 0, 0, 0, 0],
                                     [0, 0, 100, 0, 0, 0, 0, 0, 0, 0],
                                     [0, 0, 0, 0, 0, 0, -50, 0, 0, 0],
                                     [0, 0, 0, 0, 0, 0, 200, -180, -200, 0],
                                     [0, 6.08915, 6.75506575, 0, 0, -1000, 0, 0, -1000, -800]])
# X0 = X[-2, :]
X0 = [-0.02, 0.02]
# X0 = np.array([[0], [0]])

# Check this statement
theta0 = np.angle(X[-2, 0] - X[-3, 0] + 1j * (X[-2, 1] - X[-3, 1]))

# print(theta0)
Vel_ss = 1E-3
pump_vel = 1.3E-3

# print(X[0, 0])
# print(X[0, 1])

# print(I(int(X[0, 0]), int(X[0, 1])))

#permute???, I????? PLEASE FIX			I(int(X[0, 0]), int(X[0, 1]))
(T, Y, Iout, X, SpikeCount, Vel, W65, W98) = AEF_Simulation([40, 80], np.zeros((1, 2, 10)), tstep, Types, I, Syn_const_mat_2, tau, tau_s, X0, theta0, tau_v, pump_vel, Vel_ss)

print("SPIKECOUNT: " + str(SpikeCount))


TempHist = np.zeros((1, 800001))

# for a in range(0, 800000):
#     TempHist[0, a] = Temp(X[a, 0], X[a, 1])
TempHist[0, :] = Temp(X[:, 0], X[:, 1])


# PLOT SECTION

Tstart_plot = 40
Tstop_plot = 80

# Tsel = (T>Tstart_plot) and (T<Tstop_plot)
Tsel = np.ones((800001, 1))
Tsel[[0, 800000], 0] = 0

try:
    plt.scatter(X[:,0], X[:, 1])
    plt.plot(X[:,0], X[:, 1])
    plt.show()
except:
    print("BAD")
