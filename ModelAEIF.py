import numpy as np
import math
import matplotlib.pyplot as plt


def AEF_dY_dt(Iin, Yin, Params):

    # print(Yin)
    # Get 1st column (0th column)
    V = Yin[:, 0]

    # Get 2nd column (1st column) (Every variable down below has similar comment, reduced number by 1)
    U = Yin[:, 1]

    C = Params[:, 0]
    gL = Params[:, 1]
    El = Params[:, 2]
    VT = Params[:, 3]
    del_T = Params[:, 4]
    a = Params[:, 5]
    tau_w = Params[:, 6]
    b = Params[:, 7]
    Vr = Params[:, 8]

    N = np.size(V)
    dY_dt = np.zeros((N, 2))

    # Use numpy syntax
    # dY_dt[:,1] = (-gL .*(V - El) + gL.*del_T.*exp((V-VT)./del_T) - U + Iin)./C
    dY_dt[:, 0] = (-gL * (V - El) + gL * del_T * np.exp((V - VT) / del_T) - U + Iin) / C
    # print("-gL")
    # print(-gL)
    # print("V")
    # print(V)
    # print("El")
    # print(El)
    # print("Del_t")
    # print(del_T)
    # print("VT")
    # print(VT)
    # print("U")
    # print(U)
    # print("Iin")
    # print(Iin)
    # print("c")
    # print(C)

    # dY_dt(:,2) = (a.*(V - El) - U)./tau_w;
    dY_dt[:, 1] = (a * (V - El) - U) / tau_w

    return dY_dt

def Temp_Noise(x,y):
    a = 1000

    # T = sin(5*a*x)*sin(5*a*sqrt(abs(y))) + ...
    #     sin(6*a*sqrt(abs(x)))*sin(9*a*y) + ...
    #     sin(4*a*sqrt(abs(x)))*sin(6*a*sqrt(abs(y))) + ...
    #     sin(7*a*x)*sin(5*a*y);

    T = math.sin(5 * a * x) * math.sin(5 * a * math.sqrt(abs(y))) + \
        math.sin(6 * a * math.sqrt(abs(x))) * math.sin(9 * a * y) +\
        math.sin(4 * a * math.sqrt(abs(x))) * math.sin(6 * a * math.sqrt(abs(y))) + \
        math.sin(7 * a * x) * math.sin(5 * a * y)

    return T / 70

def Temperature(x, y, Mean, Weight, Sigma):

    T = 19.9
    for i in range(0, Mean.shape[0]):

        X = x - Mean[i, 0]
        Y = y - Mean[i, 1]
        # T = T + Weight[i] * exp(-(X^2+Y^2)/Sigma^2)
        T = T + Weight[i] * np.exp(-(X ** 2 + Y ** 2) / (Sigma ** 2))
    return T

# MAKE SURE TO HAVE X AND Y BE DECREMENTED BY 1
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def Iinput(Temp, x, y):

    x = x - 1
    y = y - 1
    # I = zeros(10,1);
    # Tc = 20;
    # Tnow = Temp(x,y);

    # I(1) = 600 + (Tnow-Tc)*500;
    # I(2) = 830.5;
    # I(3) = -396;
    # I(4) = 600;
    # I(5) = 800;
    # I(6) = 0;
    # I(7) = 600;
    # I(8) = 800;
    # I(9) = 0;
    # I(10) = 205;

    # I = I * 1E-12;

    I = np.zeros((10, 1))
    Tc = 20
    Tnow = Temp

    I[0] = 600 + (Tnow - Tc) * 500
    I[1] = 830.5
    I[2] = -396
    I[3] = 600
    I[4] = 800
    I[5] = 0
    I[6] = 600
    I[7] = 800
    I[8] = 0
    I[9] = 205

    I = I * 1E-12

    return I

def AEF_Simulation(Tin, Y0, stepT, type, Iin, Syn_const_mat, tau, tau_s, X0, theta0, tau_v, pump_vel, Vel_ss):
    Params = np.array([ [200, 10, -70, -50, 2, 2, 30, 0, -58],
                        [130, 18, -58, -50, 2, 4, 150, 120, -50],
                        [200, 10, -58, -50, 2, 2, 120, 100, -46]    ])

    OrderList = np.array([1E-12, 1E-9, 1E-3, 1E-3, 1E-3, 1E-9, 1E-3, 1E-12, 1E-3])
    Params = (np.ones((3, 1)) * OrderList) * Params

    # Check this statement - it's supposed to make a logical array
    Syn_is_conn = (Syn_const_mat[1, :, :] != 0)
    # print(Syn_is_conn.shape)

    # is_RS = strcmp(type, 'RS')
    # is_IB = strcmp(type, 'IB')
    # is_CH = strcmp(type, 'CH')
    is_RS = np.ones((1, 10))
    is_IB = np.zeros((1, 10))
    is_CH = np.zeros((1, 10))

    # N = length(type)
    N = 10

    ParamList = np.zeros((N, 9))
    # ParamList(is_RS,:) = ones(length(is_RS(is_RS)),1)*Params(1,:);
    # ParamList[is_IB,:] = ones(length(is_IB(is_IB)),1)*Params(2,:)
    # ParamList[is_CH,:] = ones(length(is_CH(is_CH)),1)*Params(3,:)
    ParamList = np.ones((10, 1)) * Params[0, :]

    b = ParamList[:, 7]
    Vr = ParamList[:, 8]

    VERTSIZE = math.floor((Tin[1] - Tin[0]) / stepT) + 1

    T = np.zeros((VERTSIZE, 1))
    T[0] = Tin[0]
    Y = np.zeros((2, VERTSIZE, N))
    Iout = np.zeros((VERTSIZE, N))
    # Iout[1, :] = Iin[int(X0[0]), int(X0[1])]
    Iout[1, :] = np.reshape(Iin(X0[0], X0[1]), (10))

    #Gradient Detector Parameters
    tauGrad = 3
    c = 3.373316578276000 * 2
    d = 547.3011736356295 * 2

    if Tin[0] == 0:
        Y[0, 0, :] = -70E-3 * np.ones((1, N))
        Y[1, 0, :] = 0
    else:
        Y[:, 0, :] = np.reshape(np.transpose(Y0, [1, 2, 0]), (2, 10))
        #CHECK STATEMENT
        # Y[:, 0, :] = np.swapaxes(Y0, 2, 0)
        # Y[1, :, :] = np.swapaxes(Y0, 0, 1)
        # random_statement = 0

    # Checking for Initial Step - DECREMENTED BY ONE (originally n = 2) (used in while loop and l8r)
    n = 1
    maxT = Tin[1]

    I1 = np.zeros((N, N))
    I2 = np.zeros((N, N))

    # Check to make sure / signs are correct
    decay1 = np.exp(-stepT / tau)
    decay2 = np.exp(-stepT / tau_s)
    decay_v = np.exp(-stepT / tau_v)
    temp = np.zeros((1, N))

    X = np.zeros((VERTSIZE, 2))
    X[0, :] = X0
    Vel = np.zeros((VERTSIZE, 0))
    wormVel = Vel_ss
    Vel[0] = wormVel

    # Weight arrays created, initial values provided from Syn_const_mat
    W65 = np.zeros((VERTSIZE, 1))
    W65[0] = Syn_const_mat[1, 5, 4]
    W98 = np.zeros((VERTSIZE, 1))
    W98[0] = Syn_const_mat[1, 8, 7]

    SpikeCount = 0
    while T[n - 1] + stepT < maxT:
        T[n] = T[n - 1] + stepT

        I1 *= decay1
        I2 *= decay2
        if np.any(temp):
            # print(Syn_is_conn[:, np.nonzero(temp)])
            temp1 = np.any(Syn_is_conn[:, temp], axis=1)

            # print("N")
            print(n)

            # print("\nTEMP")
            # print(temp)

            # print("TEMP1")
            # print(temp1)

            # print(I1.shape)
            cord_list = np.nonzero(np.reshape(temp1, (10, 1)) @ np.reshape(temp, (1, 10)))

            # print(np.reshape(temp1, (10,1)))
            xes = [i[0] for i in cord_list]
            yes = [i[1] for i in cord_list]

            # I1[np.nonzero(temp1), np.nonzero(temp)] = I1[np.nonzero(temp1), np.nonzero(temp)] + (Syn_const_mat[0, np.nonzero(temp1), np.nonzero(temp)] * Syn_const_mat[0, np.nonzero(temp1), np.nonzero(temp)])
            # print(I1[np.nonzero(temp1), np.nonzero(temp)])
            # print(I1[temp1, temp])
            # I1[temp1, temp] = I1[temp1, temp] + (Syn_const_mat[0, temp1, temp] * Syn_const_mat[0, temp1, temp])
            I1[xes, yes] = I1[xes, yes] + (Syn_const_mat[0, xes, yes] * Syn_const_mat[0, xes, yes])
            # I1 = I1 + (Syn_const_mat * Syn_const_mat)

            # I2[np.nonzero(temp1), np.nonzero(temp)] = I2[np.nonzero(temp1), np.nonzero(temp)] + (Syn_const_mat[0, np.nonzero(temp1), np.nonzero(temp)] * Syn_const_mat[1, np.nonzero(temp1), np.nonzero(temp)])
            I2[xes, yes] = I2[xes, yes] + (Syn_const_mat[0, xes, yes] * Syn_const_mat[1, xes, yes])
            # I2[:, :] = I2[:, :] + (Syn_const_mat[0, :, :] * Syn_const_mat[1, :, :])

        else:
            temp1 = temp

        # Iext = Iin[int(X[n - 1, 0]), int(X[n - 1, 1])]
        Iext = Iin(X[n - 1, 0], X[n - 1, 1])
        # print(Iext)
        # print(Iext.shape)
        # Iext = [ 0.136616701317192e-08, 0.08305e-08, -0.0396e-08, 0.060e-08, 0.0800e-08, 0., 0.0600e-08, 0.0800e-08, 0., 0.0205e-08]
        # print("Iext")
        # print(Iext)
        Iin_curr = np.transpose(Iext) + np.sum(I1 - I2, axis=1)
        # print(np.transpose(np.sum(I1 - I2, axis=1)))
        # print(Iin_curr.shape)
        # print("Iin_curr")
        # print(Iin_curr)

        # Iout[n, :] = np.transpose(Iin_curr)
        Iout[n, :] = np.reshape(Iin_curr, (10))

        # print("Iout[n, :]")
        # print(Iout[n, :])

        # print(Y[:, n-1, :])
        #PERMUTE?
        K1 = AEF_dY_dt(Iin_curr, np.swapaxes(Y[:, n - 1, :], 0, 1), ParamList)

        # K1 = AEF_dY_dt(Iin_curr, np.reshape(Y[n - 1, :, :], [10, 2]), ParamList)

        # K1 = np.transpose(K1, [2, 0, 1])
        # print(K1.shape)
        # # K1 = np.reshape(np.swapaxes(K1, 0, 1), (2, 10))
        # # K1.reshape((10, 2))
        # print(K1.shape)
        # print(K1)

        # print(Y[:, n, :].shape)
        # print(Y[:, n-1, :].shape)
        # print(K1.shape)

        # Y[:, n, :] = Y[:, n - 1, :] + stepT * K1
        Y[:, n, :] = Y[:, n - 1, :] + stepT * np.reshape(np.swapaxes(K1, 0, 1), (2, 10))
        #int( n prp.reshape(K1, (2, 10)))
        # print(K1)
        # print("Y[:, n, :]")
        # print(Y[:, n, :])
        # print("\nY[:, n - 1, :]")
        # print(Y[:, n - 1, :])
        # print("\nK1")
        # print(np.reshape(np.swapaxes(K1, 0, 1), (2, 10)))
        # print(stepT * np.reshape(np.swapaxes(K1, 0, 1), (2, 10)))
        # print("\nStepT")
        # print(stepT)


        # temp is 1x10 logical array
        # print("Y[0, n, :]")
        # print(Y[0, n, :])
        temp = (Y[0, n, :] > 0)
        # print("temp")
        # print(temp)

        #CHECK STATEMENT WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo
        Y[0, n, np.nonzero(temp)] = 0
        # print(Y[0, n, :])

        # Another 1x10 logical array
        # print(Y[0, n - 1, :])
        tempreset = (Y[0, n - 1, :] == 0)
        # Y[n, tempreset, 0] = np.transpose(Vr[tempreset, 0])
        # print(str(n) + "ITERATION ***")
        # if(np.any(tempreset)):

             # print(tempreset)
             # print(n)
        # print("*****")

        # print("\nTEMPRESET")
        # print(tempreset)
        # CHECK STATEMENT SCALAR MULTIPLICATION - I FILLED WITH ZEROS (Not None or NaN) OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        # Y[0, n, tempreset] = Vr[0] * Vr[tempreset]

        Y[0, n, np.nonzero(tempreset)] = Vr[0]

        # print("Y[0, n, :]")
        # print(Y[0, n, :])
        # Y[0, n , tempreset] = Vr[tempreset, 0]

        # Y[n, tempreset, 1] = Y[n - 1, tempreset, 1] + np.transpose(b[tempreset, 0])
        # Y[1, n, tempreset] = Y[1, n - 1, tempreset] + b[tempreset, 0]
        # print(b.shape)

        # REMOVED + B * TEMPRESET
        Y[1, n, np.nonzero(tempreset)] = Y[1, n - 1, np.nonzero(tempreset)]

        # Modify     Weight
        K = Syn_const_mat[1, :, :]
        K[5, 4] = K[5, 4] - (d + K[5, 4]) * stepT / tauGrad
        K[8, 7] = K[8, 7] - (d + K[8, 7]) * stepT / tauGrad

        if tempreset[4]:
            K[5, 4] = K[5, 4] + c / tauGrad
        if tempreset[7]:
            K[0, 7] = K[0, 7] + c / tauGrad
        Syn_const_mat[1, :, :] = K
        W65[n] = K[5, 4]
        W98[n] = K[8, 7]


        # TURNING
        if tempreset[5]:
            theta0 = theta0 - np.pi / 180 * 5
        if tempreset[8]:
            theta0 = theta0 + np.pi / 180 * 5
        if tempreset[9]:
            theta0 = theta0 + np.pi * (np.random.rand() - 0.5)

        # MOVEMENT
        if tempreset[1] or tempreset[2]:
            wormVel = wormVel + pump_vel
        else:
            wormVel = (wormVel - Vel_ss) * decay_v + Vel_ss
        Vel[n] = wormVel

        rstep = wormVel * stepT
        X[n, 0] = X[n - 1, 0] + rstep * math.cos(theta0)
        X[n, 1] = X[n - 1, 1] + rstep * math.sin(theta0)

        #Count total number of spikes

        # Added Conditional
        if tempreset.ndim < 2:
            SpikeCount = SpikeCount + np.sum(tempreset)
        else:
            SpikeCount = SpikeCount + np.sum(tempreset, axis=1)

        # if(n % 100000 == 0):
        #     plt.scatter(X[:,0], X[:, 1])
        #     plt.plot(X[:,0], X[:, 1])
        
        #     plt.show()
        # print(n)
        # print(tempreset)
        n = n + 1

    print("GOING TO BE RETURNED SPIKECOUNT: " + str(SpikeCount))
    print("DONE!!")
    return (T, Y, Iout, X, SpikeCount, Vel, W65, W98)
