import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
from cycler import cycler

g = 9.81                    #Gravedad
w = 1.448                   #Wheelbase
varepsilon = 26.8 * np.pi/180   #Caster angle
a_n = 0.105                 #Trail
m_0 = 195                   #Masa moto
b_0 = 0.722                 #Posición CoM
h_0 = 0.482                 #Posición CoM
I_0xx = 13.5                #Momentos inercia
I_0xz = 3                   #Momentos inercia
I_0zz = 55                  #Momentos inercia
m = 270                     #Masa total (c/ piloto)
b = 0.688                   #CoM pos (c/ piloto)
h = 0.64                    #CoM pos (c/ piloto)
I_xx = 35.5                 #Momentos inercia (c/ piloto)
I_xz = -1.7                 #Momentos inercia (c/ piloto)
I_zz = 59.3                 #Momentos inercia (c/ piloto)
m_f = 34                    #Masa eje delantero
e_f = 0.025                 #CoM delantero 
h_f = 0.6                   #CoM delantero
I_fzz = 0.83                #Momento incercia frontal
I_omega_f = 0.6             #Inercia rueda delantera (giro)
I_omega_r = 0.8             #Inercia rueda trasera (giro)
c_delta = 1                 #Amortiguador direccion

R_f = 0.294                 #Radio rueda delantera
R_r = 0.299                 #Radio rueda trasera
rho_f = 0.064               #Radio de la cross-section delantera
rho_r = 0.078               #Radio de la cross-section trasera
k_alpha_f = 16              #Rigidez de viraje normalizada
k_alpha_r = 14.5            #Rigidez de viraje normalizada
k_phi_f = 0.85              #Rigidez camber normalizada
k_phi_r = 0.95              #Rigidez camber normalizada
k_a_f = -0.2                #Rigidez auto alineante normalizada
k_a_r = -0.2                #Rigidez auto alineante normalizada
k_t_f = 0.015               #Rigidez de torsión
k_t_r = 0.018               #Rigidez de torsión
k_l_f = 160000              #Rigidez estructural transversal
k_l_r = 140000              #Rigidez estructural transversal

l_beta = 0.67               #Punto de flexión horquilla
l_b = 0.67 ## podria ser    #BUSCAR
k_beta = 38000              #Rigidez de flexión
m_b = 18                    #Masa de flexión (Masa de horquilla?)
e_b = 0                     #CoM horquilla
h_b = 0.35                  #CoM horquilla
I_bxx = 0.8                 #Momento de inercia

CdA = 0.467                 #Factor arrastre aerodinámico
h_A = 0.35                  #Altura centro aerodinámico

a_x = 0                     #Aceleración
omega_f_dot = 0             #Aceleración angular rueda
omega_r_dot = 0             #Aceleración angular rueda

N_r_0 = ((w - b)/w) * m * g
N_f_0 = (b/w) * m * g

#Asumido---------------------------------

k_gamma_f = 1
k_gamma_r = 1

k_y_f = k_l_r
k_y_r = k_l_f


I_bzz = I_bxx * 0.1       #
x_f    = 0.5                #
x_b    = 0.155              #

#---------------------------------------

b_f = w + (x_f + a_n - h_f * np.sin(varepsilon))/np.cos(varepsilon)
b_b = w + (x_b + a_n - h_b * np.sin(varepsilon))/np.cos(varepsilon)
z_b = l_b + ((a_n + x_b) * np.sin(varepsilon) - h_b)/np.cos(varepsilon)


#--------------------------------------------------------------------
E = np.zeros((10,10))
E[0,0] = m
E[0,1] = m*b
E[0,2] = m*h
E[0,3] = m_f * e_f
E[0,4] = -m_b*z_b

E[1,1] = m*b**2 + I_zz
E[1,2] = m*b*h - I_xz
E[1,3] = m_f*e_f*b_f + I_fzz * np.cos(varepsilon)
E[1,4] = -m_b*z_b*b_b -I_bxx * np.sin(varepsilon)

E[2,2] = m*h**2 + I_xx
E[2,3] = m_f*e_f*h_f + I_fzz * np.sin(varepsilon)
E[2,4] = -m_b*h_b*z_b + I_bxx * np.cos(varepsilon)

E[3,3] = m_f*e_f**2 + I_fzz
E[3,4] = -m_b*e_b*z_b

E[4,4] = m_b*z_b**2 + I_bzz

E[5,5] = k_alpha_r

E[6,6] = k_alpha_f

E[7,7] = 1
E[8,8] = 1
E[9,9] = 1

for i in range(1, 10):
    for j in range(i):
        E[i, j] = E[j, i]
        
E2 = np.linalg.inv(E)

A_total = []
A_total2 = []

Vx_rango = np.linspace(5,61,75)

for i in Vx_rango:
    V_x = i
    X_f = 0
    omega_f = V_x/R_f
    omega_r = V_x/R_r

    rho_air = 1.2041
    F_ad = 0.5 * rho_air * CdA * V_x**2
    X_r = F_ad

    N_r = N_r_0 + (h_A/w) * F_ad
    N_f = N_f_0 - (h_A/w) * F_ad


    A = np.zeros((10,10))
    A[0,1] = -m * V_x
    A[0,5] = k_alpha_r * N_r
    A[0,6] = k_alpha_f * N_f
    A[0,7] = k_gamma_f * N_f + k_gamma_r * N_r
    A[0,8] = X_f * np.cos(varepsilon) + N_f * k_gamma_f * np.sin(varepsilon)
    A[0,9] = -X_f * np.sin(varepsilon) + N_f * k_gamma_f * np.cos(varepsilon)

    A[1,1] = -m * b * V_x
    A[1,2] = I_omega_r * omega_r + I_omega_f * omega_f
    A[1,3] = I_omega_f * omega_f * np.sin(varepsilon)
    A[1,4] = I_omega_f * omega_f * np.cos(varepsilon)
    A[1,5] = k_a_r * N_r
    A[1,6] = k_alpha_f * w * N_f + k_a_f * N_f
    A[1,7] = k_t_r * N_r + (k_t_f + w * k_gamma_f) * N_f + h * m * a_x + h_A * F_ad + I_omega_r * omega_r_dot + I_omega_f * omega_f_dot - X_f * rho_f - X_r * rho_r
    A[1,8] = k_t_f * N_f * np.sin(varepsilon) + (w * np.cos(varepsilon) - rho_f * np.sin(varepsilon) + a_n) * X_f + m_f * e_f * a_x + I_omega_f * omega_f_dot * np.sin(varepsilon) + N_f * k_gamma_f * w * np.sin(varepsilon)
    A[1,9] = (l_b - rho_f * np.cos(varepsilon) - w * np.sin(varepsilon)) * X_f + I_omega_f * omega_f_dot * np.cos(varepsilon) - m_b * z_b * a_x + k_t_f * N_f * np.cos(varepsilon) + N_f * k_gamma_f * w * np.cos(varepsilon)

    A[2,1] = -m * h * V_x - I_omega_r * omega_r - I_omega_f * omega_f
    A[2,3] = -I_omega_f * omega_f * np.cos(varepsilon)
    A[2,4] = I_omega_f * omega_f * np.sin(varepsilon)
    A[2,7] = m * g * h - rho_f * N_f - rho_r * N_r
    A[2,8] = (a_n - rho_f * np.sin(varepsilon)) * N_f + m_f * e_f * g - I_omega_f * omega_f_dot * np.cos(varepsilon)
    A[2,9] = (l_b - rho_f * np.cos(varepsilon)) * N_f - m_b * z_b * g + I_omega_f * omega_f_dot * np.sin(varepsilon)

    A[3,1] = -m_f * e_f * V_x - I_omega_f * omega_f * np.sin(varepsilon)
    A[3,2] = I_omega_f * omega_f * np.cos(varepsilon)
    A[3,3] = -c_delta
    A[3,4] = I_omega_f * omega_f
    A[3,6] = (k_a_f  * np.cos(varepsilon) - a_n * k_alpha_f) * N_f
    A[3,7] = (a_n * (1 - k_gamma_f) - rho_f * np.sin(varepsilon)) * N_f - rho_f * X_f * np.cos(varepsilon) + m_f * e_f * g + N_f * k_t_f * np.cos(varepsilon)
    A[3,8] = k_gamma_f * a_n * N_f * np.sin(varepsilon)  + A[3,7] * np.sin(varepsilon) - N_f * a_n * k_gamma_f * np.sin(varepsilon)
    A[3,9] = (k_t_f * np.cos(varepsilon)**2 - rho_f * np.sin(varepsilon) * np.cos(varepsilon) + l_b * np.sin(varepsilon)) * N_f - m_b * z_b * (g * np.sin(varepsilon) + a_x * np.cos(varepsilon)) + (a_n * np.sin(varepsilon) - rho_f * np.cos(varepsilon)**2 + l_b * np.cos(varepsilon)) * X_f + I_omega_f * omega_f_dot  - N_f * a_n * k_gamma_f * np.cos(varepsilon)

    A[4,1] = m_b * z_b * V_x - I_omega_f * omega_f *np.cos(varepsilon)
    A[4,2] = -I_omega_f * omega_f
    A[4,3] = -I_omega_f * omega_f * np.sin(varepsilon)
    A[4,4] = -I_omega_f * omega_f * np.cos(varepsilon)
    A[4,6] = -(k_a_f * np.sin(varepsilon) + l_b * k_alpha_f) * N_f
    A[4,7] = ((1 - k_gamma_f) * l_b - k_t_f * np.sin(varepsilon) - rho_f * np.cos(varepsilon)) * N_f + rho_f * X_f * np.sin(varepsilon) - m_b * z_b * g
    A[4,8] =  - m_b * z_b * a_x * np.cos(varepsilon) + A[4,7] * np.sin(varepsilon)
    A[4,9] =  m_b * z_b * a_x * np.sin(varepsilon) + A[4,7] * np.cos(varepsilon) - k_beta

    A[5,0] = -k_y_r / N_r
    A[5,2] = 1 - k_gamma_r
    A[5,5] = -V_x * k_y_r / N_r

    A[6,0] = -k_y_f/N_f
    A[6,1] = (- k_y_f) / N_f * w
    A[6,2] = 1 - k_gamma_f
    A[6,3] = (1 - k_gamma_f) * np.sin(varepsilon) + a_n * k_y_f/N_f
    A[6,4] = (1 - k_gamma_f) * np.cos(varepsilon) + l_b * k_y_f/N_f + (R_f * k_l_f * np.cos(varepsilon))/N_f
    A[6,6] = -(V_x * k_y_f) / N_f
    A[6,8] = V_x*np.cos(varepsilon) * k_y_f/N_f
    A[6,9] = - V_x*np.sin(varepsilon) * k_y_f/N_f

    A[7,2] = 1
    A[8,3] = 1
    A[9,4] = 1

    A_total.append(A)
    A_total2.append(E2@A)

x_val_total = []
x_val_real  = []

x_val_total2 = []
x_val_real2 = []


# for k in A_total:
#     u = np.linalg.eigvals(np.dot(E2, k))
#     plt.scatter(np.real(u), np.imag(u), facecolors="none", edgecolors="black", s=10)


# plt.xlabel('Parte real')
# plt.ylabel('Parte imaginaria')
# plt.xlim([-30, 10])
# plt.ylim([0, 65])
# plt.text(-7, 45, 'Wobble', size=10)
# plt.text(-7, 20, 'Weave', size=10)
# plt.text(-7, 1, 'Capsize', size=10)
# plt.grid()
# plt.show()



fig, ax1 = plt.subplots()

# Scatter plot on the first y-axis
sizecounter = 0
for k in A_total:
    u = np.linalg.eigvals(np.dot(E2, k))
    ax1.scatter(np.real(u), np.imag(u), facecolors="none", edgecolors="black", s=10+sizecounter)
    sizecounter += 1

ax1.axvline(x = 0, color = 'k', linestyle = '-')
ax1.set_xlabel('Parte real [1/s]')
ax1.set_ylabel('Parte imaginaria [1/s]')
ax1.set_xlim([-30, 10])
ax1.set_ylim([0, 65])
txtwobble = (r'$Wobble$')
txtweave = (r'$Weave$')
txtcapsize = (r'$Capsize$')
ax1.text(-7, 57, txtwobble, size=15, math_fontfamily='cm') #style='italic')
ax1.text(-7, 20, txtweave, size=15, math_fontfamily='cm') #style='italic')
ax1.text(-10, 1, txtcapsize, size=15, math_fontfamily='cm') #style='italic')
ax1.text(-25, 60, 'Estable', size=15, fontweight='bold')
ax1.text(1, 60, 'Inestable', size=15, fontweight='bold')
ax1.fill_betweenx(np.arange(0, 66), 0, 10, where=(np.arange(0, 66) >= 0), color='lightgray', alpha=0.5)
ax1.grid()


# Create a second y-axis
ax2 = ax1.twinx()
ax2.set_ylim([0, 10])  # Set the desired limits for the second y-axis
ax2.set_ylabel('Frecuencia [Hz]')


plt.show()


for j in A_total2:
    x_val = eigvals(j)
    x_val_total.append(x_val)
    x_val_real.append(x_val.real)

# default_cycler = (cycler(color=['r', 'g', 'b',]) +
#                   cycler(linestyle=['-', '--', ':']))

# default_cycler = (cycler(color=['k', 'k', 'k',]) +
#                   cycler(linestyle=['-', '--', ':']))


# plt.rc('lines', linewidth=2)  
# plt.rc('axes', prop_cycle=default_cycler) 

plt.plot(Vx_rango, x_val_real, 'ok', markersize=3)#, fillstyle="none")
plt.hlines(0, 0, 60, colors='k', linestyles='solid')
plt.xlabel('Velocidad [m/s]')
plt.ylabel('Parte real')
plt.xlim([0, 60])
plt.ylim([-10, 5])
#plt.legend(['Weave', 'Wobble', 'Capsize'])
plt.grid()
plt.show()
