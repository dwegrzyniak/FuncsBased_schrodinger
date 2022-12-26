import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from IPython.display import HTML, Image
from matplotlib.animation import FuncAnimation

data_file = open("./dane/output_bez_15.txt", "r+")

V = []
x = []
psi_norm = []
psi_real = []
psi_imag = []

x_line = data_file.readline()
data_file.readline()
psi_norm_line = data_file.readline()
data_file.readline()
psi_real_line = data_file.readline()
data_file.readline()
psi_imag_line = data_file.readline()
data_file.readline()
V_line = data_file.readline()
data_file.readline()

while x_line:
    x.append(x_line.split(',')[:-1])
    psi_norm.append(psi_norm_line.split(',')[:-1])
    psi_real.append(psi_real_line.split(',')[:-1])
    psi_imag.append(psi_imag_line.split(',')[:-1])
    V.append(V_line.split(',')[:-1])
    
    x[-1] = [float(tmp) for tmp in x[-1]]
    psi_norm[-1] = [float(tmp) for tmp in psi_norm[-1]]
    psi_real[-1] = [float(tmp) for tmp in psi_real[-1]]
    psi_imag[-1] = [float(tmp) for tmp in psi_imag[-1]]
    V[-1] = [float(tmp) for tmp in V[-1]]
    
    x_line = data_file.readline()
    data_file.readline()
    psi_norm_line = data_file.readline()
    data_file.readline()
    psi_real_line = data_file.readline()
    data_file.readline()
    psi_imag_line = data_file.readline()
    data_file.readline()
    V_line = data_file.readline()
    data_file.readline()
    
def animate(i):
    ax.clear()
    
    ax.set_xlabel('x')
    ax.set_ylabel('psi')
    ax.set_title('psi(x)')
    
    plt.ylim((0, 1))
    ax.plot(x[i], psi_real[i], label = "psi real")
    ax.plot(x[i], psi_imag[i], label = "psi imag")
    ax.plot(x[i], psi_norm[i], label = "psi norm")
    ax.plot(x[i], V[i], label = "V")
    
    ax.legend()
    
fig, ax = plt.subplots()
ax.plot(x[0], psi_norm[0])

start, end = ax.get_xlim()
ax.xaxis.set_ticks(np.arange(start, end, 10))



ani = FuncAnimation(fig, animate, frames=len(x), interval=50)
#plt.show()
ani.save("bez_15.gif")


