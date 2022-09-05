import numpy as np
from matplotlib import pyplot as plt
from scipy.fft import fft2, ifft2
# Author: Elvis do A. Soares
# Github: @elvissoares
# Date: 2022-09-05

"""
 The python script to solve the Phase Field Crystal 
 equation using an semi-implicit pseudo-spectral algorithm
"""
# PFC model constants
r = -0.25
M = 1.0 # mobility

# Size of the system
N = 2**8 # 2**8 = 256
L = 16*np.pi
x = np.linspace(0,L,N)
dx = x[1]-x[0]

# The time step definition
dt = 0.1
T = 3000
Nsteps = int(T/dt)
Nframes = Nsteps//100 #frames to the output

# The vector of variables
n_hat = np.empty((N,N), dtype=np.complex64)
n = np.empty((Nframes,N,N), dtype=np.float32)

# Initial condition
rng = np.random.default_rng(12345)
noise = 0.02*0.285
n0 = -0.285
n[0] = n0 + noise*rng.standard_normal(n[0].shape)

# The Fourier variables and dealising vector
kx = np.fft.fftfreq(N, d=dx)*2*np.pi
k = np.array(np.meshgrid(kx , kx ,indexing ='ij'), dtype=np.float32)
k2 = np.sum(k*k,axis=0, dtype=np.float32)
kmax_dealias = kx.max()*2.0/3.0 # The Nyquist mode
dealias = np.array((np.abs(k[0]) < kmax_dealias )*(np.abs(k[1]) < kmax_dealias ),dtype =bool)

# The linear terms of PDE
L_operator = -M*k2*(k2**2-2*k2+1+r)

# The non-linear terms of PDE (with dealising)
def Noperator_func(n):
    return -(k2*M*fft2(n**3))*dealias

# auxiliary variables
lineardenominator_hat = 1.0/(1.0-dt*L_operator) # can be calculated once
Noperator_hat = n_hat.copy()

# time evolution loop
nn = n[0].copy()
n_hat[:] = fft2(n[0]) # FT initial condition
for i in range(1,Nsteps):
    Noperator_hat[:] = Noperator_func(nn) # calculate the non-linear term
    n_hat[:] = (n_hat+dt*Noperator_hat)*lineardenominator_hat # updating in time
    nn[:] = ifft2(n_hat).real # IFT to next step
    if (i//100): 
        n[i//100] = nn
    
N0 = n[0].sum()*dx**2/L**2
Nlast = n[-1].sum()*dx**2/L**2

print('The relative difference between the first and last frame mean densities: ',np.abs(Nlast/N0-1))

# plot the last frame
plt.imshow(n[-1],cmap='viridis', vmin=-0.6, vmax=0.4)
plt.colorbar(label=r'$n(x,y)$', shrink=0.8)
plt.title('$n_0=%.3f$'% n0)
plt.savefig('pfc2d-crystal.png')
plt.show()

# plot the total free-energy
F = np.zeros(Nframes)
t = np.linspace(0.0,Nsteps*dt,Nframes)
lapn = np.empty((N,N), dtype=np.float32)
laplapn = np.empty((N,N), dtype=np.float32)
for i in range(t.size):
    n_hat[:] = fft2(n[i])
    lapn[:] = ifft2(-k2*n_hat).real
    laplapn[:] = ifft2(k2**2*n_hat).real
    F[i] = np.sum(n[i]*(lapn+0.5*laplapn)+0.5*(1+r)*n[i]**2 + 0.25*n[i]**4)*dx**2
plt.xlim(0,3000)
plt.ylim(0.03,0.034)
plt.plot(t,F/L**2,'k')
plt.xlabel('$t$')
plt.ylabel('$F[n(t)]/L^2$')
plt.text(2000,0.032,'$n_0=-0.285$')
plt.savefig('pfc2d-freenergy-crystal.png')
plt.show()

# Export animation in gif
from matplotlib import animation
from matplotlib.animation import PillowWriter

fig, ax = plt.subplots(1,1,figsize=(5,4))
im = ax.imshow(n[0],cmap='viridis', vmin=-0.6, vmax=0.4)
cb = fig.colorbar(im,ax=ax, label=r'$n(x,y)$', shrink=0.8)
tx = ax.text(200,20,'t={:.1f}'.format(0.0),
         bbox=dict(boxstyle="round",ec='white',fc='white'))
ax.set_title(r'$n_0=%.3f$'% n0)

def animate(i):
    im.set_data(n[i])
    tx.set_text('t={:.0f}'.format(i*dt*100))
    return fig,

ani = animation.FuncAnimation(fig, animate, frames= Nframes,
                               interval = 50)
ani.save('pfc2d-crystal.gif',writer='pillow',fps=10,dpi=100)

