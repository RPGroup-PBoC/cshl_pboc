# the usual imports
import numpy as np
import matplotlib.pyplot as plt

# In this tutorial, we will investigate the dynamics of mRNA expression from an
# unregulated (constitutive) promoter.

# As we derived in class, the differential equation we wish to integrate is
#               dm/dt = r - y*m
# where r is the production rate, y is the degradation rate, and m is the
# number of mRNAs in the cell.
# To start, we will define some parameter values.

# Assign values to model parameters.
r = 1 # mRNA production rate in 1/minute
gamma = 1/3 # mRNA decay rate in 1/minute

# Define parameters for our Euler integration
time = 20 # minute
dt = 0.1 # time step in min
num_steps = int(time/dt)
time_vec = np.linspace(0, time, num_steps)

# Set up storage for mRNA copy #
m_t = np.empty_like(time_vec)
# initial condition
m_t[0] = 0

# Loop over time points and update mRNA copy number
for i in range(1, num_steps):
    m_t[i] = m_t[i-1] + r*dt - gamma*dt*m_t[i-1]

# Plot results of integration
plt.figure()
plt.plot(time_vec, m_t, '-', label='Euler integration')
plt.xlabel('time (min)')
plt.ylabel('mRNA copy #')
# Plot a horizontal line to show steady-state value
plt.hlines(r/gamma, xmin=0, xmax=time, label='steady state')
plt.legend()
plt.show()
