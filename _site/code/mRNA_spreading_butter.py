# Import the usual modules.
import numpy as np
import matplotlib.pyplot as plt

# Import the course utilities.
import pboc_utils as pboc

# In this tutorial, we will expand upon our numerical integration of the mean
# mRNA copy number by solving for the complete mRNA copy number distribution as
# a function of time.

# As is dervied in class, the change in the probability of observing a given
# mRNA copy number can be described as
#
#  dP(m,t)/dt = rP(m -1 ,t) + y(m + 1)P(m+1, t) - rP(m, t) - ymP(m, t)
#
# where m is the number of mRNAs, t is the timepoint, r is the productionj
# rate of mRNA, y is the degradation rate of the mRNA, and P(m,t) is the
# probability of mRNAs at time t.

# First define parameters for problem
r = 2 # mRNA production rate in min^(-1)
gamma = 1/3 # mRNA decay rate in min^(-1)
time = 20 # Total simulation time to run, in units of min
dt = 0.05 # Time step in units of min
num_steps = int(time/dt) # total # of time steps
upper_bound = 20 # max copy # to simulate

# We want to keep track of the probability of each number of mRNAs at each
# individual timepoint. To do so, we'll keep them in a two-dimensional array
# in which each row corresponds to a possible mRNA copy number and
# each column is a specific time point.

# Construct a 2D storage vector
prob = np.zeros((upper_bound+1, num_steps))

# Now just set initial condition and start integration
prob[0, 0] = 1

# Run the Euler integration! As with diffusion we should consider the boundaries separately
for t in range(1, int(num_steps)):
    # First handle the m=0 case
    prob[0, t] = prob[0, t-1] + gamma*dt*prob[1, t-1] - r*dt*prob[0, t-1]
    # now handle the upper boundary condition
    prob[upper_bound, t] = prob[upper_bound, t-1] \
                            + r*dt*prob[upper_bound-1, t-1] \
                            - gamma*dt*upper_bound*prob[upper_bound, t-1]
    # loop over mRNA copy number to cover everything else
    for m in range(1, upper_bound):
        # Compute p(m,t)
        prob[m, t] = prob[m, t-1] + r*dt*prob[m-1, t-1] \
                    + gamma*(m+1)*dt*prob[m+1, t-1] - gamma*m*dt*prob[m, t-1] \
                    - r*dt*prob[m, t-1]

# We can now show the distribution of copy numbers at the beginning and end
# of the integration.

# Make a vector of mRNAs.
mRNA_vec = np.arange(0, upper_bound + 1, 1)

# Make a new figure with two plotting axes.
fig, ax = plt.subplots(2, 1)

# Make a bar plot of the time 0 distribution
ax[0].bar(mRNA_vec, prob[:, 0], width=1)
ax[0].set_ylabel('probability')
ax[0].set_title('t=0 min')
# and another at the end of the sim
ax[1].bar(mRNA_vec, prob[:, -1], width=1)
ax[1].set_title('t=5 min')
ax[1].set_ylabel('probability')
ax[1].set_xlabel('number of mRNAs')
plt.show()

# To get a better sense of how the distribution changes with time, let's make
# a bar plot for each point in time and plot it in three dimensions. To save
# ourselves a bit of headache with syntax, we'll use the function 'bar3' in
# the course utilities file to generate this plot.
# Generate a time vector.
time_vec = np.linspace(0, time, num_steps)

pboc.bar3(prob, xlabel='time', ylabel='number of mRNA',
          zlabel='probability',  y_vec=time_vec, bin_step=5)
plt.show()
