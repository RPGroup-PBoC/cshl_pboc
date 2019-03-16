# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Import our class utilities
import pboc_utils as pboc


# In this script, we will examine the process of diffusion by performing
# numerical integration of a master equation. As we discussed in lecture
# the diffusion master equation can be written as
#
# P(x, t + dt) = P(x, t) + k*dt*P(x-1, t-1) + k*dt*P(x+1,t-1) - 2*k*dt*P(x, t-1), #
# where P(x,t) is the probability of finding a particle in box x at time
# t, k is the diffusion rate of the particle, and dt is the time difference.
# We'll examine the behavior of this master equation under a number of
# different situtations.

# To begin, let's  examine a relatively simple problem. Imagine that you have
# a single particle in the middle (x) of a one dimensional cell. What is the
# probaility of finding the particle in a position x + dx at time t + dt? To
# perform this integration, we will use the Euler-forward mehtod where we
# numerically compute the probability  at each position using the knowldege of
# the probabilities in the previous time step. To start, we'll define some
# parameters. D = 10 um^2/sec is order of magnitude for protein in cytoplasm.
D = 10 # diffusion constant in um^2/sec
dx = 0.01 # lattice spacing in microns
k = D/dx**2 # hopping probability per unit time (sec^-1)
dt = 1/(10*k) # time step. note k*dt should be small compared to 1
                # for numerical stability

# To do this in the 'continuous' limit numerically, we'll set up an array with
# many boxes for the particle to diffuse into. This means that we will have to
# be careful of boundaries (the first box and the last) when computing the
# probabilities. Since we'll do this for each example in this script, let's
# write the master equation as a function.
def master_eq(prob, k, dt):
    """
    Computes the master equation for diffusion in one dimension

    Parameters
    ----------
    prob : 2d-array
        Array in which the probabilities will be calculated. This should be in
        the shape of N boxes by M time points. This should have a preset
        initial condition.

    k : float
        Diffusion rate of the particles in units of 1/s
    dt : float
        Time step for the integration.

    Returns
    -------
    prob : 2d-array
        The probability vector supplied populated with the calculated
        probabilities.
    """

    # We'll first figure out the number of boxes and the number of time
    # steps.
    num_boxes, time_points = np.shape(prob)

    # We need to integrate over each time step. We'll start at the second
    # time point since we will already have the initial condition set.
    for t in range(1, time_points):

        # Now we will deal with the boundary conditions. At the first box,
        # there is no way to calculate the probability of the previous box.
        # We will have this boundary be reflecting, meaning that nothing can
        # diffuse past this point.
        prob[0, t] = prob[0, t-1] - k * dt * prob[0, t-1] \
                        + k * dt * prob[1, t-1]

        # We will also set the boundary of the last box to be reflecting.
        prob[-1, t] = prob[-1, t-1] - k * dt * prob[-1, t-1] \
                        + k * dt * prob[-2, t-1]

        # With the boundary conditions met, we can now do the integration over
        # the other boxes.
        for x in range(1, num_boxes - 1):
            prob[x, t] = prob[x, t-1] + k * dt * prob[x+1, t-1] \
                        - 2 * k * dt * prob[x, t-1] + k * dt * prob[x-1, t-1]
    return prob

# With this function in place, let's generate the vectors needed for our
# infite plane diffusion.
num_boxes = 100
time_points = 100
time_conversion = np.linspace(0, time_points*dt, time_points)
prob_inf = np.zeros((num_boxes, time_points))

# Now we'll set the starting position. We'll put our particle in the middle
# (origin) and give it a probability of one.
prob_inf[49, 0] = 1.0   # Remember that Python indexing begins at 0!

# Now we can simply pass this to our function!
prob_inf = master_eq(prob_inf, k, dt)

# And that's it! Let's look at the diffusion as a three dimensional bar plot.
# For simplicity, we'll use the 'bar3' function provided in the class utilities
fig1 = pboc.bar3(prob_inf, xlabel='time (sec)', ylabel='box number',
                 zlabel='probability', bin_step=3, y_vec=time_conversion)

# As we pan around in the above figure, we can see that the result makes sense
# intuitive sense. As time goes on, the particles will diffuse outward with
# the highest probability in the middle of the distribution How does this
# change when we put the diffusing particles in a finite box? Let's take a
# look. Redefine system parameters for so particles can reach walls.
num_boxes = 15
time_points = 200
time_conversion = np.linspace(0, time_points*dt, time_points)
prob_box = np.zeros((num_boxes, time_points))

# For fun, we'll set the initial condition as the very corner of the box.
prob_box[0, 0] = 1.0

# Let 'er rip!
prob_box = master_eq(prob_box, k, dt)

# Again, let's plot this in three dimensions
fig2 = pboc.bar3(prob_box, xlabel='time (sec)', ylabel='box number',
                zlabel='probability', bin_step=3, y_vec=time_conversion)

# This matches our physical intuition. Notice now that the probability is not
# zero anywhere in the box. it has filled up the box and the distribution has
# now become uniform.

# Let's look at a more complicated case. Fluorescence Recovery After
# Photobleaching is a common method in biology of measuring the
# diffussivity of a particular molecule in a cell. This includes labeling
# your molecule of interest with a fluorescent molecule and allowing it to
# unformly diffuse in the cell. By photobleaching a small section (breaking
# the fluorophores so they are no longer fluorescent) and watching how rapidly
# the fluorescence diffuses into the hole, you can measure the mobility of
# your molecule. Let's work this into our numerical integration. We'll
# generate a two-dimensional vector as we've done above  except the initial
# condition will have a high probability except for a few boxes in the
# middle where there is no probability at all.
num_boxes = 15
prob_frap = np.zeros((num_boxes, time_points))

# Set the initial condition. the oddball definition is to preserve
# normalization of probability after bleaching
prob_frap[:, 0] = 1 / (num_boxes - 7)
# We'll bleach six boxes in this example.
prob_frap[4:11, 0] = 0   # In Python, we can index from [start:stop)

# Perform the integration.
prob_frap = master_eq(prob_frap, k, dt)

# Let's take a look.
fig3 = pboc.bar3(prob_frap, xlabel='time (sec)', ylabel='box_number',
                    zlabel='probability', bin_step=3, y_vec=time_conversion)
plt.show()

# Notice how the probability within the hole filled, but it never returned to
# the initial probability. This is because we actively removed molecules by
# bleaching them.
