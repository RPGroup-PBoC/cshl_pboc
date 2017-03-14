# the usual imports
import numpy as np
import matplotlib.pyplot as plt

# for reading comma separated data files
import pandas as pd

# In this exercise we analyze mRNA copy # data. The data comes from smFISH
# experiments and consists of the probability distribution
# over mRNA copy #. We will compare the distributions to the Poisson
# distribution expected for an unregulated constitutive promoter.

# First we need to load data from provided comma-separated data files.
# To do this we'll use a package called pandas. pandas is extremely powerful
# and widely used for data science in python, but we will have no need of it
# beyond this tutorial. If you're curious about the syntax or would like
# to learn more, we refer you to the internet.
df_mdn1 = pd.read_csv('/Users/muir/datasets/2017cshl_pboc/MDN1.csv')
df_pdr5 = pd.read_csv('/Users/muir/datasets/2017cshl_pboc/PDR5.csv')

# To avoid delving into pandas, let us simply extract the mRNA
# probability distributions and error bars to numpy arrays.
# Note mdn1_prob and pdr5_prob are 1D arrays; the i-th element contains
# the probability of having i copies of the mRNA in the cell.
# *_counts enumerate possible mRNA copy numbers, determined by the upper
# extent of the data.
# first do mdn1...
mdn1_prob = np.array(df_mdn1['Probability'])
mdn1_errbar = np.array(df_mdn1['Error in probability'])
mdn1_counts = np.arange(len(mdn1_prob))
# ...then same thing for pdr5
pdr5_prob = np.array(df_pdr5['Probability'])
pdr5_errbar = np.array(df_pdr5['Error in probability'])
pdr5_counts = np.arange(len(pdr5_prob))

# It's always a good idea to plot the data to see what we're dealing with.
# The plt.errorbar syntax takes x data, then y data, then keyword args,
# including the array describing errorbars for each datum.
plt.figure(1)
plt.errorbar(mdn1_counts, mdn1_prob, yerr=mdn1_errbar,
                fmt='go', label='mdn1')
plt.xlabel('mRNA copy number')
plt.ylabel('probability')
plt.legend()
# Here we're creating separate figs so we can overlay the theory curves
# on the data, once we've calculated them below, without overcrowding.
plt.figure(2)
plt.errorbar(pdr5_counts, pdr5_prob, yerr=pdr5_errbar,
                fmt='bo', label='pdr5')
plt.xlabel('mRNA copy number')
plt.ylabel('Probability')
plt.legend()

# By eye, mdn1 could plausibly follow a Poisson dist, but pdr5 appears to
# have a much longer tail. How can we verify this? The Poisson distribution
# is 1-parameter: if you know its mean, you know everything about it. So
# let's calculate the means of each distribution and overlay the Poisson
# distribution that that implies.
# Recall mean(x) = sum( x*p(x) ). Therefore...
mdn1_mean = np.sum(mdn1_counts * mdn1_prob)
pdr5_mean = np.sum(pdr5_counts * pdr5_prob)
print('mdn1 mean is '+str(mdn1_mean))
print('pdr5 mean is '+str(pdr5_mean))

# The means look plausible. So let's compute the theory curves. Let's do
# a little extra work to be lazy and write a function to compute the
# Poisson distribution.
def poisson_dist(mean, length):
    """
    A function to calculate the Poisson distribution. Note the distribution is
    not normalized in the sense that if you sum over the output_array, the
    result will not be 1; instead it respects the actual values the function
    should take on the finite range of support.
    Recall a Poisson distribution with mean mu has the form
    p(n) = exp(-mu) * mu**n / n!

    Inputs:
    -------
    mean: the mean of the distribution to calculate.
    length: length of support to cover (from 0 to length-1)

    Returns:
    output_array: the array to populate with probabilities
    """
    output_array = np.zeros(length)
    # Loop over copy number
    for m in range(length):
        # If you implement the distribution as above, you will have
        # round-off errors that wreck the computation. One solution is to
        # compute in logs, then apply exponential at end.
        # First need log n! = log1 + log2 + ... + log(n-1) + logn.
        # Note np.arange(1,n+1) returns 1,2,...,n-1,n as we need.
        log_fact = np.sum(np.log(np.arange(1, m+1)))
        output_array[m] = np.exp(-mean + m*np.log(mean) - log_fact)
    return output_array

# Now let's compute the theory distributions for both cases.
mdn1_theory = poisson_dist(mdn1_mean, len(mdn1_counts))
pdr5_theory = poisson_dist(pdr5_mean, len(pdr5_counts))

# Finally we can add the distributions to our data plots. Note how we can
# choose which plot to add stuff to by calling plt.figure(x) first.
plt.figure(1)
plt.plot(mdn1_counts, mdn1_theory, 'g-')
plt.figure(2)
plt.plot(pdr5_counts, pdr5_theory, 'b-')
plt.show()

# We don't need any fancy curve fitting to see the story here! mdn1
# fits quite well to the Poisson distribution, but pdr5 does not fit at all!
# It has a much longer tail, which suggests something else is at play,
# perhaps regulation, but that is a matter for other experiments to settle.
