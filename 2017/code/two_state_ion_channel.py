# the usual imports
import numpy as np
import matplotlib.pyplot as plt

# This short script demonstrates how the open probability for an ion
# channel depends on the energy difference between the open and
# closed conformations.

# First define range to plot over.
delta_E = np.linspace(-5, 5, 100)
# The open probability as derived in class
p_open = 1/(1+np.exp(-delta_E))

# Generate the plot
plt.plot(delta_E, p_open)
# Some fancy formatting for matplotlib: you can include math using
# standard LaTeX syntax if you precede the string with an r.
plt.xlabel(r'$\Delta_E$ ($k_BT$ units)')
plt.ylabel(r'$p_{open}$')
plt.show()
