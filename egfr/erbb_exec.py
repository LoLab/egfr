"""
Model M1a: Extrinsic apoptosis model with "embedded together" model of MOMP.
"""


from pysb import *
#from egfr import shared
from egfr import chen_modules

Model()
from pysb import *
from pysb.integrate import odesolve
from pylab import linspace, plot, xlabel, ylabel, show

# Declare monomers
chen_modules.rec_monomers()
chen_modules.mapk_monomers()

# Generate the upstream and downstream sections
chen_modules.rec_events()
chen_modules.mapk_events()

# Initial protein concentrations
chen_modules.rec_initial()
chen_modules.mapk_initial()


# Declare observables
Observable('ERKPP', ERK(b=None, st='PP'))

time = linspace(0,2400,100)
print "Simulating..."
x = odesolve(model, time)
# Plot the trajectory of ERKPP
plot(time, x['ERKPP'])
xlabel('Time (seconds)')
ylabel('Amount of ERKPP')
show()

