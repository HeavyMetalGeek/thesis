import matplotlib.pyplot as plt
from meanfields import MeanFields

# Process time average of mean values for 100% run
mfA = MeanFields(meanPath='./postProcessing/averaging/all/')
mfA.time_mean(min_time=275.5)

# Process time average of mean values for 60% run
mfB = MeanFields(meanPath='./postProcessing_60/averaging/all/')
mfB.time_mean(min_time=20)

# Calculate maximum mean velocity errors
mfA_max_err = (mfA.U['log-law']-mfA.U['simulation']).abs().max()
print("100% simulation max error = {}".format(mfA_max_err))
mfB_max_err = (mfB.U['log-law']-mfB.U['simulation']).abs().max()
print(" 60% simulation max error = {}".format(mfB_max_err))

# Plots
figs = []
axes = []
fig, ax = plt.subplots(1)
xA = mfA.U.values
yA = mfA.U.index
xB = mfB.U.values
yB = mfB.U.index

ax.plot(xA[:,0], 
	    yA, 
	    marker='s', 
	    color='red', 
	    ls='', 
	    label='100% simulation', 
	    markeredgecolor='none')

ax.plot(xB[:,0], 
	    yB, 
	    marker='o', 
	    color='blue', 
	    ls='', 
	    label=' 60% simulation', 
	    markeredgecolor='none')

ax.plot(xA[:,1], 
	    yA, 
	    marker='', 
	    color='black', 
	    ls='--', 
	    label='log-law', 
	    markeredgecolor='none')



ax.set_xlabel('Mean Streamwise Velocity')
ax.set_ylabel('z/H')
ax.legend(loc='best')
fig.tight_layout()
ax.set_yscale('log')
ax.set_xlim([4.875, 12])
figs.append(fig)
axes.append(ax)

fig, ax = plt.subplots(1)
xA = mfA.phi.values
yA = mfA.phi.index
xB = mfB.phi.values
yB = mfB.phi.index

ax.plot(xA[:,0], 
	    yA, 
	    marker='s', 
	    color='red', 
	    ls='', 
	    label='100% simulation', 
	    markeredgecolor='none')

ax.plot(xB[:,0], 
	    yB, 
	    marker='o', 
	    color='blue', 
	    ls='', 
	    label=' 60% simulation', 
	    markeredgecolor='none')

ax.plot(xA[:,1], 
	    yA, 
	    marker='', 
	    color='black', 
	    ls='--', 
	    label='log-law', 
	    markeredgecolor='none')



ax.set_xlabel('$\Phi$')
ax.set_ylabel('z/H')
ax.legend(loc='best')
fig.tight_layout()
ax.set_xlim([0, 2.5])

mfA.plot_val('velocity_variance')
mfB.plot_val('velocity_variance')

plt.show()