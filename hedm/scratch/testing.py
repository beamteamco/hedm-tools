# coding: utf-8
# %load testing.py
import numpy as np
import hedm.plot
from matplotlib import pyplot as plt
grains = np.genfromtxt('Grains.csv', skip_header=8, names=True, dtype=None)
xyz = grains[np.array(['X', 'Y', 'Z'])].view(np.float64).reshape((-1, 3))
orient = grains[np.array(['O11', 'O12', 'O13',
	'O21', 'O22', 'O23',
	'O31', 'O32', 'O33'])].view(np.float64).reshape((-1, 3, 3))
x, y, z = xyz.T
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='hedm')
scalars = grains['Confidence'].view(np.float64)
sizes = grains['GrainRadius']
sizes = 80*(sizes - sizes.min())/(sizes.max() - sizes.min()) + 40
_ = ax.grains3D(x, y, z, orient, c=scalars, s=sizes)
