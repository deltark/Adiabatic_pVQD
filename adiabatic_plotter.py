import numpy as np
import matplotlib.pyplot as plt
import json

plt.rcParams.update({'font.size': 16, "text.usetex": True})
mksize = 4

exact  = json.load(open('data/exact_result_J0.25_B1.dat'))
adiab = json.load(open('trial.dat'))

times = exact['times'][:41]

plt.figure()
plt.plot(times, adiab['E'], 'o', markersize=mksize)
plt.xlabel('t')
plt.ylabel('E')
plt.tight_layout()
plt.show()
