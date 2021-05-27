import numpy as np
import matplotlib.pyplot as plt
import json

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exact = np.load('data/J_param_080/exactGS.npy')
# trotter = json.load(open('data/trotter_dt005.dat'))
VQD = json.load(open('data/VQD_J080_dt005_ths5.dat'))

times = VQD['times']
# times = trotter['times']

plt.figure()
# plt.plot(times, trotter['E'])
plt.plot(times, VQD['E'])
# plt.plot(times_shots, shots['E'])
plt.hlines(exact[0], 0, times[len(times)-1], linestyle = '--', color = 'black')
plt.xlabel('t')
plt.ylabel('E')
plt.tight_layout()
# plt.savefig('graphs/J_150.eps')
plt.show()
