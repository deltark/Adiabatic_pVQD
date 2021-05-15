import numpy as np
import matplotlib.pyplot as plt
import json

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exact = np.load('data/J_param_025/exactGS.npy')
# exact = np.load('data/exactGS_J150.npy')
# trotter = json.load(open('data/J_param_025/trotter.dat'))
# VQD = json.load(open('data/J_param_025/VQD.dat'))
VQD = json.load(open('data/sv_J025_T4.dat'))
# shots = json.load(open('data/trial_shots_J025.dat'))

# times = [i*0.01 for i in range(len(VQD['E']))]
# times = [i*0.05 for i in range(len(VQD['E']))]
times = VQD['times']

plt.figure()
# plt.plot(times, trotter['E'])
plt.plot(times, VQD['E'])
# plt.plot(times_shots, shots['E'])
plt.hlines(exact[0], 0, times[len(times)-1], linestyle = '--', color = 'black')
# plt.hlines(exact[1], 0, 2.0, linestyle = '--', color = 'black')
# plt.legend(['Trotter statevector','pVQD statevector', 'GS energy'])
plt.xlabel('t')
plt.ylabel('E')
plt.tight_layout()
# plt.savefig('graphs/J_080.eps')
plt.show()
