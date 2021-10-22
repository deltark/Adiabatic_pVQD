import numpy as np
import matplotlib.pyplot as plt
import json

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exactGS = np.load('data/exactGS.npy')
# exactFE = np.load('data/exactFE.npy')

VQD = json.load(open('data/VQD/noisy_shots8k.dat'))
times = VQD["times"]

plt.plot(times, exactGS, linestyle='--', color='black')
plt.plot(times, VQD['E(t)'])


# plt.errorbar(times, old_ansatz['E'])
# plt.errorbar(times, VQD['E'])
# plt.plot(times_shots, shots['E'])
# plt.hlines(exact[0], 0, 1, linestyle = '--', color = 'black')
# plt.legend(['Exact GS', r'T=10.0', r'T=1.0', r'T=0.1'])
# plt.legend(['Exact GS','depth=2','depth=3','depth=4'])
plt.xlabel('t')
plt.ylabel('E(t)')
plt.tight_layout()
# plt.savefig('graphs/T_compar_VQD.eps')
plt.show()

# plt.plot(timesmag[1:], VQDmagd4['final_F'], 'o')
# plt.show()
