import numpy as np
import matplotlib.pyplot as plt
import json

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exactGS = np.load('data/exactGS.npy')
exactFE = np.load('data/exactFE.npy')
trotter10 = json.load(open('data/trotter/T10.dat'))
trotter1 = json.load(open('data/trotter/T1.dat'))
trotter01 = json.load(open('data/trotter/T01.dat'))
# VQD = json.load(open('data/ansatz_compar/new.dat'))

# times = old_ansatz['times']
times10 = trotter10['times']
times1 = trotter1['times']
# times02 = trotter02['times']

plt.figure()
plt.plot(times1, exactGS, linestyle='--', color='black')
# plt.plot(times1, exactFE, linestyle='--', color='grey')
plt.plot(times1, trotter01['E(t)'])
plt.plot(times1, trotter1['E(t)'])
plt.plot(times10, trotter10['E(t)'])
# plt.errorbar(times, old_ansatz['E'])
# plt.errorbar(times, VQD['E'])
# plt.plot(times_shots, shots['E'])
# plt.hlines(exact[0], 0, 1, linestyle = '--', color = 'black')
plt.xlabel('t/T')
plt.ylabel('E')
# plt.legend(['GS energy','Local Ansatz','Long-range Ansatz'])
plt.tight_layout()
# plt.savefig('graphs/ansatz_compar.eps')
plt.show()
