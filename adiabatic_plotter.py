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

VQD10 = json.load(open('data/VQD/T10.dat'))
VQD1 = json.load(open('data/VQD/T1.dat'))
VQD01 = json.load(open('data/VQD/T01.dat'))

# times = old_ansatz['times']
# times10 = np.array(trotter10['times'])
times10 = np.array(VQD10['times'])/10
times1 = trotter1['times']
# times02 = trotter02['times']

plt.figure()
plt.plot(times1, exactGS, linestyle='--', color='black')
plt.plot(times1, exactFE, linestyle='--', color='grey')

# plt.plot(times10, trotter10['E(t)'])
# plt.plot(times1, trotter1['E(t)'])
# plt.plot(times1, trotter01['E(t)'])

plt.plot(times10, VQD10['E(t)'])
plt.plot(times1, VQD1['E(t)'])
plt.plot(times1, VQD01['E(t)'])


# plt.errorbar(times, old_ansatz['E'])
# plt.errorbar(times, VQD['E'])
# plt.plot(times_shots, shots['E'])
# plt.hlines(exact[0], 0, 1, linestyle = '--', color = 'black')
plt.legend(['Exact GS', r'T=10.0', r'T=1.0', r'T=0.1'])
plt.xlabel('t/T')
plt.ylabel('E')
plt.tight_layout()
# plt.savefig('graphs/T_compar_VQD.eps')
plt.show()
