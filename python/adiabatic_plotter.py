import numpy as np
import matplotlib.pyplot as plt
import json

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exactGS = np.load('data/exactGS.npy')
exactFE = np.load('data/exactFE.npy')

VQD1 = json.load(open('data/VQD/T1.dat'))
VQD10 = json.load(open('data/VQD/T10.dat'))
VQDmag = json.load(open('data/VQD/Magnus2T10_Nt50.dat'))
VQDmagd3 = json.load(open('data/VQD/Magnus2T10_Nt50_depth3.dat'))
VQDmagd4 = json.load(open('data/VQD/Magnus2T10_Nt50_depth4.dat'))

VQD10d3 = json.load(open('data/VQD/T10_Nt50_depth3.dat'))
VQD10d4 = json.load(open('data/VQD/T10_Nt50_depth4.dat'))

trottermag = json.load(open('data/trotter/MagnusT1_Nt50.dat'))
trottermag2 = json.load(open('data/trotter/Magnus2T1_Nt50.dat'))

# times = old_ansatz['times']
# times10 = np.array(trotter10['times'])
times10 = np.array(VQD10['times'])/10
timesmag = np.array(VQDmag['times'])/10
times1 = VQD1['times']
times = trottermag['times']
# times02 = trotter02['times']

plt.figure()
plt.plot(times1, exactGS, linestyle='--', color='black')
# plt.plot(times1, exactFE, linestyle='--', color='grey')

# plt.plot(times10, trotter10['E(t)'])
# plt.plot(times1, trotter1['E(t)'])
# plt.plot(times1, trotter01['E(t)'])
# plt.plot(times, trottermag['E(t)'])
# plt.plot(times, trottermag2['E(t)'])

plt.plot(times10, VQD10['E(t)'])
# plt.plot(times1, VQD1['E(t)'])
# plt.plot(timesmag, VQDmag['E(t)'])
# plt.plot(timesmag, VQDmagd3['E(t)'])
# plt.plot(timesmag, VQDmagd4['E(t)'])

plt.plot(timesmag, VQD10d3['E(t)'])
plt.plot(timesmag, VQD10d4['E(t)'])


# plt.errorbar(times, old_ansatz['E'])
# plt.errorbar(times, VQD['E'])
# plt.plot(times_shots, shots['E'])
# plt.hlines(exact[0], 0, 1, linestyle = '--', color = 'black')
# plt.legend(['Exact GS', r'T=10.0', r'T=1.0', r'T=0.1'])
plt.legend(['Exact GS','depth=2','depth=3','depth=4'])
plt.xlabel('t/T')
plt.ylabel('E')
plt.tight_layout()
# plt.savefig('graphs/T_compar_VQD.eps')
plt.show()

# plt.plot(timesmag[1:], VQDmagd4['final_F'], 'o')
# plt.show()