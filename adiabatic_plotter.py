import numpy as np
import matplotlib.pyplot as plt
import json

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exact = np.load('data/J_param_080/exactGS.npy')
# trotter = json.load(open('data/trotter_dt005.dat'))
VQD = json.load(open('data/VQD_J080_dt005_ths6.dat'))
# VQD = json.load(open('data/trial_J025_oldansatz.dat'))
shots5 = json.load(open('data/cluster/VQD_shots5_J080.dat'))
shots6 = json.load(open('data/VQD_shots6_J080.dat')) ###KEEP

times = VQD['times']
# times = trotter['times']

plt.figure()
# plt.plot(times, trotter['E'])
plt.hlines(exact[0], 0, times[len(times)-1], linestyle = '--', color = 'black')
plt.plot(times, VQD['E'], marker='o', markersize=mksize, linestyle='')
# plt.errorbar(times, shots5['E'], yerr=shots5['err_E'], marker='o', linestyle='', markersize=mksize, label='shots5')
# plt.errorbar(times, shots6['E'], yerr=shots6['err_E'], marker='o', linestyle='', markersize=mksize, label='shots6')
plt.plot(times, shots5['E'], marker='o', markersize=mksize, linestyle='')
plt.plot(times, shots6['E'], marker='o', markersize=mksize, linestyle='')
plt.xlabel('t')
plt.ylabel('E')
plt.legend(['Statevector', r'$10^5$ shots', r'$10^6$ shots', 'GS energy'])
plt.tight_layout()
plt.savefig('graphs/J_080_shots.eps')
plt.show()
