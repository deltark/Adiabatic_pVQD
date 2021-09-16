import numpy as np
import matplotlib.pyplot as plt
import json

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exact = np.load('data/J_param_080/exactGS.npy')
sv = json.load(open('data/trotter/sv_T2_J080.dat'))
trotter = json.load(open('data/trotter/T2_J080.dat'))

times = trotter['times']

plt.figure()
plt.errorbar(times, sv['E'])
plt.errorbar(times, trotter['E'], yerr=trotter['err_E'], marker = 'o', markersize = mksize, linestyle='')
plt.hlines(exact[0], 0, times[len(times)-1], linestyle = '--', color = 'black')
# plt.hlines(exact[1], 0, times[len(times)-1], linestyle = '--', color = 'grey')
plt.legend(['GS energy', 'Statevector', '10000 shots'])
plt.xlabel('t')
plt.ylabel('E')
plt.tight_layout()
plt.savefig('graphs/final_trotter_T2_J080.eps')
plt.show()
