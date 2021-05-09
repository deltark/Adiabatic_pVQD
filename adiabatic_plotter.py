import numpy as np
import matplotlib.pyplot as plt
import json

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exact = np.load('data/J_param_050/exactGS.npy')
trotter = json.load(open('data/J_param_050/trotter.dat'))
VQD = json.load(open('data/J_param_050/VQD.dat'))

times = [i*0.01 for i in range(len(VQD['E']))]

plt.figure()
plt.plot(times, trotter['E'])
plt.plot(times, VQD['E'])
plt.hlines(exact[0], 0, 2.0, linestyle = '--', color = 'black')
# plt.hlines(exact[1], 0, 2.0, linestyle = '--', color = 'black')
plt.legend(['Trotter statevector','pVQD statevector', 'GS energy'])
plt.xlabel('t')
plt.ylabel('E')
plt.tight_layout()
plt.savefig('graphs/J_050.eps')
plt.show()
