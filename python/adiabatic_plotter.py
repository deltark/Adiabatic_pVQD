import numpy as np
import matplotlib.pyplot as plt
import json

spins = 3
depth = 2
ex_params = np.zeros(depth*((spins-1) + (spins-2) + spins))

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exactGS = np.load('data/exactGS.npy')
# exactFE = np.load('data/exactFE.npy')

svVQD = json.load(open('data/VQD/T3_dt05_sv.dat'))
times = svVQD["times"]

onestep = json.load(open('data/VQD/single_step_test.dat'))

optsteps = onestep["iter_number"][0]
fidelities = onestep["interm_F"][0]
shifts = onestep["shifts"][0]
shifts = [list(x) for x in zip(*shifts)]
grads = onestep["gradients"][0]
grads = [list(x) for x in zip(*grads)]
norm_grad = onestep["norm_grad"][0]
njobs = onestep["njobs"]
print(njobs[0])

plt.plot(range(optsteps), norm_grad)
# for i in range(len(ex_params)):
#     plt.plot(range(optsteps), grads[i])

# noisesims = [json.load(open('data/VQD/noisy_shots8k_run'+str(i)+'.dat'))['E(t)'] for i in range(1,13)]

# mean = np.mean(noisesims, 0)
# std  = np.std(noisesims, 0)

# dico = {}
# dico["mean_energy"] = list(mean)
# dico["std"] = list(std)
# dico["times"] = times

# json.dump(dico, open("data/stataverage_12_noisy_8k.dat", 'w+'))


# plt.plot(times, exactGS, linestyle='--', color='black')
# plt.plot(times, svVQD['E(t)'], 'o')
# plt.errorbar(times, mean, yerr = std, marker = 'o', markersize = mksize, linestyle = '')

# plt.legend(['exact', 'statevector', '8k shots + noise'])


# plt.errorbar(times, old_ansatz['E'])
# plt.errorbar(times, VQD['E'])
# plt.plot(times_shots, shots['E'])
# plt.hlines(exact[0], 0, 1, linestyle = '--', color = 'black')
# plt.legend(['Exact GS', r'T=10.0', r'T=1.0', r'T=0.1'])
# plt.legend(['Exact GS','depth=2','depth=3','depth=4'])
# plt.xlabel('t')
# plt.ylabel('E(t)')
# plt.tight_layout()
# plt.savefig('graphs/T_compar_VQD.eps')
plt.show()
