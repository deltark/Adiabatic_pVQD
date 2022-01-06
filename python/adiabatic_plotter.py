import numpy as np
import matplotlib.pyplot as plt
import json

spins = 3
depth = 1
ex_params = np.zeros(depth*spins + depth*(spins-1)) #NN1
# ex_params = np.zeros(depth*((spins-1) + (spins-2) + spins)) #NN2

plt.rcParams.update({'font.size': 18, "text.usetex": True})
mksize = 4

exactGS = np.load('data/exactGS.npy')
# exactFE = np.load('data/exactFE.npy')

svVQD = json.load(open('data/VQD/hybrid_test_remote.dat'))
times = svVQD["times"]

# onestep = json.load(open('data/VQD/runtime_ibm_lagos_NN1_iter30_shots2000.dat'))
# optsteps = onestep["iter_number"][0]
# fidelities = onestep["interm_F"][0]
# loss = [1-f for f in fidelities]
# shifts = onestep["shifts"][0]
# shifts = [list(x) for x in zip(*shifts)]
# grads = onestep["gradients"][0]
# grads = [list(x) for x in zip(*grads)]
# err_grad = onestep["err_grad"][0]
# err_grad = [list(x) for x in zip(*err_grad)]
# norm_grad = onestep["norm_grad"][0]
# njobs = onestep["njobs"]
# print(njobs[0])
# print(grads)
# print(err_grad)

# plt.plot(range(optsteps), norm_grad)
# for i in range(len(ex_params)):
#     plt.errorbar(range(optsteps), grads[i], yerr=err_grad[i], capsize=2)
# plt.xlabel("opt step")
# plt.ylabel("gradient")
# plt.tight_layout()
# plt.show()

# plt.semilogy(range(optsteps), loss, 'o')
# plt.xlabel("opt step")
# plt.ylabel("loss")
# # plt.grid(True)
# plt.tight_layout()
# plt.show()

# noisesims = [json.load(open('data/VQD/noisy_shots8k_depth1_NN1_run'+str(i)+'.dat'))['E(t)'] for i in range(1,13)]

# mean = np.mean(noisesims, 0)
# std  = np.std(noisesims, 0)

# dico = {}
# dico["mean_energy"] = list(mean)
# dico["std"] = list(std)
# dico["times"] = times

# json.dump(dico, open("data/stataverage_12_noisy_8k_depth1_NN1.dat", 'w+'))


plt.plot(times, exactGS, linestyle='--', color='black')
plt.plot(times, svVQD['E(t)'], 'o')
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