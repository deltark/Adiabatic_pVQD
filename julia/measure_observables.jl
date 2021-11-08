using Yao,YaoExtensions,PyCall,YaoPlots, JSON, Statistics, PlutoUI, Yao.AD, Compose, LaTeXStrings
using QuAlgorithmZoo: Adam, update!


# include all the useful functions for VpVQD
include("pvqd_functions.jl")
import PyPlot; const plt = PyPlot

# Now import useful subroutines from python
pushfirst!(PyVector(pyimport("sys")."path"), "")



## First of all, load data from pVQD file

## Sys data
n_spins   = 3
NN        = 2
depth     = 1
dt        = 0.05
tmax      = 3.0
opt_steps = 50
J         = -1.0
B         = -1.0

step      = "trotter"
shots     = nothing
q_circ    = alternate_timedep_ansatz(n_spins,depth,NN,zeros(500))


data     = JSON.parse(open("data/p-VQD/nqubits"*string(n_spins)*"_depth"*string(depth)*"_T"*string(tmax)*"_dt"*string(dt)*"_nshots"*string(shots)*"_opt"*string(opt_steps)*"_NN"*string(NN)*"_J"*string(J)*".dat","r"))
exact    = JSON.parse(open("data/exact/nqubits"*string(n_spins)*"_T"*string(tmax)*"_dt0.05_J"*string(J)*".dat","r"))
noise    = JSON.parse(open("data/p-VQD/stataverage_12_noisy_8k_depth1.dat","r"))
noiseNN1 = JSON.parse(open("data/p-VQD/stataverage_12_noisy_8k_depth1_NN1.dat","r"))


# Now for every timestep mesure the obs and save the values
obs_values = []

times  = data["times"]
params = data["parameters"]
# obs    = hamiltonian(n_spins, J, B)

for (t,time) in enumerate(times)
	p = params[t]
	#println("Time: $time, parameters $p \n")
	dispatch!(q_circ,p)
	obs = hamiltonian(n_spins, time/tmax*J, B) #time-dependent observable
	value = real.(expect(obs, zero_state(n_spins) |> q_circ))
	push!(obs_values,value)
end


# println(times)
# println(typeof(obs_values[1]))

# Now plot the data

plt.hlines(last(exact["energies"]), 0, (tmax+0.1), linestyle="--", color="grey")
plt.plot(exact["times"],exact["energies"],linestyle="dashed",color="black")
plt.plot(times,obs_values,marker="o",linestyle="",markersize=4)
plt.errorbar(times,noiseNN1["mean_energy"],yerr=noiseNN1["std"],marker="^",linestyle="",markersize=4,capsize=2)
plt.errorbar(times,noise["mean_energy"],yerr=noise["std"],marker="d",linestyle="",markersize=4,capsize=2)

#Plot to show how much worse Trotter is
# time_steps = [0.4, 0.6]
# legend = ["exact"]
# for dt in time_steps
# 	trotter  = JSON.parse(open("data/trotter/T"*string(tmax)*"_dt"*string(dt)*"_nshots"*string(shots)*".dat","r"))
# 	plt.plot(trotter["times"],trotter["energies"],"o",markersize=4,linestyle="-")
# 	push!(legend, ("dt = "*string(dt)))
# end

# trotter_steps = [7, 8, 10]
# # trotter_steps = [100]
# legend = ["exact evol"]
# for nsteps in trotter_steps
# 	trotter  = JSON.parse(open("data/trotter/T"*string(tmax)*"_dt"*string(dt)*"_nsteps"*string(nsteps)*"_nshots"*string(shots)*"_fixed.dat","r"))
# 	plt.plot(trotter["times"],trotter["energies"],"o",markersize=4,linestyle="-")
# 	push!(legend, ("nsteps = "*string(nsteps)))
# end
# push!(legend, "exact final")

plt.xlabel(L"t")
plt.ylabel(L"$\langle H(t) \rangle$")
# plt.legend(legend)
#
plt.gcf()


############
# data1     = JSON.parse(open("data/p-VQD/nqubits"*string(n_qubits)*"_depth"*string(depth)*"_T"*string(tmax)*"_dt"*string(dt)*"_nshots"*string(shots)*"_opt"*string(opt_steps)*"_NN1_J"*string(J)*".dat","r"))
# data2     = JSON.parse(open("data/p-VQD/nqubits"*string(n_qubits)*"_depth"*string(depth)*"_T"*string(tmax)*"_dt"*string(dt)*"_nshots"*string(shots)*"_opt"*string(opt_steps)*"_NN2_J"*string(J)*".dat","r"))

# params1 = data1["parameters"]
# params2 = data2["parameters"]

# obs = hamiltonian(n_spins, J, B)

# p1 = last(params1)
# #println("Time: $time, parameters $p \n")
# q_circ    = alternate_timedep_ansatz(n_spins,depth,1,zeros(500))
# dispatch!(q_circ,p1)

# value1 = real.(expect(obs, zero_state(n_spins) |> q_circ))

# p2 = last(params2)
# #println("Time: $time, parameters $p \n")
# q_circ    = alternate_timedep_ansatz(n_spins,depth,2,zeros(500))
# dispatch!(q_circ,p2)

# value2 = real.(expect(obs, zero_state(n_spins) |> q_circ))

# print(value1)
# print(value2)

# print(value2-value1)