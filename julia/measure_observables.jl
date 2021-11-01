using Yao,YaoExtensions,PyCall,YaoPlots, JSON, Statistics, PlutoUI, Yao.AD, Compose, LaTeXStrings
using QuAlgorithmZoo: Adam, update!


# include all the useful functions for VpVQD
include("pvqd_functions.jl")
import PyPlot; const plt = PyPlot

# Now import useful subroutines from python
pushfirst!(PyVector(pyimport("sys")."path"), "")



## First of all, load data from pVQD file

## Sys data
n_spins   = 7
NN        = n_spins-1
depth     = 3
dt        = 0.05
tmax      = 10.0
opt_steps = 50
J         = -1.0
B         = -1.0

step      = "trotter"
shots     = nothing
q_circ    = alternate_timedep_ansatz(n_spins,depth,zeros(500))


# data     = JSON.parse(open("data/p-VQD/depth"*string(depth)*"_T"*string(tmax)*"_dt"*string(dt)*"_nshots"*string(shots)*"_opt"*string(opt_steps)*".dat","r"))
exact    = JSON.parse(open("data/exact/nqubits"*string(n_spins)*"_T"*string(tmax)*"_dt0.05.dat","r"))
# noise    = JSON.parse(open("data/p-VQD/stataverage_12_noisy_8k_depth1.dat","r"))


# data     = JSON.parse(open("data/p-VQD/nqubits"*string(n_spins)*"_depth"*string(depth)*"_T"*string(tmax)*"_dt"*string(dt)*"_nshots"*string(shots)*"_opt"*string(opt_steps)*"_NN"*string(NN)*".dat","r"))
# times  = data["times"]
# params = data["parameters"]

# Now for every timestep mesure the obs and save the values
obs_values = []

# times  = data["times"]
# params = data["parameters"]
# obs    = spin_z(n_spins, 3)
obs    = hamiltonian(n_spins, J, B)
# print(obs)

# for (t,time) in enumerate(times)
# 	p = params[t]
# 	#println("Time: $time, parameters $p \n")
# 	dispatch!(q_circ,p)
# 	obs = hamiltonian(n_spins, time/tmax*J, B) #time-dependent observable
# 	value = real.(expect(obs, zero_state(n_spins) |> q_circ))
# 	push!(obs_values,value)
# end

# Measure correlation function as a function of distance

legend = []
markers = ["o","^","s","d","P"]
for (k,tmax) in enumerate([1.0, 2.0, 4.0, 7.0, 10.0])
	if tmax==7.0
		global dt = 0.07
	elseif tmax==10.0
		global dt = 0.1
	end
	data     = JSON.parse(open("data/p-VQD/nqubits"*string(n_spins)*"_depth"*string(depth)*"_T"*string(tmax)*"_dt"*string(dt)*"_nshots"*string(shots)*"_opt"*string(opt_steps)*"_NN"*string(NN)*".dat","r"))
	times  = data["times"]
	params = data["parameters"]

	t = length(times)
	p = last(params)
	dispatch!(q_circ,p)

	obs_corr = []
	xrange = 1:NN
	for x in xrange
		site = 1
		obs = spin_zz(n_spins, site, site+x)
		value = real.(expect(obs, zero_state(n_spins) |> q_circ))
		push!(obs_corr, value)
	end
	push!(legend, "tmax="*string(tmax))
	plt.loglog(xrange, obs_corr, marker=markers[k], markersize=4)

	# plt.plot(xrange, obs_corr, "o")
end

plt.xlabel(L"x")
plt.ylabel("ZZ Correlation")
plt.legend(legend)
plt.gcf()

# println(times)
# println(typeof(obs_values[1]))

# Now plot the data

# plt.hlines(last(exact["energies"]), 0, (tmax+0.1), linestyle="--", color="grey")
# plt.plot(exact["times"],exact["energies"],linestyle="dashed",color="black")
# plt.plot(times,obs_values,marker="o",linestyle="",markersize=4)
# plt.errorbar(times,noise["mean_energy"],yerr=noise["std"],marker="o",linestyle="",markersize=4)


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

# plt.xlabel(L"t")
# plt.ylabel(L"$\langle H(t) \rangle$")
# plt.legend(legend)
# #
# plt.gcf()