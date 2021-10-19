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
depth     = 2
dt        = 0.05
tmax      = 2.0
J         = -1.0
B         = -1.0

step      = "magnus"
shots     = 2000
q_circ    = alternate_timedep_ansatz(n_spins,depth,zeros(500))


data     = JSON.parse(open("data/p-VQD/depth"*string(depth)*"_T"*string(tmax)*"_dt"*string(dt)*"_nshots"*string(shots)*"_opt"*string(opt_steps)*"_"*step*".dat","r"))
exact    = JSON.parse(open("data/exact/T"*string(tmax)*"_dt0.05.dat","r"))

# Now for every timestep mesure the obs and save the values
obs_values = []

times  = data["times"]
params = data["parameters"]
# obs    = spin_z(n_spins, 3)
# obs    = hamiltonian(n_spins, J, B)
# print(obs)

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

# plt.plot(exact["times"],exact["energies"],linestyle="dashed",color="black")
plt.plot(times,obs_values)
plt.xlabel(L"t")
plt.ylabel(L"$\langle H(t) \rangle$")
#
plt.gcf()