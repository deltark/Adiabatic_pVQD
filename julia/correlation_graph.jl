using Yao,YaoExtensions,PyCall,YaoPlots, JSON, Statistics, PlutoUI, Yao.AD, Compose, LaTeXStrings
using QuAlgorithmZoo: Adam, update!


# include all the useful functions for VpVQD
include("pvqd_functions.jl")
import PyPlot; const plt = PyPlot

# Now import useful subroutines from python
pushfirst!(PyVector(pyimport("sys")."path"), "")

## Sys data
n_spins   = 7
NN        = n_spins-1
depth     = 3
dt        = 0.05
tmax      = 1.0
opt_steps = 50
J         = -1.0
B         = -1.0

step      = "trotter"
shots     = nothing
q_circ    = alternate_timedep_ansatz(n_spins,depth,NN,zeros(500))

xrange = 1:NN

# exact = JSON.parse(open("data/exact/nqubits"*string(n_spins)*"_T1.0_dt"*string(dt)*"_J"*string(J)*".dat","r"))
exact    = JSON.parse(open("data/exact/nqubits"*string(n_spins)*"_T"*string(tmax)*"_dt0.05_J"*string(J)*".dat","r"))
plt.plot(xrange,exact["corr"],linestyle="--",color="black",marker="x")

legend = ["exact"]
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
plt.ylabel(L"C(T,t=T,x)")
plt.legend(legend)
plt.gcf()