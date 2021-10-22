using Yao,YaoExtensions,PyCall,YaoPlots, JSON, Statistics, Yao.AD, Compose
using QuAlgorithmZoo: Adam, update!

# include all the useful functions for VpVQD
include("pvqd_functions.jl")
import PyPlot; const plt = PyPlot

# Now import useful subroutines from python
pushfirst!(PyVector(pyimport("sys")."path"), "")

begin
	n_spins   = 3
	J          = -1.0
	B          = -1.0
	dt         = 0.05
	tmax       = 3.0
	n_dt       = tmax/dt
	print("n_dt: ",n_dt)
	# n_dt       = 40
	# tmax       = dt*n_dt
	shots      = nothing
end

t_steps      = [dt*i for i in 0:n_dt]

obs_values = []

state = zero_state(n_spins)

for (k,t_step) in enumerate(t_steps)

    println("----------------------------------------------------")
	println("Time: $t_step")

    trott_circ = ising_trotter_step(n_spins,dt,t_step/tmax*J,B)
    obs = hamiltonian(n_spins, t_step/tmax*J, B) #time-dependent observable
    global state = (state |> trott_circ)
	value = real.(expect(obs, state))
     
	push!(obs_values,value)
end

### Now save the data on external file

save_data = true


if save_data
	res = Dict("energies"=>obs_values,"spins"=> [n_spins],"dt"=>[dt],"times"=>t_steps)
	j_res = JSON.json(res)
	open("data/trotter/T"*string(tmax)*"_dt"*string(dt)*"_nshots"*string(shots)*".dat","w") do j
		write(j,j_res)
	end
end