using Yao,YaoExtensions,PyCall,YaoPlots, JSON, Statistics, Yao.AD, Compose
# using QuAlgorithmZoo: Adam, update!

# include all the useful functions for VpVQD
include("pvqd_functions.jl")
import PyPlot; const plt = PyPlot

# Now import useful subroutines from python
pushfirst!(PyVector(pyimport("sys")."path"), "")

begin
	n_spins   = 3
	J          = -1.0
	B          = -1.0
	tmax       = 3.0
	# n_dt       = 40
	# tmax       = dt*n_dt
	shots      = nothing
end

time_steps = [0.4, 0.6]

for dt in time_steps

	n_dt       = tmax/dt
	print("n_dt: ",n_dt)

	t_steps      = [dt*i for i in 0:n_dt]

	obs_values = []

	circ = chain(n_spins, put(i=> H) for i in 1:n_spins) #Hadamard layer
	state = zero_state(n_spins) |> circ


	for (k,t_step) in enumerate(t_steps)

		println("----------------------------------------------------")
		println("Time: $t_step")

		trott_circ = ising_trotter_step(n_spins,dt,t_step/tmax*J,B)
		obs = hamiltonian(n_spins, t_step/tmax*J, B) #time-dependent observable
		state = (state |> trott_circ)
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
end