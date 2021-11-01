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
	dt         = 0.3
	# n_dt       = 40
	# tmax       = dt*n_dt
	shots      = nothing
end


n_dt       = tmax/dt
print("n_dt: ",n_dt)

t_steps      = [dt*i for i in 0:n_dt]

circ = chain(n_spins, put(i=> H) for i in 1:n_spins) #Hadamard layer
# init_state = zero_state(n_spins) |> circ


trotter_steps = [3, 5, 7, 8, 10, 100]

for nsteps in trotter_steps

	obs_values = []
	
	for (k,t_step) in enumerate(t_steps)

		println("----------------------------------------------------")
		println("Time: $t_step")

		subdt = t_step/nsteps
		substeps = [subdt*i for i in 1:nsteps]

		trott_circ = chain(ising_trotter_step(n_spins,subdt,substep/tmax*J,B) for substep in substeps)
		tot_circ = chain(circ, trott_circ)
		# YaoPlots.plot(trott_circ) |> SVG("trott_circ_nsteps"*string(nsteps)*"_tstep"*string(k)*".svg")
		obs = hamiltonian(n_spins, last(substeps)/tmax*J, B) #time-dependent observable
		state = (zero_state(n_spins) |> tot_circ)
		value = real.(expect(obs, state))
		
		push!(obs_values,value)
	end

	### Now save the data on external file

	save_data = true


	if save_data
		res = Dict("energies"=>obs_values,"spins"=> [n_spins],"dt"=>[dt],"times"=>t_steps)
		j_res = JSON.json(res)
		open("data/trotter/T"*string(tmax)*"_dt"*string(dt)*"_nsteps"*string(nsteps)*"_nshots"*string(shots)*"_fixed.dat","w") do j
			write(j,j_res)
		end
	end
end