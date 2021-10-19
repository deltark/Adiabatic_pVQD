using Yao,YaoExtensions,PyCall,YaoPlots, JSON, Statistics, Yao.AD, Compose
using QuAlgorithmZoo: Adam, update!

# include all the useful functions for VpVQD
include("pvqd_functions.jl")
import PyPlot; const plt = PyPlot

function ansatz(n, depth, params)
    count = 1
	circ = chain(n, put(i=> H) for i in 1:n)
	
	for d in 1:depth
		for i in 1:n
			#push!(circ, Rx(params[count]))
			push!(circ,chain(n,put(i=>Ry(params[count]))))
			count = count +1
		end

        for i in 1:n-1
            push!(circ,chain(n,cnot(i,i+1)))
        end

		for i in 1:n
			push!(circ,chain(n,put(i=>Ry(params[count]))))	
			count = count+1
		end
	end
    # print("count= ",count)
	return circ
end

n = 3; depth = 1;
J = -1; B = -1;
# n_params = depth*((n-1) + (n-2) + n)
n_params = depth*2*n
# print("nparams =", n_params)
params = zeros(n_params)
circuit = ansatz(n, depth, params)
# circuit = dispatch!(variational_circuit(n, depth),:random);

h = hamiltonian(n, J, B);
# h = heisenberg(n)

energies = []

opt_steps = 10
eta = 1e-3
for i = 1:opt_steps
    grad = faithful_grad(h, zero_state(n)=>circuit; nshots=100)
    # _, grad = expect'(h, zero_state(n)=>circuit)
    energy = expect(h, zero_state(n)=>circuit)
    push!(energies, energy)
    global params -= eta*grad
    dispatch!(circuit, params)
    # global circuit = alternate_timedep_ansatz(n, depth, params - eta*grad)
    # println("Step $i, energy = $(energy)")
end

# gatecount(circuit)

# close("all")
# exact    = last(JSON.parse(open("data/exact/T2.0_dt0.05.dat","r"))["energies"])
# plt.hlines(exact, 1, opt_steps, linestyle="dashed", color="black")
# plt.plot(1:opt_steps, energies)
# plt.xlabel("step")
# plt.ylabel("energy")
# plt.gcf()