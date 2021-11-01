using LinearAlgebra, JSON

## Functions

function ising_X(nqubits,site)
    oplist = [Matrix(1I,2,2) for j in 1:nqubits]

    matX = [0 1;
            1 0]

    oplist[site] = matX

    mat = kron(oplist...)

    return mat
end

function ising_ZZ(nqubits,site)
    oplist = [Matrix(1I,2,2) for j in 1:nqubits]

    matZ = [1 0;
            0 -1]

    oplist[site], oplist[site+1] = matZ, matZ

    mat = kron(oplist...)

    return mat
end

function ising_Hamiltonian(nqubits,J,B)
    Ham = sum(B*ising_X(nqubits,i) + J*ising_ZZ(nqubits,i) for i in 1:nqubits-1)
    Ham += B*ising_X(nqubits,nqubits)
    return Ham
end

ground_energy(matrix) = minimum(eigvals(matrix))

##

## Sys data
n_spins   = 7
dt        = 0.05
tmax      = 4.0
n_dt      = tmax/dt
J         = -1.0
B         = -1.0

t_steps      = [dt*i for i in 0:n_dt]

energies = []

for (k, t_step) in enumerate(t_steps)
    hamt = ising_Hamiltonian(n_spins, t_step/tmax*J, B)
    push!(energies, ground_energy(hamt))
end

res = Dict("spins"=> [n_spins],"dt"=>[dt],"times"=>t_steps,"energies"=>energies)
j_res = JSON.json(res)
open("data/exact/nqubits"*string(n_spins)*"_T"*string(tmax)*"_dt"*string(dt)*".dat","w") do j
	write(j,j_res)
end