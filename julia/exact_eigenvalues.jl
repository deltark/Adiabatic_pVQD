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

function ising_ZZ(nqubits,site1,site2)
    oplist = [Matrix(1I,2,2) for j in 1:nqubits]

    matZ = [1 0;
            0 -1]

    oplist[site1], oplist[site2] = matZ, matZ

    mat = kron(oplist...)

    return mat
end

function ising_Hamiltonian(nqubits,J,B)
    Ham = sum(B*ising_X(nqubits,i) + J*ising_ZZ(nqubits,i,i+1) for i in 1:nqubits-1)
    Ham += B*ising_X(nqubits,nqubits)
    return Ham
end

ground_energy(matrix) = findmin(eigen(matrix).values)


##

## Sys data
n_spins   = 20
dt        = 0.05
tmax      = 1.0
n_dt      = tmax/dt
J         = -1.0
B         = -1.0

t_steps      = [dt*i for i in 0:n_dt]

energies = []
for (k, t_step) in enumerate(t_steps)
    global hamt = ising_Hamiltonian(n_spins, t_step/tmax*J, B)
    global ground = ground_energy(hamt)
    push!(energies, ground[1])
end

crit_eigvec = eigen(hamt).vectors[:,ground[2]]

correlations = []
for x in 1:n_spins-1
    correlator = ising_ZZ(n_spins, 1, 1+x)
    expectation = crit_eigvec'*correlator*crit_eigvec
    push!(correlations, expectation)
end

# print(energies,"\n")
# print(correlations)


res = Dict("spins"=> [n_spins],"dt"=>[dt],"times"=>t_steps,"energies"=>energies,"corr"=>correlations)
j_res = JSON.json(res)
open("data/exact/nqubits"*string(n_spins)*"_T"*string(tmax)*"_dt"*string(dt)*"_J"*string(J)*".dat","w") do j
	write(j,j_res)
end