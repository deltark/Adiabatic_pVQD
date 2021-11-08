## pVQD meets Julia! 

# This file will contain all the useful functions for pVQD and variational circuit in general

# using LinearAlgebra #for testing

###################################################################################################
## GATES

function Rzz(n,i,j,theta)
	circ = chain(n,[cnot(i,j),put(j=>Rz(theta)),cnot(i,j)])
	return circ
end

function Rxx(n,i,j,theta)
	circ = chain(n,[put(i=>Ry(pi/2)),cnot(i,j),put(j=>Rz(theta)),cnot(i,j),put(i=>Ry(-pi/2))])
	return circ
end


function Ryy(n,i,j,theta)
	circ = chain(n,[put(i=>Rx(pi/2)),put(j=>Rx(pi/2)),cnot(i,j),put(j=>Rz(theta)),cnot(i,j),put(i=>Rx(-pi/2)),put(j=>Rx(-pi/2))])
	return circ
end

function Rzy(n,i,j,theta)
	circ = chain(n,[put(j=>Rx(pi/2)),cnot(i,j),put(j=>Rz(theta)),cnot(i,j),put(j=>Rx(-pi/2))])
	return circ
end

###################################################################################################
# ANSATZE

function trotter_ansatz(n,depth,params)
	count = 1
	circ = chain(n, put(i=> I2) for i in 1:n)
	
	

	for d in 1:depth
		# Rx layer
		for i in 1:n
			#push!(circ, Rx(params[count]))
			push!(circ,chain(n,put(i=>Rx(params[count]))))
			count = count +1
		end

		# Rzz layer
		for j in 1:n-1
			push!(circ,Rzz(n,j,j+1,params[count]))	
			count = count+1
		end
	end
	return circ
end

function alternate_trotter_ansatz(n,depth,params)
	count = 1
	circ = chain(n, put(i=> I2) for i in 1:n)
	

	for d in 1:depth
		# Rx/Ry layer

		for i in 1:n
			if mod(d,2) == 1
				#push!(circ, Rx(params[count]))
				push!(circ,chain(n,put(i=>Rx(params[count]))))
				count = count +1
			else 
				#push!(circ, Ry(params[count]))
				push!(circ,chain(n,put(i=>Ry(params[count]))))
				count = count +1
			end	
		end

		# Rzz layer
		for j in 1:n-1
			push!(circ,Rzz(n,j,j+1,params[count]))	
			count = count+1
		end
	end
	return circ
end

function alternate_timedep_ansatz(n,depth,NN,params)
	neigh = NN #nth nearest neighbor for the Rzz layer

	count = 1
	# circ = chain(n, put(i=> I2) for i in 1:n)
	circ = chain(n, put(i=> H) for i in 1:n) #Hadamard layer
	

	for d in 1:depth
		# Rzz layer
		for j in 1:n-1, k in j+1:j+neigh
			if k <= n
				push!(circ,Rzz(n,j,k,params[count]))	
				count = count+1
			end
		end

		# Rx/Ry layer

		for i in 1:n
			if mod(d,2) == 1
				#push!(circ, Rx(params[count]))
				push!(circ,chain(n,put(i=>Rx(params[count]))))
				count = count +1
			else 
				#push!(circ, Ry(params[count]))
				push!(circ,chain(n,put(i=>Ry(params[count]))))
				count = count +1
			end	
		end
	end
	return circ
end


###################################################################################################
# PROJECTORS AND OBSERVABLES

function projector_zero(n)
	prj    = kron([0.5*(I2+Z) for i in 1:n]...)
	#id 	   = kron([I2 for i in 1:n]...)
	return  prj
end	

function projector_site(n,i)
	op_list = [0.5*(I2+Z) for i in 1:n]
	op_list[i] = 0.5*(I2-Z)

	prj = kron(op_list...)

	return prj
end

function spin_x(n,site)
	op_list = []

	for i in 1:n
		if i==site 
			push!(op_list,X)
		else 
			push!(op_list,I2)
		end
	end

	obs = kron(op_list...)

	return obs
end	

function spin_y(n,site)
	op_list = []

	for i in 1:n
		if i==site 
			push!(op_list,Y)
		else 
			push!(op_list,I2)
		end
	end 

	obs = kron(op_list...)

	return obs
end	

function spin_z(n,site)
	op_list = []

	for i in 1:n
		if i==site 
			push!(op_list,Z)
		else 
			push!(op_list,I2)
		end
	end
	
	obs = kron(op_list...)

	return obs
end

function spin_zz(n,site1,site2)
	op_list = []

	for i in 1:n
		if i==site1 || i==site2
			push!(op_list,Z)
		else 
			push!(op_list,I2)
		end
	end
	
	obs = kron(op_list...)

	return obs
end

function hamiltonian(n,J,B)
	op_list = []

	for i in 1:n-1
		push!(op_list,J*spin_zz(n,i,i+1))
	end
	for i in 1:n
		push!(op_list,B*spin_x(n,i))
	end
	
	obs = sum(op_list)

	return obs
end

###################################################################################################
# TROTTER/MAGNUS

function ising_trotter_step(n,dt,J,B)

	circ = chain(n, put(i=> I2) for i in 1:n)
	# Rx layer
	for i in 1:n
		#push!(circ, Rx(B*dt))
		push!(circ,chain(n,put(i=>Rx(2*B*dt))))
		
	end

	# Rzz layer
	for j in 1:n-1
		push!(circ,Rzz(n,j,j+1,2*J*dt))	
	end

	return circ
end

function magnus2_step(n,dt,J,B,tmax)

	circ = chain(n, put(i=> I2) for i in 1:n)

	for i in 1:n-1
		push!(circ,Rzy(n,i,i+1,2*J*B*dt^3/(6*tmax)))
		push!(circ,Rzy(n,i+1,i,2*J*B*dt^3/(6*tmax)))
	end

	return circ
end

#testing

# Zmat = [1 0;
# 	    0 -1]
# Ymat = [0 -1im;
# 	    1im 0]

# RZYalgebra = exp(-1im*kron([Matrix(I,2,2),Zmat,Ymat]...))
# zerostatealg = reshape([1 0 0 0 0 0 0 0], (8,1))

# print("comparison","\n")
# print(RZYalgebra*zerostatealg,"\n")

# zerostate = zero_state(3)
# stateq = zerostate |> Ryz(3,1,2,2)
# print(state(stateq),"\n")
# print((RZYalgebra*zerostatealg)-state(stateq))
