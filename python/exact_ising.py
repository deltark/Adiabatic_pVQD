from qutip import *
from qutip.correlation import correlation_2op_1t

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import pickle

import json


PI=np.pi
idmat = qeye(2)


#parameter definition
Lx=3
Ly=1
V= -1
hx= -1
hz = 0.0


pbc = False 

#Hilbert space dimension
d=2


H=0
diss=[]
ii=-1
for i in range(Lx):
	for j in range(Ly):
		ii=ii+1
		op = [idmat for i in range(Lx*Ly)]
		op[ii]=hx*sigmax()
		H=H+tensor(op)
		
		if ii<(Lx-1):
			op = [idmat for i in range(Lx*Ly)]
			op[ii] = V*sigmaz()
			op[ii+1] = sigmaz()
			H=H+tensor(op)
		
            

if pbc:     
	#enforcing periodic BC
	op = [idmat for i in range(Lx*Ly)]
	op[0]=V*sigmaz()
	op[Lx-1]=sigmaz()
	H=H+tensor(op)
            
#observables
obs=[]
Stotx=0
Stoty=0
Stotz=0
ii=-1
for i in range(Lx):
	for j in range(Ly):
		ii=ii+1
		op = [idmat for i in range(Lx*Ly)]
		op[ii]=sigmax()
		Stotx=Stotx+tensor(op)
		op[ii]=sigmay()
		Stoty=Stoty+tensor(op)
		op[ii]=sigmaz()
		Stotz=Stotz+tensor(op)
obs.append(Stotx)
obs.append(Stoty)
obs.append(Stotz)


#observables on single spin

S0x=0
S0y=0
S0z=0
ii=-1
for i in range(1):
	for j in range(Ly):
		ii=ii+1
		op = [idmat for i in range(Lx*Ly)]
		op[ii]=sigmax()
		S0x=S0x+tensor(op)
		op[ii]=sigmay()
		S0y=S0y+tensor(op)
		op[ii]=sigmaz()
		S0z=S0z+tensor(op)
obs.append(S0x)
obs.append(S0y)
obs.append(S0z)

## Correlators
'''
corr12 = [idmat,sigmaz(),sigmaz()]
corr12 = tensor(corr12)
corr02 = [sigmaz(),idmat,sigmaz()]
corr02 = tensor(corr02)
'''
## Auto correlators
'''
ac1   = [idmat,sigmaz()]
ac1   = tensor(ac1)

ac2   = [sigmaz(),idmat]
ac2   = tensor(ac2)
'''

#obs.append(corr12)
#obs.append(corr02)


print("Hamiltonian")
print(H)

s2 = 1/np.sqrt(2)
#initial state uparrow
op = [basis(2, 0) for i in range(Lx)]
psi0up=tensor(op)


    
times = [0.05*i for i in range(400)]
opts  = Options(store_final_state=True,store_states=False)
data  = sesolve(H, psi0up, times, obs, options=opts)

'''
states = data.states
stat_list = []

for i in range(len(states)):
	ts = []
	for j in range(4):
		ts.append(states[i][j][0][0])

	stat_list.append(ts)
'''


res = data.expect

log_data = {}
log_data['times'] = times
#log_data['states'] = stat_list
#print(stat_list)

### AUTOCORRELATORS
'''
corr_data_z1z1_re = correlation_2op_1t(H,psi0up,times,[],ac1,ac1)
corr_data_z1z1_im = correlation_2op_1t(H,psi0up,times,[],ac1,-1j*ac1)
corr_data_z2z2_re = correlation_2op_1t(H,psi0up,times,[],ac2,ac2)
corr_data_z2z2_im = correlation_2op_1t(H,psi0up,times,[],ac2,-1j*ac2)

log_data['auto_Sz1Sz1_re'] = list(np.real(corr_data_z1z1_re)) 
log_data['auto_Sz1Sz1_im'] = list(np.real(corr_data_z1z1_im))
log_data['auto_Sz2Sz2_re'] = list(np.real(corr_data_z2z2_re)) 
log_data['auto_Sz2Sz2_im'] = list(np.real(corr_data_z2z2_im))
'''



log_data['Sx']   = list((1/int(Lx))*res[0])
log_data['Sy']   = list((1/int(Lx))*res[1])
log_data['Sz']   = list((1/int(Lx))*res[2])

log_data['Sx_0']   = list(res[3])
log_data['Sy_0']   = list(res[4])
log_data['Sz_0']   = list(res[5])
#log_data['SzSz'] = list(res[3])
#print(list(res[3]))


#json.dump(log_data, open( 'data/open/V_0.25/OBS_exact_result_2spin_J0.25_B1.dat','w'))
json.dump(log_data, open( 'data/long_exact/OBS_exact_ising_'+str(Lx)+'spin_J'+str(V)+'_B'+str(hx)+'.dat','w'))
#fp =  open( 'STATES_exact_ising_'+str(Lx)+'spin_J'+str(V)+'_B'+str(hx)+'.dat','wb')
#pickle.dump(log_data,fp,0)

#fp =  open( 'STATES_exact_ising_'+str(Lx)+'spin_J'+str(V)+'_B'+str(hx)+'.dat',"rb")
#data = pickle.load(fp)
#print(data)

