
import numpy as np
m=5
g=(1,0,1,0,0,1)  
zero_vector=np.zeros(m)
sr_init=np.zeros(m)
sr=np.zeros(m)
sr_init[0]=1
for i in range(1,m): sr_init[i]=0
GFE={}
GFE[-1]=list(zero_vector) # array value in Dictionary make very strange result!!! Thus list type is used!!!
GFE[0]=list(sr_init)
for i in range(1,2**m-1):
    feedback=sr_init[m-1]
    sr[0]=feedback
    for j in range(1,m): sr[j]=(sr_init[j-1]+feedback*g[j])%2
    for j in range(m):   sr_init[j]=sr[j]
    GFE[i]=list(sr_init)
for i in range(-1,2**m-1):    print('GF[',i,']=', GFE[i])