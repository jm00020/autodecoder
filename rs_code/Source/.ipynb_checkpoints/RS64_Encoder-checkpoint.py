
# coding: utf-8

# In[1]:


import numpy as np
from Source.GF64 import GFE
from Source import GF_FUNC as GF_F


# In[2]:


def Generation_Polynomial(Z_RS, order_alpha):
    a_x=[0]
    a_x.append(Z_RS[0])
    b_x=[0]
    for i in range(1,len(Z_RS)):
        b_x.append(Z_RS[i])
#     print('a(x)=',a_x)
#     print('b(x)=',b_x)
        c_x=GF_F.POL_MUL(a_x,b_x, order_alpha)  ## c(x)=a(x)b(x)
        a_x=c_x
        b_x=[0]

#     print('c(x):',c_x)
    g_x=c_x  
    g_x.reverse()   ## g(x)=g_0+g_1*x+...+g_n-k*x^n-1
#     print('g(x)=g_0+g_1*x+...+g_n-k*x^n-1: ',g_x)
#     print('length of g:',len(g_x))
    
    return g_x


# In[12]:


def RS_encoding(msg_RS, g_x, order_alpha, codeword_size):
    parity_symbols=GF_F.Parity_RS(msg_RS, g_x, order_alpha, codeword_size)
    codeword_RS=np.hstack((parity_symbols, msg_RS))
#     print('msg:',msg_RS)
#     print('parity symbols=',parity_symbols)
#     print('codeword:',codeword_RS)
#     print('length of codeword:',len(codeword_RS))
    
    return codeword_RS


# In[14]:


def Encoder(order_alpha, k_RS,codeword_size, g_x):
    msg_RS = np.random.randint(-1,order_alpha, size = (codeword_size,k_RS))
    
    codeword_RS = RS_encoding(msg_RS, g_x, order_alpha, codeword_size)
    
    return msg_RS, codeword_RS


def trans_vector(codeword_RS, codeword_size, n_RS, m):
    vector_codeword_RS = GFE[codeword_RS.astype(int)+1].reshape(codeword_size, n_RS * m)
    
    
    return vector_codeword_RS
