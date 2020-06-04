
# coding: utf-8

# In[1]:


import numpy as np
from GF64 import GFE
import GF_FUNC as GF_F


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


def RS_encoding(msg_RS, g_x, order_alpha):
    parity_symbols=GF_F.Parity_RS(msg_RS, g_x, order_alpha)
    codeword_RS=list(parity_symbols)+list(msg_RS)
#     print('msg:',msg_RS)
#     print('parity symbols=',parity_symbols)
#     print('codeword:',codeword_RS)
#     print('length of codeword:',len(codeword_RS))
    
    return codeword_RS


# In[14]:


def Encoder(m, t_RS, b_RS, order_alpha, k_RS):

    Z_RS=[]
    for i in range(2*t_RS):    
        Z_RS.append((b_RS+i)%order_alpha)
        
    g_x = Generation_Polynomial(Z_RS, order_alpha)
    
    msg_RS=[]
    for i in range(k_RS):
        msg_RS.append(np.random.randint(-1,order_alpha))

#     np.random.seed(700)
    codeword_RS = RS_encoding(msg_RS, g_x, order_alpha)
    
    return codeword_RS, msg_RS


# In[19]:


def Add_Noise(Eb_No, codeword_RS, k_RS, m):
    r_vector_RS=[]
    req_Eb_No = 10**(Eb_No/10)
    Es_No = req_Eb_No * (k_RS/len(codeword_RS))
    Sigma=np.sqrt(1/(2*Es_No))
    for i in range(len(codeword_RS)):
        m_bits=GFE.get(codeword_RS[i])
        tx=[]
        rx=[]
        rx_Hard=[]
        rn=np.random.normal(loc=0.0,scale=Sigma,size=m)
        for j in range(m): 
            tx.append(np.power(-1,m_bits[j]))
        for j in range(m):
            Bi_AWGN=tx[j]+rn[j]
            rx.append(Bi_AWGN)
            if Bi_AWGN >=0 :
                rx_Hard.append(0)
            else:
                rx_Hard.append(1)
#    print('m bits:',m_bits,'Index:',codeword_RS[i])
#    print('tx m bits:',tx,'Index:',codeword_RS[i])
#    print('rx m bits:',rx)
#    print('rx Hard m bits:',rx_Hard)
        r_vector_RS.append(GF_F.GF_Index(rx_Hard))

#     print('rx symbol:',r_vector_RS)
    
    return r_vector_RS

