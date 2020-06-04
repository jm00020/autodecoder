
# coding: utf-8

# In[2]:


import numpy as np
from GF64 import GFE

# In[3]:


def GF_Index(cc):
    for GFindex, GFvector in GFE.items():
        if GFvector == list(cc): 
            return GFindex
    print("Invalid vector of m bits")


# In[4]:


def GF_Add(a,b): # alpha^c = alpha^a + alpha^b
    aa=GFE.get(a)
    bb=GFE.get(b)
    cc=np.add(aa,bb)%2
    return GF_Index(cc)


# In[5]:


def GF_MUL(a,b, order_alpha): 
    if a >= 0 and b>=0: # if a and b are greater than -1, alpha^a * alpha^b =alpha^(c), c=a+b (mod n)
        c=(a+b)%order_alpha
        return c
    else:
        return -1  # if a and/or b is minus, the result is zero vector that is, alpha^-1


# In[6]:


def GF_MUL_Inv(a, order_alpha): # alpha^a * alpha^b = alpha^0
    b=-1*a%order_alpha    
    return b


# In[7]:


def POL_MUL(a_coeff,b_coeff, order_alpha):
    conv_window=len(b_coeff)
    sr_size=conv_window-1
    num_computation=len(a_coeff)+sr_size
    buff=-1*np.ones(conv_window) ## setting zero vectors
    c_coeff=[]
    input_coeff=list(a_coeff)
    for i in range(sr_size): input_coeff.append(-1)
    
#   print('input coeff:',input_coeff)
    
    for i in range(num_computation):
        for j in range(conv_window-2,-1,-1): buff[j+1]=buff[j]
        buff[0]=input_coeff[i]
#       print('Buffer:',buff)
        
        branch_mul=[]
        for j in range(conv_window):
            dum1=GF_MUL(b_coeff[j], buff[j], order_alpha)
            branch_mul.append(dum1)
#       print('Branch Multip:',branch_mul)
        
        branch_add=-1
        for j in range(len(branch_mul)):
            branch_add=GF_Add(branch_add, branch_mul[j])
        
#        print('Branch add:',branch_add)
        c_coeff.append(branch_add)
        
    return c_coeff


# In[8]:


def Parity_RS(message_RS,generator_RS, order_alpha):
    sr_size=len(generator_RS)-1
    sr_content=-1*np.ones(sr_size)
    for i in range(len(message_RS)-1,-1,-1):
        feedback_add=GF_Add(sr_content[-1],message_RS[i])
        feedback=GF_MUL(feedback_add, generator_RS[-1], order_alpha)
#        print('feedback:',feedback)
        branch_output=[]
        for j in range(sr_size):
            dum2=GF_MUL(feedback, generator_RS[j], order_alpha)
            branch_output.append(dum2)
#        print('branch output:',branch_output) 
        for j in range(sr_size-1,0,-1):
            sr_content[j]=GF_Add(sr_content[j-1], branch_output[j])
        sr_content[0]=branch_output[0]
#        print('sr:',sr_content)
    return sr_content


# In[9]:


def POL_DIV(Dividend,Divisor):
    Dividend_coeff=[]
    for i in range(len(Dividend)): Dividend_coeff.append(Dividend[i])
    Divisor_coeff=[]
    for i in range(len(Divisor)): Divisor_coeff.append(Divisor[i])
    Divisor_coeff[-1]=(-1*Divisor_coeff[-1])%(2**m-1)
#    print("Div coeff:",Divisor_coeff)
    Quotient_coeff=[]
    Remainder_coeff=[]
    sr_size=len(Divisor_coeff)-1
    sr_content=-1*np.ones(sr_size)
    for i in range(len(Dividend)-1,-1,-1):
        feedback=GF_MUL(sr_content[-1],Divisor_coeff[-1], order_alpha)
#        print('feedback:',feedback)
        Quotient_coeff.append(feedback)
        
        branch_output=[]
        for j in range(sr_size):
            dum2=GF_MUL(feedback, Divisor_coeff[j], order_alpha)
            branch_output.append(dum2)
#        print('branch output:',branch_output)
        for j in range(sr_size-1,0,-1):
            sr_content[j]=GF_Add(sr_content[j-1], branch_output[j])
        sr_content[0]=GF_Add(branch_output[0], Dividend_coeff[-1])
        del  Dividend_coeff[-1]
#        print('Shift reg content:',sr_content)
    for i in range(sr_size): Remainder_coeff.append(sr_content[i])
    for i in range(sr_size):
        if Remainder_coeff[-1] < 0: 
            del Remainder_coeff[-1]
        else: break
            
    Quotient_coeff.reverse()
    for i in range(len(Quotient_coeff)):
        if Quotient_coeff[-1] < 0:
            del Quotient_coeff[-1]
        else: break
    
    return (Quotient_coeff,Remainder_coeff) # sr_content = coefficients of remainder #


# In[10]:


def POL_Add(a_coeff, b_coeff):
    max_a_b_length=max(len(a_coeff),len(b_coeff))
    min_a_b_length=min(len(a_coeff),len(b_coeff))
    a_b_GF_addition_coeff=[]
    for i in range(min_a_b_length):
        dum=GF_Add(a_coeff[i],b_coeff[i])
        a_b_GF_addition_coeff.append(dum)
    if len(a_coeff) == len(b_coeff): return a_b_GF_addition_coeff
    if len(a_coeff) > len(b_coeff):
        for i in range(min_a_b_length,max_a_b_length):
            a_b_GF_addition_coeff.append(a_coeff[i])
        return a_b_GF_addition_coeff
    else:
        for i in range(min_a_b_length,max_a_b_length):
            a_b_GF_addition_coeff.append(b_coeff[i])
        return a_b_GF_addition_coeff


# In[11]:


# GFE_Index*j mod Order_GFE, j=0,1,..,Vector_length-1 #  
def GFE_Power_vector(GFE_Index,Order_GFE,Vector_length):
    vector_list=[]
    if GFE_Index < 0: 
        for i in range(Vector_length): vector_list.append(-1)
    else:
        for i in range(Vector_length): vector_list.append((GFE_Index*i)%Order_GFE)
    return vector_list


# In[12]:


# Innder Product of Two GF Vectors #
def GF_InnerProduct(A_GF_Vector,B_GF_Vector):
    Vector_Length=min(len(A_GF_Vector),len(A_GF_Vector))
    GF_Sum=-1
    GF_Product_Vector=[]
    for i in range(Vector_Length):
        GF_Product_Vector.append(GF_MUL(A_GF_Vector[i],B_GF_Vector[i], order_alpha))
    for i in range(Vector_Length):
        GF_Sum=GF_Add(GF_Sum,GF_Product_Vector[i])
    return GF_Sum


# In[13]:


def GF_Function_Value(GFE_Index,Order_GFE,GF_Function_coeff):
    A_GF_Vector=GFE_Power_vector(GFE_Index,Order_GFE,len(GF_Function_coeff))
    Function_Value=GF_InnerProduct(A_GF_Vector,GF_Function_coeff)
    return Function_Value


# In[14]:


# From a0+a1*x+a2*x^2+... ==> To: a1+a2*2*x+... #
def GF_Derivative(GFVector):
    GFVector_Prime=[]
    A_GF_Vector=[]
    for i in range(len(GFVector)):
        if i%2 == 0:
            A_GF_Vector.append(-1)
        else:
            A_GF_Vector.append(0)
    for i in range(len(GFVector)):
        GFVector_Prime.append(GF_MUL(A_GF_Vector[i],GFVector[i], order_alpha))
    
    del GFVector_Prime[0]
    return(GFVector_Prime)

