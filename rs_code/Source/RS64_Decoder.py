
# coding: utf-8

# In[ ]:


import numpy as np
from Source.GF64 import GFE
from Source import GF_FUNC as GF_F


# In[ ]:


def Syndrome_Check(t_RS, r_vector_RS, order_alpha, Z_RS):
    syndrome_RS=-1*np.ones(2*t_RS)
    for i in range(2*t_RS):
        for j in range(len(r_vector_RS)): 
            dum1 = GF_F.GF_MUL(r_vector_RS[j], Z_RS[i]*j, order_alpha)
            syndrome_RS[i] = GF_F.GF_Add(dum1, syndrome_RS[i])[0]

    syndrome_check = 0
    for i in range(2*t_RS):
        if syndrome_RS[i] >= 0:
            syndrome_check+=1
    
    return syndrome_RS, syndrome_check


# In[ ]:

def Receive_Bit(Eb_No, cx, k_RS, n_RS):
    req_Eb_No = 10**(Eb_No/10)
    Es_No = req_Eb_No * (k_RS/n_RS)
    Sigma=np.sqrt(1/(2*Es_No))
        
    rx = cx + np.random.normal(loc=0.0,scale=Sigma,size=cx.shape)
    
    return rx

def Trans_symbol(rx, m):
    r_RS = []
    rx_Hard = np.where(rx > 0.5, 1, 0)
    for i in range(0, len(rx), m):
        r_RS.append(GF_F.GF_Index(rx_Hard[i:i+m]))
        
    return r_RS

def SER_Calculation(Estimated_code, msg_RS, length_parity):
    Estimated_msg = Estimated_code[:,length_parity:]
    index = np.where((Estimated_msg - msg_RS) != 0)
    
    error_count = index[0].shape[0]
    
    SER = error_count/(msg_RS.shape[0]*msg_RS.shape[1])
    
    return SER
    

def Decoder(syndrome_RS, syndrome_check, r_vector_RS, t_RS, order_alpha, m):
    if syndrome_check == 0:
        return r_vector_RS
    else:
        A_coeff=[]
        for i in range(2*t_RS): 
            A_coeff.append(-1) 
        A_coeff.append(0)
        
        B_coeff=list(syndrome_RS)
        
        REMAINDER=[]
        REMAINDER.append(A_coeff)
        REMAINDER.append(B_coeff)
    
        T_coefficient=[]
        T_coefficient.append([-1])  ## alpha^-1 = 0 ##
        T_coefficient.append([0])  ## alpha^0 = 1 ##
        
        QUOTIENT=[]
        cnt_i=0
        while len(REMAINDER[-1])-1 >= t_RS:
            QUO, REM = GF_F.POL_DIV(REMAINDER[-2],REMAINDER[-1], order_alpha, m)
            QUOTIENT.append(QUO)
            REMAINDER.append(REM)
            Dummy_A=T_coefficient[-2]
            Dummy_B=GF_F.POL_MUL(QUOTIENT[-1],T_coefficient[-1], order_alpha)
            Dummy_C=GF_F.POL_Add(Dummy_A,Dummy_B)
            T_coefficient.append(Dummy_C)
            cnt_i+=1
        
        Kappa=(-1*T_coefficient[-1][0])%order_alpha 
        Error_Location_coeff=[]
        for i in range(len(T_coefficient[-1])):
            Dummy=T_coefficient[-1][i]
            Error_Location_coeff.append(GF_F.GF_MUL(Kappa,Dummy, order_alpha))
    
        for i in range(len(T_coefficient[-1])):
            if T_coefficient[-1][-1] < 0:
                del T_coefficient[-1][-1]
            else: 
                break
                
        Error_Evaluation_coeff=[]
        for i in range(len(REMAINDER[-1])):
            Dummy=REMAINDER[-1][i]
            Error_Evaluation_coeff.append(GF_F.GF_MUL(Kappa,Dummy, order_alpha))
            
        error_positions=[]
        Chien_Search_Values=[]
        num_errors=0
        for i in range(order_alpha):
            alpha_minus_i=(-1*i)%order_alpha ## alpha^(-i mod order_alpha), i=0,1,...,order_alpha ##
            Error_Location_Function_Value=GF_F.GF_Function_Value(alpha_minus_i,order_alpha,Error_Location_coeff)
            Chien_Search_Values.append(Error_Location_Function_Value)
            if Error_Location_Function_Value == -1: 
                error_positions.append(i)
                num_errors+=1

    ## Error Evaluation after finding error locations ##
    # Derivative of Error Location Polynomial #
        Error_Location_Der=GF_F.GF_Derivative(Error_Location_coeff, order_alpha)
    
    # Numerator of Error Evaluation function #
        Error_Numerator=[]
        Error_Denominator=[]
        Error_Estimated_Values=[]
        Estimated_Error_Vector=-1*np.ones(len(r_vector_RS)) # Error Vector #
        for i in range(num_errors):
            Inv_Xi_Index=(-1*error_positions[i])%order_alpha
            Numerator_i=GF_F.GF_Function_Value(Inv_Xi_Index,order_alpha,Error_Evaluation_coeff)
            Error_Numerator.append(Numerator_i)
            Denominator_i=GF_F.GF_Function_Value(Inv_Xi_Index,order_alpha,Error_Location_Der)
            Error_Denominator.append(Denominator_i)
            Estim_Error_Value_i=GF_F.GF_MUL(Numerator_i,GF_F.GF_MUL_Inv(Denominator_i, order_alpha), order_alpha)
            Error_Estimated_Values.append(Estim_Error_Value_i)
            Estimated_Error_Vector[error_positions[i]]=Estim_Error_Value_i
        

    # Estimated Codeword #
        Estimated_Codeword=GF_F.POL_Add(r_vector_RS, Estimated_Error_Vector)
        
        return Estimated_Codeword

