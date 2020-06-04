import numpy as np
import scipy.special as sp
from math import gcd

class GF_function:
    def __init__(self, GF_type):
        self.GF_type = GF_type
        if self.GF_type == 32:
            from Source.GF32 import GFE
            self.GFE = GFE
            self.create_parameter()

        elif self.GF_type == 64:
            from Source.GF64 import GFE
            self.GFE = GFE
            self.create_parameter()

    def create_parameter(self):
        if self.GF_type == 32:         
            self.m = 5  
            self.t_RS = 2
        elif  self.GF_type == 64:
            self.m = 6
            self.t_RS = 3

        self.order_alpha = 2**self.m-1 
        
        b_RS = 1
        alpha_power = 1
        self.code_length = int(self.order_alpha/gcd(self.order_alpha, alpha_power))
        self.Z_RS = []
        for i in range(2*self.t_RS): self.Z_RS.append((b_RS+i)%self.order_alpha)

        Roots_BCH=[]
        for i in range(len(self.Z_RS)):
            Init_Power=self.Z_RS[i]
            Roots_BCH.append(Init_Power)
            while 1:
                Next_power=Init_Power*2%self.order_alpha
                if Next_power == self.Z_RS[i]:
                    break
                else:
                    Roots_BCH.append(Next_power)
                Init_Power=Next_power
            Roots_BCH=list(set(Roots_BCH))

        a_x = [0]
        a_x.append(Roots_BCH[0])
        b_x = [0]
        for i in range(1, len(Roots_BCH)):
            b_x.append(Roots_BCH[i])
            c_x = self.POL_MUL(a_x,b_x)
            a_x = c_x
            b_x = [0]
        self.g_x = c_x
        self.g_x.reverse()
        pairty_length=len(Roots_BCH)
        self.msg_length=int(self.code_length-pairty_length)

    def encoding(self, msg_size):
        msg_RS = np.random.randint(-1, self.order_alpha, size=(msg_size, self.msg_length))
        parity_symbols = self.Parity_RS(msg_RS, self.g_x, msg_size)
        codeword_RS = np.hstack((parity_symbols, msg_RS))
        return msg_RS, codeword_RS
    
    def send_code(self, codeword_RS):
        send_code = self.GFE[(codeword_RS+1).astype(np.int)].reshape(codeword_RS.shape[0], self.order_alpha*self.m)
        return np.where(send_code > 0, -1, 1)

    def add_noise(self, EbNo, code):
        CodeRate = self.msg_length/self.code_length
        ebno = 10**(EbNo/10)
        Sigma = np.sqrt(1/(2*ebno*CodeRate))
        noise = np.random.normal(loc=0,scale=Sigma, size=(code.shape[0],code.shape[1]))
        noise_code = code + noise
        return noise_code

    def syndrome_check(self, ):
        syndrome_RS = -1*np.ones(2*self.t_RS)
        for i in range(2*self.t_RS):
            for j in range(len(codeword_RS)):
                deum1 = GF_MUL(r_vector_RS[j], self.Z_RS[i]*j)
                syndrome_RS[i] = GF_Add(dum1, syndrome_RS[i])
        syndrome_check = 0
        for i in range(2*self.t_RS):
            if syndrome_RS[i] >=0:
                syndrome_check+=1

    def GF_Index(self, cc):
        equal_check = (self.GFE + np.array(cc))%2
        index = np.argmin(equal_check.sum(axis=1))
        return index-1

    def GF_Add(self, a,b): # alpha^c = alpha^a + alpha^b
        a = np.array(a)
        b = np.array(b)
        aa=self.GFE[(a+1).astype(np.int)]
        bb=self.GFE[(b+1).astype(np.int)]
        cc=np.add(aa,bb)%2
        cc = np.reshape(cc, (-1,self.m))
        index = []
        for i in range(cc.shape[0]):
            index.append(self.GF_Index(cc[i]))
        return index

    def GF_MUL(self, a,b):
        a = np.array([a])
        c = -1*np.ones(a.shape)
        index = np.where(a>=0)

        if b >= 0:
            c[index] = (a[index]+b)%self.order_alpha
            
        return c
        
        
#         if a.all() >= 0 and b>=0: # if a and b are greater than -1, alpha^a * alpha^b =alpha^(c), c=a+b (mod n)
#             c=(a+b)%self.order_alpha
#             return c
#         else:
#             return np.array(-1)  # if a and/or b is minus, the result is zero vector that is, alpha^-1

    def GF_MUL_Inv(self, a): # alpha^a * alpha^b = alpha^0
        b=-1*a%self.order_alpha    
        return b


    def POL_MUL(self, a_coeff,b_coeff):
        conv_window=len(b_coeff)
        sr_size=conv_window-1
        num_computation=len(a_coeff)+sr_size
        buff=-1*np.ones(conv_window) ## setting zero vectors
        c_coeff=[]
        input_coeff=list(a_coeff)
        for i in range(sr_size): input_coeff.append(-1)
    
#       print('input coeff:',input_coeff)
    
        for i in range(num_computation):
            for j in range(conv_window-2,-1,-1): buff[j+1]=buff[j]
            buff[0]=input_coeff[i]
#           print('Buffer:',buff)
        
            branch_mul=[]
            for j in range(conv_window):
                dum1=self.GF_MUL(b_coeff[j], buff[j])
                branch_mul.append(dum1)
#           print('Branch Multip:',branch_mul)
        
            branch_add=-1
            for j in range(len(branch_mul)):
                branch_add=self.GF_Add(branch_add, branch_mul[j])
        
#           print('Branch add:',branch_add)
            c_coeff.extend(branch_add)
        
        return c_coeff

    def POL_DIV(self, Dividend,Divisor):
        Dividend_coeff=[]
        for i in range(len(Dividend)): Dividend_coeff.append(Dividend[i])
        Divisor_coeff=[]
        for i in range(len(Divisor)): Divisor_coeff.append(Divisor[i])
        Divisor_coeff[-1]=(-1*Divisor_coeff[-1])%(2**self.m-1)
#    print("Div coeff:",Divisor_coeff)
        Quotient_coeff=[]
        Remainder_coeff=[]
        sr_size=len(Divisor_coeff)-1
        sr_content=-1*np.ones(sr_size)
        for i in range(len(Dividend)-1,-1,-1):
            feedback=self.GF_MUL(sr_content[-1],Divisor_coeff[-1])
#        print('feedback:',feedback)
            Quotient_coeff.append(feedback)
        
            branch_output=[]
            for j in range(sr_size):
                dum2=self.GF_MUL(feedback, Divisor_coeff[j])
                branch_output.append(dum2)
#        print('branch output:',branch_output)
            for j in range(sr_size-1,0,-1):
                sr_content[j]=self.GF_Add(sr_content[j-1], branch_output[j])
            sr_content[0]=self.GF_Add(branch_output[0], Dividend_coeff[-1])
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

    def POL_Add(self, a_coeff, b_coeff):
        max_a_b_length=max(len(a_coeff),len(b_coeff))
        min_a_b_length=min(len(a_coeff),len(b_coeff))
        a_b_GF_addition_coeff=[]
        for i in range(min_a_b_length):
            dum=self.GF_Add(a_coeff[i],b_coeff[i])
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

    def Parity_RS(self,message_RS,generator_RS, msg_size):
        sr_size=len(generator_RS)-1
        sr_content=-1*np.ones((msg_size, sr_size))
        for i in range(message_RS.shape[1]-1,-1,-1):
            feedback_add=self.GF_Add(sr_content[:,-1],message_RS[:,i])
            feedback=self.GF_MUL(feedback_add, generator_RS[-1])
#        print('feedback:',feedback)
            for j in range(sr_size):
                dum2=self.GF_MUL(feedback, generator_RS[j])
                if j == 0:
                    branch_output = dum2.reshape(-1,1)
                else:
                    branch_output = np.hstack((branch_output, dum2.reshape(-1,1)))
#            print('branch output:',branch_output)
            for j in range(sr_size-1,0,-1):
                sr_content[:,j]=self.GF_Add(sr_content[:,j-1], branch_output[:,j])
            sr_content[:,0]=branch_output[:,0]
#        print('sr:',sr_content)
        return sr_content

 
# GFE_Index*j mod Order_GFE, j=0,1,..,Vector_length-1 #          
    def GFE_Power_vector(self, GFE_Index,Order_GFE,Vector_length):
        vector_list=[]
        if GFE_Index < 0: 
            for i in range(Vector_length): vector_list.append(-1)
        else:
            for i in range(Vector_length): vector_list.append((GFE_Index*i)%Order_GFE)
        return vector_list


# Innder Product of Two GF Vectors #
    def GF_InnerProduct(self,A_GF_Vector,B_GF_Vector):
        Vector_Length=min(len(A_GF_Vector),len(A_GF_Vector))
        GF_Sum=-1
        GF_Product_Vector=[]
        for i in range(Vector_Length):
            GF_Product_Vector.append(self.GF_MUL(A_GF_Vector[i],B_GF_Vector[i]))
        for i in range(Vector_Length):
            GF_Sum=self.GF_Add(GF_Sum,GF_Product_Vector[i])
        return GF_Sum

    def GF_Function_Value(self,GFE_Index,Order_GFE,GF_Function_coeff):
        A_GF_Vector=self.GFE_Power_vector(GFE_Index,Order_GFE,len(GF_Function_coeff))
        Function_Value=self.GF_InnerProduct(A_GF_Vector,GF_Function_coeff)
        return Function_Value
    
   
# From: a0+a1*x+a2*x^2+... ==> To: a1+a2*2*x+... #
    def GF_Derivative(self, GFVector):
        GFVector_Prime=[]
        A_GF_Vector=[]
        for i in range(len(GFVector)):
            if i%2 == 0:
                A_GF_Vector.append(-1)
            else:
                A_GF_Vector.append(0)
        for i in range(len(GFVector)):
            GFVector_Prime.append(self.GF_MUL(A_GF_Vector[i],GFVector[i]))
    
        del GFVector_Prime[0]
        return(GFVector_Prime)
