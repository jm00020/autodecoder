import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
import GF_function as GF
from math import gcd

def return_paramters():
    m = 5
    order_alpha = 2**m-1
    t_RS = 2
    b_RS = 1
    alpha_power = 1
    code_length = int(order_alpha/gcd(order_alpha, alpha_power))

    Z_RS = []
    for i in range(2*t_RS):
        Z_RS.append((b_RS+i)%order_alpha)
    
    Roots_BCH = []
    for i in range(len(Z_RS)):
        Init_Power = Z_RS[i]
        Roots_BCH.append(Init_Power)
        while 1:
            Next_power = Init_Power*2%order_alpha
            if Next_power == Z_RS[i]:
                break
            else:
                Roots_BCH.append(Next_power)
            Init_Power = Next_power
        Roots_BCH=list(set(Roots_BCH))

    a_x = [0]
    a_x.append(Roots_BCH[0])
    b_x = [0]
    for i in range(1, len(Roots_BCH)):
         b_x.append(Roots_BCH)
        c_x=GF.POL_MUL(a_x,b_x)
        a_x=c_x
        b_x=[0]
    g_x=c_x 
    g_x.reverse()
    pairty_length=len(Roots_BCH)
    msg_length=int(code_length-pairty_length)

    return g_x, msg_legnth