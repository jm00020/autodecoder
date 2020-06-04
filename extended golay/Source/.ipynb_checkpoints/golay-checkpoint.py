import numpy as np
import pickle
import matplotlib.pyplot as plt
from itertools import repeat

P = np.array([
    [1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1],
    [1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1],
    [1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1],
    [1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1],
    [1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0],
    [0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0],
    [0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1],
    [1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0],
    [0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0],
    [0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1],
    [1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0],
    [0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1],
], dtype='int')

def ex_syndrome_decoding(codes):
    with open('syndrome.pickle', 'rb') as syn:
        syndrome = pickle.load(syn)
        
    H = np.concatenate((np.eye(11, dtype="int"), P.T), axis=1)
    fix_code = []
    for code in codes:
        code_syndrome = tuple(np.dot(code[1:], H.T) % 2)
        error_code = syndrome[code_syndrome]
        if sum(error_code) == 3:
            if code[0] == 1:
                fix_code.append((code[1:] + error_code) % 2)
            elif code[0] == 0:
                fix_code.append(list(repeat(-1, 23)))
        elif not error_code:
            fix_code.append(list(repeat(-1, 23)))
        else:
            fix_code.append((code[1:] + error_code) % 2)
    return np.array(fix_code)

def syndrome_decoding(codes):
    with open('syndrome.pickle', 'rb') as syn:
        syndrome = pickle.load(syn)
        
    H = np.concatenate((np.eye(11, dtype="int"), P.T), axis=1)
    fix_code = []
    for code in codes:
        code_syndrome = tuple(np.dot(code, H.T) % 2)
        error_code = syndrome[code_syndrome]
        fix_code.append((code + error_code) % 2)
    return np.array(fix_code)

def make_G():
    G = np.concatenate((P, np.eye(12, dtype='int')), axis=1)
    return G
