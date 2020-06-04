import numpy as np

def G_matrix(length, m, r):
    G = np.ones(length)
    for i in range(m):
        v = np.zeros((int(length/(2**(i+1)))))
        v = np.hstack((v, np.ones((int(length/(2**(i+1)))))))
        while v.shape[0] < length :
            v = np.hstack((v, np.zeros((int(length/(2**(i+1)))))))
            v = np.hstack((v, np.ones((int(length/(2**(i+1)))))))
        G = np.vstack((G,v))
    if r == 1:
        return G
    elif r == 2:
        for i in range(1,m):
            for j in range(i+1,m+1):
                G = np.vstack((G,(G[i]*G[j])))
        return G
    return G

def G_prime(G,m):
    G_p = np.flip(G[1:,:],0)
    return G_p

def Encoding(msg, G, Eb_No, length):
    send_code = np.dot(msg,G)%2
    send_code = np.where(send_code > 0, -1, 1)

    Eb_No = 10**(Eb_No/10)
    Es_No = Eb_No*(msg.shape[1]/length)
    noise_power = np.sqrt(1/(2*Es_No))
#     noise_power = 0

    received_code = send_code + (noise_power * np.random.normal(size = send_code.shape))

    return send_code, received_code

def FHT(send_code, m, length, r, masking=None):
    if r == 1:
        recived_code = send_code
        for j in range(m):
            index = np.array(range(length))
            for l in range(int(length/2)):
                k = index[0]
                
                sum_code = recived_code[:,k] + recived_code[:,k+1*(2**j)]
                subtract_code = recived_code[:,k] - recived_code[:,k+1*(2**j)]
                recived_code[:,k] = sum_code
                recived_code[:,k+1*(2**j)] = subtract_code
                
                index = np.delete(index, [np.where(index == k), np.where(index == k+1*(2**j))])
    
        recived_index = np.zeros((2, send_code.shape[0]))
        max_index = np.absolute(recived_code).argmax(axis=1)
        value = recived_code[range(recived_code.shape[0]), max_index]
        value = np.where(value > 0, 0, 1)
        
        recived_index[0] = value
        recived_index[1] = max_index
        
        return recived_index
    
    elif r == 2:
        max_value = np.zeros((send_code.shape[0]))
        recived_index = np.zeros((3, send_code.shape[0]))
        
        for i in range(masking.shape[0]):
            recived_code = send_code * masking[i]
            
            for j in range(m):
                index = np.array(range(length))
                for l in range(int(length/2)):
                    k = index[0]
                
                    sum_code = recived_code[:,k] + recived_code[:,k+1*(2**j)]
                    subtract_code = recived_code[:,k] - recived_code[:,k+1*(2**j)]
                    recived_code[:,k] = sum_code
                    recived_code[:,k+1*(2**j)] = subtract_code
                
                    index = np.delete(index, [np.where(index == k), np.where(index == k+1*(2**j))])
    
            current_recived_index = np.zeros((3, send_code.shape[0]))
            current_max_value = np.absolute(recived_code).max(axis=1)
            max_index = np.absolute(recived_code).argmax(axis=1)
            value = recived_code[range(recived_code.shape[0]), max_index]
            value = np.where(value > 0, 0, 1)
            
            current_recived_index[0] = value
            current_recived_index[1] = max_index
            current_recived_index[2] = i
            
            max_check = max_value - current_max_value
            change_index = np.where(max_check < 0)
            max_value[change_index] = current_max_value[change_index]
            recived_index[:, change_index] = current_recived_index[:, change_index]
            
        return recived_index
            
    return recived_index

def Decoding(send_code, m, length, r, G, G_p):
    if r == 1:
        index = FHT(send_code, m, length, r)
        Estimated_code = np.vstack((index[0].reshape(1, send_code.shape[0]), np.flip(G_p[:, index[1].astype(int)], 0)))
        
        return Estimated_code.T
    
    elif r == 2:
        masking_len = int(np.math.factorial(m)/(np.math.factorial(2)*np.math.factorial(m-2)))
        masking_list = []
        for i in range(2**masking_len):
            masking_list.extend(format(i, 'b').zfill(masking_len+(m+1)))
        masking_list = np.reshape(masking_list,(-1,masking_len+(m+1))).astype(int)
        masking_code = np.dot(masking_list,G)%2
        masking_code = np.where(masking_code > 0, -1, 1)                        
        
        index = FHT(send_code, m, length, r, masking_code)
        Estimated_code = np.vstack((index[0].reshape(1, send_code.shape[0]), np.flip(G_p[:, index[1].astype(int)], 0)))
        masking = masking_list[index[2].astype(int)]
        Estimated_code[m+1:,:] = masking[:,m+1:].T
        
        return Estimated_code.T
    
    return Estimated_code
    