"""
Created on Thu Sep 19 08:25:17 2019
SCD_BEC.py
@author: Young Joon Song
"""

import numpy as np
from matplotlib import pyplot as plt

np.random.seed(103)
n_power=7
CodeRate=0.5
N_length=int(2**n_power)
msg_length=int(CodeRate*N_length)
frozen_length=N_length-msg_length
TotalNumMessage=4

Prob_start=0.15
Prob_delta=0.05
BEC_count=8
BECs=np.arange(Prob_start,Prob_start+BEC_count*Prob_delta,Prob_delta)


file=open("SCDBEC.txt","w")
file.write("BER for polar code over BEC channel of code length ")
file.write(str(N_length))
file.write(", Total number of messages: ")
file.write(str(TotalNumMessage))
file.write("\n")
file.write("Delta:  BER:   \n ")


def BhanttacharyyaParameterBEC():
    prob_erase=0.5
    Z_before=[]
    Z_after=[]
    CapacitySplit=np.ones(N_length)
    Z_index_sorted=np.arange(N_length)
    Z_before.append(prob_erase)  # Initial value of Bhattacharyya parameter #
    i=0
    while (i < n_power):
        N_half=len(Z_before)
        for j in range(N_half):
            Z_Value = 2*Z_before[j] - Z_before[j]**2
            if Z_Value < 10**-30: 
                Z_Value=10**-30
            elif Z_Value >= (1 - 10**-30): 
                Z_Value=1 - 10**-30
            Z_after.append(Z_Value)
            
            Z_Value=Z_before[j]**2
            if Z_Value < 10**-30: 
                Z_Value=10**-30
            elif Z_Value >= (1 - 10**-30): 
                Z_Value=1 - 10**-30
            Z_after.append(Z_Value)
#        print("Z:",Z_after)    
        Z_before=Z_after
        Z_after=[]
        i+=1
    
    Z_sorted=sorted(Z_before,reverse=True)
    for i in range(N_length):
        Z_index_sorted[i]=Z_before.index(Z_sorted[i])
        Z_before[Z_index_sorted[i]]=-99
    return (Z_after, Z_index_sorted)


def LLR_Value_BEC(R_bit):
    return{0:77,1:-77,-1:0}[R_bit]
    
def AllGeneratorMatrixPolarCode():
    GenAll={}
    g_before=np.ones([1,1])
    GenAll[0]=g_before
    n_variable=0
    while (n_variable < n_power):
        n_col_before=2**n_variable # number of columns for g_before matrix
        n_row_before=2**n_variable # number of rows for g_before matrix
        n_variable+=1
        g_after=np.ones([2**n_variable,2**n_variable])
        for i in range(n_col_before):
            g_after[0:n_col_before,2*i]=g_before[:,i]
            g_after[0:n_col_before,2*i+1]=np.zeros(n_col_before)
            g_after[n_row_before::,2*i]=g_before[:,i]
            g_after[n_row_before::,2*i+1]=g_before[:,i]
#        print(g_after)
        g_before=g_after
        GenAll[n_variable]=g_before
    return GenAll

'''
from Decimal to Bits in List type of lenght BitLength.
From the left to right bit positios, the order is from MSB to LSB.
'''
def D2B(BitLength,DecimalNum): 
    BinList=list(bin(DecimalNum))
    del BinList[0:2]
    Len_BinList=len(BinList)
    for i in range(Len_BinList): BinList[i]=int(BinList[i]) 
    if Len_BinList>BitLength:
        print('Not enough BitLength , DecimalNum is too big. I will return empty List []. Try again!!\n')
        return []
    elif Len_BinList==BitLength: return BinList
    else:
        diff_len=BitLength-Len_BinList
        while(diff_len>0):
            BinList.insert(0,0)
            diff_len-=1
        return BinList

def B2D(Bin_vector): ## binary to decimal conversion
    Decimal_num=0
    for i in range(len(Bin_vector)):
        Decimal_num += Bin_vector[i]*2**(len(Bin_vector)-1-i)
    return int(Decimal_num)
    
def NodeOutValue(L_In,R_In):
    LEFT=np.exp(L_In)
    RIGHT=np.exp(R_In)
    N_Out=np.log(1+LEFT*RIGHT)-np.log(LEFT+RIGHT)
    if N_Out >= 77: N_Out=77
    elif N_Out < -77: N_Out=-77
    return round(N_Out,5)

def NodeOutValueRef(L_In,R_In,Ref):
    N_Out=((-1)**Ref)*L_In + R_In
    if N_Out >= 77: N_Out=77
    elif N_Out < -77: N_Out=-77
    return round(N_Out,5)

G_All=AllGeneratorMatrixPolarCode()
#print('All Gen:\n',G_All)

(Z_BEC,Z_decending_index_BEC)=BhanttacharyyaParameterBEC()
#print("Z:\n",Z_BEC)
print("index sorted:",Z_decending_index_BEC)

SelectedBitsPositions=np.ones(msg_length)
for i in range(msg_length):
    SelectedBitsPositions[i]=Z_decending_index_BEC[frozen_length+i]
#print('SelectedBitsPositions',SelectedBitsPositions)
SelectedBitsPositions.sort()
print('Sorted SelectedBitsPositions',SelectedBitsPositions)

FrozenBitsPositions=np.ones(frozen_length)
for i in range(N_length-msg_length):
    FrozenBitsPositions[i]=Z_decending_index_BEC[i]
#print('Frozen Bits Positions:',FrozenBitsPositions)

FrozenBitsPositions.sort()
print('sorted Frozen Bits Positions:',FrozenBitsPositions)

print('Generator Matrix:\n',G_All[n_power])

BERs=[]
for EraseProb_cnt, EraseProb in enumerate(BECs):
    TotalErrorBits=0
    for MessageNo in range(TotalNumMessage):
        msg=np.random.randint(0,2,size=msg_length)
     #   print('message bits:',msg)
        uncodedbits=np.zeros(N_length)
        decodedbits=[]
       
        for i in range(msg_length):
            uncodedbits[int(SelectedBitsPositions[i])]=msg[i]
     #   print("uncoded bits with frozen",uncodedbits)  
        
        codeword=np.matmul(uncodedbits,G_All[n_power])%2
     #   print("polar codeword:",codeword)
    
        ReceviedBEC=[]
        for i in range(N_length):
            if(np.random.uniform(low=0,high=1) > EraseProb):
                ReceviedBEC.append(codeword[i])
            else:
                ReceviedBEC.append(-1)
     #   print('Received BEC:',ReceviedBEC)
        ''' -1 means earased '''
        LLR_BEC=np.ones(N_length)
        for i in range(N_length):
            LLR_BEC[i]=LLR_Value_BEC(ReceviedBEC[i])
     #   print('LLR BEC',LLR_BEC)
        Remain_FBPositions=list(FrozenBitsPositions)
    #    print("Remaing FB positions",Remain_FBPositions)
        
        for CurrentBitNumber in range(N_length):
            Binary_CurrentBitNumber=D2B(n_power,CurrentBitNumber)
            length_Remain_FBPositions=len(Remain_FBPositions)
            FrozenBitAssigned = False
            if CurrentBitNumber in Remain_FBPositions:
                FrozenBitAssigned = True
                decodedbits.append(0)
                Remain_FBPositions.remove(CurrentBitNumber)
      #          print("Remaing FB positions",Remain_FBPositions)
        
            if (FrozenBitAssigned==False):
       #         print('Decoded bits until now:',decodedbits)
        #        print('Bit No. to be decoded:',CurrentBitNumber)
         #       print('Binary of bit index:',Binary_CurrentBitNumber)
                TreeNodeBitNumbers={}
                for TreeLevelNo in range(n_power):
                    NodeBitSize=Binary_CurrentBitNumber[(n_power-1)-TreeLevelNo]*(2**TreeLevelNo)
                    TreeNodeBitNumbers[TreeLevelNo]=NodeBitSize
          #      print('Tree Bit numbers',TreeNodeBitNumbers)
                cnt_TreeAssignedBits=0
                TreeAssignedBits={}
                TreeLevelNo=n_power-1
                while (cnt_TreeAssignedBits < CurrentBitNumber):
                    TreeAssignedBits[TreeLevelNo]=decodedbits[cnt_TreeAssignedBits:cnt_TreeAssignedBits+TreeNodeBitNumbers[TreeLevelNo]]
                    cnt_TreeAssignedBits+=TreeNodeBitNumbers[TreeLevelNo]
                    TreeLevelNo-=1
           #     print('Assigned bits:',TreeAssignedBits)
                NodeInputs=LLR_BEC
                for TreeLevelNo in range(n_power-1,-1,-1):
                    NodeNum=int(2**TreeLevelNo)
                    NodeOutputs=np.ones(NodeNum)
                    if (TreeNodeBitNumbers[TreeLevelNo] == 0):
                        for NodeCnt in range(NodeNum):
                            LeftInput=NodeInputs[2*NodeCnt]
                            RightInput=NodeInputs[2*NodeCnt+1]
                            NodeOutputs[NodeCnt]=NodeOutValue(LeftInput,RightInput)
                    else:
                        NodeControlValues=np.matmul(TreeAssignedBits[TreeLevelNo],G_All[TreeLevelNo])%2
        #                print('Node control bits at level no:',str(TreeLevelNo),'is',NodeControlValues)
                        for NodeCnt in range(NodeNum):
                            LeftInput=NodeInputs[2*NodeCnt]
                            RightInput=NodeInputs[2*NodeCnt+1]
                            Ref_Value=NodeControlValues[NodeCnt]
                            NodeOutputs[NodeCnt]=NodeOutValueRef(LeftInput,RightInput,Ref_Value)
                    NodeInputs=NodeOutputs
       #             print("Node Outputs:",NodeOutputs)
                if (NodeOutputs[0] >= 0):
                    decodedbits.append(0)
                else:
                    decodedbits.append(1)
        #        print('Decoded bits:',decodedbits)
        EstimatedMessage=[]
        for i in range(msg_length):
            EstimatedMessage.append(decodedbits[int(SelectedBitsPositions[i])])
      #  print("Estimated message is:",EstimatedMessage)
        NumErrors=0
        for i in range(msg_length):
            if (EstimatedMessage[i] != msg[i]): NumErrors+=1
    #    print('Error count:',NumErrors)
        TotalErrorBits+=NumErrors
        ber=NumErrors/msg_length
    #    print('ber:',ber)
    BER=TotalErrorBits/(msg_length*TotalNumMessage)
    print('Erasure Prob:',EraseProb,'  BER:',BER)
    BERs.append(BER)
    file.write(str(BECs[EraseProb_cnt]));file.write("        "); file.write(str(BER))
    file.write("\n")
print('Erasure Probs. over BECs:',BECs)
print('BERs for each BEC:',BERs)

plt.title('BER over BEC')
plt.plot(BECs,BERs,'m*')
plt.yscale('log')
plt.xlabel('Prob. of Erasure',size=12)
plt.ylabel('BER',size=12)
plt.legend(('Polar Code',),loc=0, shadow=True)
plt.grid(True)
plt.show()
file.close()
'''
plt.legend() 
location string and corresponding location number   
'best' 0 
'upper right' 1 
'upper left' 2 
'lower left' 3 
'lower right' 4 
'right' 5 
'center left' 6 
'center right' 7 
'lower center' 8 
'upper center' 9 
'center' 10 
'''

