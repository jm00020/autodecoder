"""
Created on Tue Aug 27 19:46:45 2019
Edited on Mon Sept 2, 10:20, 2019
Tue Sept 10,2019
SCD_AWGN.py
@author: Young Joon Song
"""
import numpy as np
import scipy.special as sp
from matplotlib import pyplot as plt

np.random.seed(103)
n_power=7
CodeRate=0.5
N_length=int(2**n_power)
msg_length=int(CodeRate*N_length)
frozen_length=N_length-msg_length

TotalNumMessage=10
EbNo_Start=1.0  # in dB
EbNo_delta=0.5 # in dB
EbNo_Numbers=10 
EbNo=np.arange(EbNo_Start,EbNo_Start+EbNo_Numbers*EbNo_delta,EbNo_delta)
ebno=10.**(EbNo/10.)          # dB to real value 변환
BER_BPSK_Uncoded=0.5*sp.erfc(np.sqrt(ebno))      # BPSK BER 공식

file=open("SCDAWGN.txt","w")
file.write("BER for polar code over AWGN channel of code length ")
file.write(str(N_length))
file.write("\n")
file.write(", Total number of messages: ")
file.write(str(TotalNumMessage))
file.write("\n")
file.write("EbNo(dB):  BER:   \n ")


def BhanttacharyyaParameterAWGN():
    EbNo_Real=1
    InitialValue=np.exp(-1*CodeRate*EbNo_Real)
    Z_before=[]
    Z_after=[]
    CapacitySplit=np.ones(N_length)
    Z_index_sorted=np.arange(N_length)
    Z_before.append(InitialValue)  # Initial value of Bhattacharyya parameter #
    i=0
    while (i < n_power):
        N_half=len(Z_before)
        for j in range(N_half):
            Z_Value=2*Z_before[j]-Z_before[j]**2
            if Z_Value < 10**-10: Z_Value=10**-10
            elif Z_Value > (1 - 10**-10): Z_Value=1 - 10**-10
            Z_after.append(Z_Value)
            
            Z_Value=Z_before[j]**2
            if Z_Value < 10**-10: Z_Value=10**-10
            elif Z_Value > (1 - 10**-10): Z_Value=1 - 10**-10
            Z_after.append(Z_Value)
#        print("Z:",Z_after)    
        Z_before=Z_after
        Z_after=[]
        i+=1
        
    Z_after=np.array(Z_before)
    CapacitySplit=CapacitySplit-Z_after
    plt.plot(np.arange(N_length),CapacitySplit,'m*')
    plt.xlabel("Bit positions",size=12)
    plt.ylabel("Channel capacity",size=12)
    plt.title('Split Channel Capacity over AWGN EbNo=1')
    plt.grid(True)
    plt.show()
    #######################################################


    Z_after=np.array(Z_before)
    Z_sorted=sorted(Z_before,reverse=True)
    for i in range(N_length):
        Z_index_sorted[i]=Z_before.index(Z_sorted[i])
        Z_before[Z_index_sorted[i]]=-99
    return (Z_after, Z_index_sorted)
   
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

(Z_AWGN,Z_decending_index_AWGN)=BhanttacharyyaParameterAWGN()
print("Z:\n",Z_AWGN)
#print("index sorted:",Z_decending_index_AWGN)

SelectedBitsPositions=np.ones(msg_length)
for i in range(msg_length):
    SelectedBitsPositions[i]=Z_decending_index_AWGN[frozen_length+i]
#print('SelectedBitsPositions',SelectedBitsPositions)
SelectedBitsPositions.sort()
print('Sorted SelectedBitsPositions',SelectedBitsPositions)

FrozenBitsPositions=np.ones(frozen_length)
for i in range(N_length-msg_length):
    FrozenBitsPositions[i]=Z_decending_index_AWGN[i]
#print('Frozen Bits Positions:',FrozenBitsPositions)

FrozenBitsPositions.sort()
print('sorted Frozen Bits Positions:',FrozenBitsPositions)


BERs=[]
for EbNo_cnt, EbNo_real in enumerate(ebno):
    Sigma=np.sqrt(1/(2.0*EbNo_real*CodeRate))
    TotalErrorBits=0
    for _ in range(TotalNumMessage):
        msg=np.random.randint(0,2,size=msg_length)
#        print('message bits:',msg)
        uncodedbits=np.zeros(N_length)
        decodedbits=[]
        for i in range(msg_length):
            uncodedbits[int(SelectedBitsPositions[i])]=msg[i]
    #    print("uncoded bits with frozen",uncodedbits)  
        G=G_All[n_power]
    #    print('Generator Matrix:\n',G)
        codeword=np.matmul(uncodedbits,G)%2
    #    print("polar codeword:",codeword)
        ReceviedAWGN=(-1)**codeword+np.random.normal(0,Sigma,N_length)
        LR_AWGN=np.exp(4*EbNo_real*CodeRate*ReceviedAWGN)
        LLR_AWGN=4*EbNo_real*CodeRate*ReceviedAWGN
#        print('LLR AWGN',LLR_AWGN)
        
        Remain_FBPositions=list(FrozenBitsPositions)
        
        for CurrentBitNumber in range(N_length):
            Binary_CurrentBitNumber=D2B(n_power,CurrentBitNumber)
            length_Remain_FBPositions=len(Remain_FBPositions)
            FrozenBitAssigned = False
            if (length_Remain_FBPositions>0):
                if (Remain_FBPositions[0] == CurrentBitNumber):
                    FrozenBitAssigned = True
                    decodedbits.append(0)
                    del Remain_FBPositions[0]
    #                print("Frozen bit '0' is added and the current decoded bits are:\n",decodedbits)
                    
            if (FrozenBitAssigned==False):
    #            print('index of Bit to be decoded:',CurrentBitNumber)
    #            print('Binary of bit index:',Binary_CurrentBitNumber)
                TreeNodeBitNumbers={}
                for TreeLevelNo in range(n_power):
                    NodeBitSize=Binary_CurrentBitNumber[(n_power-1)-TreeLevelNo]*(2**TreeLevelNo)
                    TreeNodeBitNumbers[TreeLevelNo]=NodeBitSize
    #                print('Tree Bit numbers',TreeNodeBitNumbers)
          
                cnt_TreeAssignedBits=0
                TreeAssignedBits={}
                TreeLevelNo=n_power-1
                while (cnt_TreeAssignedBits < CurrentBitNumber):
                    TreeAssignedBits[TreeLevelNo]=decodedbits[cnt_TreeAssignedBits:cnt_TreeAssignedBits+TreeNodeBitNumbers[TreeLevelNo]]
                    cnt_TreeAssignedBits+=TreeNodeBitNumbers[TreeLevelNo]
                    TreeLevelNo-=1
    #                print('Assigned bits:',TreeAssignedBits)
                
                NodeInputs=LLR_AWGN
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
    #                    print('Node control bits:',NodeControlValues)
                        for NodeCnt in range(NodeNum):
                            LeftInput=NodeInputs[2*NodeCnt]
                            RightInput=NodeInputs[2*NodeCnt+1]
                            Ref_Value=NodeControlValues[NodeCnt]
                            NodeOutputs[NodeCnt]=NodeOutValueRef(LeftInput,RightInput,Ref_Value)
                            
                    NodeInputs=NodeOutputs
        #            print("Node Outputs:",NodeOutputs)
                if (NodeOutputs[0] >= 0):
                    decodedbits.append(0)
                else:
                    decodedbits.append(1)
    #            print('Decoded bits:',decodedbits)
                
        EstimatedMessage=[]
        for i in range(msg_length):
            EstimatedMessage.append(decodedbits[int(SelectedBitsPositions[i])])
#        print("Estimated message is:",EstimatedMessage)
        NumErrors=0
        for i in range(msg_length):
            if (EstimatedMessage[i] != msg[i]): NumErrors+=1
#        print('Error count:',NumErrors)
        TotalErrorBits+=NumErrors
        ber=NumErrors/msg_length
    #    print('ber:',ber)
    BER=TotalErrorBits/(msg_length*TotalNumMessage)
    print('EbNo(dB):',EbNo[EbNo_cnt],'  BER::',BER)
    BERs.append(BER)
    file.write(str(EbNo[EbNo_cnt]));file.write("        "); file.write(str(BER))
    file.write("\n")

print('EbNo:',EbNo)
print('BERs for each EbNo:',BERs)
plt.title('BER over AWGN Channel')
plt.plot(EbNo,BER_BPSK_Uncoded,'r--',EbNo,BERs,'m-')
plt.yscale('log')
plt.xlabel('EbNo(dB)',size=12)
plt.ylabel('BER',size=12)
plt.legend(('Uncoded BPSK','Polar Code'),loc=0, shadow=True)
plt.grid(True)
plt.show()
file.close()
   
