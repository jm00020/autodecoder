import numpy as np

def AllGeneratorMatrixPolarCode(n_power):
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
        g_before=g_after
        GenAll[n_variable]=g_before
    return GenAll

def BhanttacharyyaParameterBEC(N_length, n_power):
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

def LLR_Value_BEC(ReceviedBEC):
    return np.where(ReceviedBEC == 0, 77, np.where(ReceviedBEC == 1, -77, 0))

def D2TreeNodeBitNumber(BitLength,DecimalNum): ###n_power, currentbitnumber
    BinList = []
    BinList.extend(format(DecimalNum,'b').zfill(BitLength))
    BinArray = np.array(BinList).astype(int)
    
    BinPower = [2**i for i in range(BitLength-1,-1,-1)]
    TreeNodeBitNumber = np.flip(BinArray * BinPower)
    
    return TreeNodeBitNumber

def NodeOutValue(L_In,R_In):
    LEFT=np.exp(L_In)
    RIGHT=np.exp(R_In)
    N_Out=np.log(1+LEFT*RIGHT)-np.log(LEFT+RIGHT)
    N_Out = np.where(N_Out >= 77, 77, np.where(N_Out < -77, -77, N_Out))
    return np.round(N_Out,5)

def NodeOutValueRef(L_In,R_In,Ref):
    N_Out= L_In*((-1)**Ref) + R_In
    N_Out = np.where(N_Out >= 77, 77, np.where(N_Out < -77, -77, N_Out))
    return np.round(N_Out,5)

def Encoding(codeNum, N_length, n_power, msg_length, SelectedBitsPositions, G_All):
    msg = np.random.randint(0,2,size=(codeNum, msg_length)) 
    
    uncodedbits = np.zeros((codeNum, N_length))
    uncodedbits[:, SelectedBitsPositions[:msg_length]]=msg
    
    codeword = np.matmul(uncodedbits, G_All)%2
    
    return msg, codeword
    
def ReceivedCode(EraseProb, codeNum, N_length, codeword):
    ReceviedBEC = np.where(np.random.uniform(low=0, high=1, size=(codeNum, N_length)) > EraseProb, codeword, -1)
    
    return ReceivedBEC

def Decoding(codeNum, N_length, n_power, SelectedBitsPositions, FrozenBitsPositions, G_All, LLR_BEC):
    decodedbits = np.zeros((codeNum, N_length))
    decodedbits[:,FrozenBitsPositions] = 0

    for CurrentBitNumber in SelectedBitsPositions:
        TreeNodeBitNumbers = D2TreeNodeBitNumber(n_power,CurrentBitNumber)
    
        cnt_TreeAssignedBits = 0
        TreeAssignedBits = {}
        TreeLevelNo = n_power-1
        while (cnt_TreeAssignedBits < CurrentBitNumber):
            TreeAssignedBits[TreeLevelNo] = decodedbits[:,cnt_TreeAssignedBits:cnt_TreeAssignedBits+TreeNodeBitNumbers[TreeLevelNo]]
            cnt_TreeAssignedBits += TreeNodeBitNumbers[TreeLevelNo]
            TreeLevelNo -= 1
        
        NodeInputs = LLR_BEC
        for TreeLevelNo in range(n_power-1,-1,-1):
            NodeNum = int(2**TreeLevelNo)
            NodeOutputs = np.ones((codeNum, NodeNum))
            if TreeNodeBitNumbers[TreeLevelNo] == 0:
                for NodeCnt in range(NodeNum):
                    LeftInput = NodeInputs[:,2*NodeCnt]
                    RightInput = NodeInputs[:,2*NodeCnt+1]
                    NodeOutputs[:,NodeCnt] = NodeOutValue(LeftInput, RightInput)
            else :
                NodeControlValues = np.matmul(TreeAssignedBits[TreeLevelNo],G_All[TreeLevelNo])%2
                for NodeCnt in range(NodeNum):
                    LeftInput = NodeInputs[:,2*NodeCnt]
                    RightInput = NodeInputs[:,2*NodeCnt+1]
                    Ref_Value = NodeControlValues[:,NodeCnt]
                    NodeOutputs[:,NodeCnt] = NodeOutValueRef(LeftInput, RightInput, Ref_Value)
            NodeInputs = NodeOutputs

        decodedbits[NodeOutputs[:,0] >= 0, CurrentBitNumber] = 0
        decodedbits[NodeOutputs[:,0]  < 0, CurrentBitNumber] = 1
    
    EstimatedMessage = decodedbits[:,SelectedBitsPositions]
    
    return EstimatedMessage

def BERCalculation(EstimatedMessage, msg, codeNum, msg_length):
    error_index = np.where(EstimatedMessage - msg != 0)
    NumErrors = error_index[0].shape[0] 
    BER = NumErrors/(codeNum*msg_length)
    
    return BER