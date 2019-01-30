import numpy as np
import numpy.linalg as la
#Global Variables
identity = np.identity(2)
hadmard = (1/np.sqrt(2))*np.array([[1,1],[1,-1]])
CNOT = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
CNOT_flip = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]])
#Vector2State
def prettyPrintBinary(state):
    output="";
    for i in range(len(state)):
        if i == len(state)-1:
            output += str(state[i][0]) + "|"+state[i][1]+">";
        else:
            output += str(state[i][0]) + "|"+state[i][1]+">"+"+ ";
    print(output);
def prettyPrintInteger(state):
    output="";
    for i in range(len(state)):
        if i == len(state)-1:
            output += str(state[i][0]) + "|"+str(int(state[i][1],2))+">";
        else:
            output += str(state[i][0]) + "|"+str(int(state[i][1],2))+">"+"+ ";
    print(output);
def stateToVec(state):
    size = 2**len(state[0][1])
    vec = np.zeros(size,dtype=complex);
    for i in range(len(state)):
        vec[int(state[i][1],2)] = state[i][0]
    return vec;
def vecToState(vec):
    state = []
    for i in range(len(vec)):
        state.append((vec[i],str(bin(i)[2:])))
    return state

#Quantum Circuit Components
#
#
#

#Build the phase matrix for 2-qbit system
def phase(theta):
    j = np.sqrt(-1+0j)
    return np.array([[1,0],[0,np.exp(j*theta)]],dtype = 'complex')

#Return the Hadmard Matrix on the i-th wire with k wires
def HadmardArray(i,k):
    if(i==0):
        if k==1:
            return hadmard
        return np.kron(HadmardArray(i,k-1),identity)
    return np.kron(identity,HadmardArray(i-1,k-1))

#Return the CNOT matrix with Given controlWire & otherWire
def CNOTArray(controlWire,otherWire,totalWire):
    if controlWire < otherWire:
        if controlWire == 0:
            if totalWire == 2:
                return CNOT
            return np.kron(CNOTArray(controlWire,otherWire,totalWire-1),identity)
        return np.kron(identity,CNOTArray(controlWire-1,otherWire,totalWire-1))
    else:
        if otherWire == 0:
            if totalWire == 2:
                return CNOT_flip
            return np.kron(CNOTArray(controlWire,otherWire,totalWire-1),identity)
        return np.kron(identity,CNOTArray(controlWire-1,otherWire-1,totalWire-1))

#Return the Phase matrix on the i-th wire with k-wires and theta
def PhaseArray(i,k,theta):
    if i==0:
        if k==1:
            return phase(theta)
        return np.kron(PhaseArray(i,k-1,theta),identity)
    return np.kron(identity,PhaseArray(i-1,k-1,theta))

#Read Input from a file and return the numberOfWires and input Circuit as a tuple (int,list)
def ReadInput(filename):
    lines = open(filename).readlines()
    myInput = []
    numberOfWires = int(lines[0])
    for line in lines[1:]:
        myInput.append(line.split())
    return (numberOfWires,myInput)


#Return a list of Matrixs from the input circuit
def CircuitMatrixList(numberOfWires,myInput):
    mats = []
    myInitState = []
    for component in myInput:
        if component[0] == 'H':
            mats.append(HadmardArray(int(component[1]),numberOfWires))
        elif component[0] == 'CNOT':
            mats.append(CNOTArray(int(component[1]),int(component[2]),numberOfWires))
        elif component[0] == 'P':
            mats.append(PhaseArray(int(component[1]),numberOfWires,float(component[2])))
        elif component[0] == 'INITSTATE':
            if component[1] == 'BASIS':
                myInitState.append(component[2])
        else:
            myInitState.append(component)

    return (mats,myInitState)


#Return the measurement as a vector
def measure(matrix_stack,stateList):
    mat=matrix_stack.pop()
    while(len(matrix_stack)!=0):
        mat = mat @ matrix_stack.pop()
    if len(stateList) == 1:
        basis = stateList[0].strip('|>')
        vector = np.zeros(2**len(basis),dtype='complex')
        vector[int(basis,2)] = 1.0+0.j
    else:
        vector = np.zeros(len(stateList),dtype='complex')
        for i in range(len(stateList)):
            vector[i] = complex(float(stateList[i][0]),float(stateList[i][1]))
    result = mat @ vector
    for i in range(len(result)):
        result[i] = result[i].conjugate() * result[i]
    return result

#################################################################
#################################################################
#################################################################
#################################################################
#########################TEST CODES##############################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
matrix_stack,stateList = CircuitMatrixList(numberOfWires,myInput)
numberOfWires,myInput = ReadInput('test_input.txt')
prettyPrintBinary(vecToState(measure(matrix_stack,stateList)))
