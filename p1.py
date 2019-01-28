import numpy as np
import numpy.linalg as la
#Global Variables
myState2=[(np.sqrt(0.1)*1.j, '101'),(np.sqrt(0.5), '000') ,(-np.sqrt(0.4), '010' )]
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
            output += str(state[i][0]) + "|"+str(int(state[i][1],len(state[i][1])))+">";
        else:
            output += str(state[i][0]) + "|"+str(int(state[i][1],len(state[i][1])))+">"+"+ ";
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
print(stateToVec(myState2))
print(vecToState(stateToVec(myState2)))


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
numberOfWires,myInput = ReadInput('p1/test_input.txt')
print(numberOfWires,myInput)

#Return a list of Matrixs from the input circuit
def CircuitMatrixList(numberOfWires,myInput):
    mats = []
    for gate in myInput:
        if gate[0] == 'H':
            mats.append(HadmardArray(int(gate[1]),numberOfWires))
        elif gate[0] == 'CNOT':
            mats.append(CNOTArray(int(gate[1]),int(gate[2]),numberOfWires))
        else:
            mats.append(PhaseArray(int(gate[1]),numberOfWires,float(gate[2])))
    return mats
matrix_stack = CircuitMatrixList(numberOfWires,myInput)
