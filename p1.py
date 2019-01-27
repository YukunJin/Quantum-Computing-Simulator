import numpy as np
import numpy.linalg as la
#Global Variables
myState2=[(np.sqrt(0.1)*1.j, '101'),(np.sqrt(0.5), '000') ,(-np.sqrt(0.4), '010' )]
identity = np.identity(2)
hadmard = (1/np.sqrt(2))*np.array([[1,1],[1,-1]])
CNOT = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
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
    vec = np.zeros(len(state),dtype=complex);
    for i in range(len(state)):
        vec[i] = state[i][0]
    return vec;
prettyPrintBinary(myState2);
prettyPrintInteger(myState2);
print(len(stateToVec(myState2)));

#Quantum Circuit
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
    #TODO: Figure out case otherWire>controlWire
    if controlWire == 0:
        if totalWire == 2:
            return CNOT
        return np.kron(CNOTArray(controlWire,otherWire,totalWire-1),identity)
    return np.kron(identity,CNOTArray(controlWire-1,otherWire-1,totalWire-1))

#Return the Phase matrix on the i-th wire with k-wires and theta
def PhaseArray(i,k,theta):
    if i==0:
        if k==1:
            return phase(theta)
        return np.kron(PhaseArray(i,k-1,theta),identity)
    return np.kron(identity,PhaseArray(i-1,k-1,theta))
