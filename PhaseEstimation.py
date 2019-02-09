
from mySimulator1c import *
import numpy as np
import numpy.linalg as la
import scipy.sparse as ss
import matplotlib.pyplot as plt
#Global Variables
initState = np.array([0,1,0,0])

circuitList = []
phi = np.linspace(0,2*np.pi,100,endpoint=True)
initState = np.array([0,0,0,1],dtype='complex')
def measure(mat,state):
    mat.reverse()
    while(len(mat)!=0):
        cur = mat.pop()
        state = cur @ state

    for i in range(len(state)):
        state[i] = state[i].conjugate() * state[i]
    return state

def output(numberOfWires,phi_array):
    matList = []
    finalState = []
    for i in phi_array:
        matList.append(HadmardArray(0,2))
        matList.append(CPhaseArray(0,1,2,i))
        matList.append(HadmardArray(0,2))
        finalState.append(measure(matList,initState))
    return finalState
def thetaEstimation(numberOfWires,state):
    base = 2**numberOfWires
    theta = 0
    for i in state:
        theta += i[0] * int(i[1],2)/base
    return theta

finalState = output(2,phi)
theta = []
for i in finalState:
    theta.append(thetaEstimation(2,vecToState(i)))

theta = np.array(theta)
plt.plot(phi/(2*np.pi),theta)
plt.xlabel('$\phi$/2$\pi$')
plt.ylabel('$\Theta$')
plt.title('$\Theta$ Estimation for a 2-Wire Circuit')
plt.show()

phi_2 = 0.1432394487827058*np.pi*2
print(phi_2)
matList = []
matList.append(HadmardArray(0,2))
matList.append(CPhaseArray(0,1,2,0.9))
matList.append(HadmardArray(0,2))
print(thetaEstimation(2,vecToState(measure(matList,initState))))
p2 = plt.plot()
vec = measure(matList,initState)
print(vec)
n,bins,patches=plt.hist(vec,len(vec))


print(n)

plt.plot(n)

plt.show()
