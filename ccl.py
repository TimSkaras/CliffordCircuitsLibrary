# canonical ordering of symplectic group elements
# from "How to efficiently select an arbitrary clifford group element"
#       by Robert Koenig and John A. Smolin
#

import numpy as np
import math
import matplotlib.pyplot as plt
import time
import random, copy, itertools
import qiskit.quantum_info as qisq

########### Smolin Code Begin #############

def directsum(m1,m2):
    n1=len(m1[0])
    n2=len(m2[0])
    output = np.zeros((n1+n2,n1+n2),dtype=np.int8)
    for i in range(0,n1):
        for j in range(0,n1):
            output[i,j]=m1[i,j]
    for i in range(0,n2):
        for j in range(0,n2):
            output[i+n1,j+n1]=m2[i,j]
    return output
    ######### end directsum

def inner(v,w):  # symplectic inner product
    t=0
    for i in range(0, np.size(v)>>1):
        t+=v[2*i]*w[2*i+1]
        t+=w[2*i]*v[2*i+1]
    return t%2

def transvection(k,v): # applies transvection Z_k to v
    return (v+inner(k,v)*k)%2

def int2bits(i,n):  
    # converts integer i to an length n array of bits
    # rightmost bit is the most significant bit
    # This is the opposite convention used in the simulator part of code
    output=np.zeros(n,dtype=np.int8)
    for j in range(0,n):
        output[j]=i&1
        i>>=1
    return output

def findtransvection(x,y):  # finds h1,h2 such that y = Z_h1 Z_h2 x
    # Lemma 2 in the text
    # Note that if only one transvection is required output[1] will be
    #         zero and applying the all-zero transvection does nothing.
    output=np.zeros((2,np.size(x)),dtype=np.int8)
    if np.array_equal(x,y):
        return output
    if inner(x,y)==1:
        output[0]=(x+y)%2
        return output
    #
    # find a pair where they are both not 00
    z=np.zeros(np.size(x))
    for i in range(0,np.size(x)>>1):
        ii=2*i
        if ((x[ii]+x[ii+1]) != 0) and ((y[ii]+y[ii+1]) != 0):  # found the pair
            z[ii]=(x[ii]+y[ii])%2
            z[ii+1]=(x[ii+1]+y[ii+1])%2
            if (z[ii]+z[ii+1])==0:  # they were the same so they added to 00
                z[ii+1]=1
                if x[ii]!=x[ii+1]:
                    z[ii]=1
            output[0]=(x+z)%2
            output[1]=(y+z)%2
            return output
    # didn't find a pair
    # so look for two places where x has 00 and y doesn't, and vice versa
    #
    # first y==00 and x doesn't
    for i in range(0,np.size(x)>>1):
        ii=2*i
        if ((x[ii]+x[ii+1]) != 0) and ((y[ii]+y[ii+1]) == 0):  # found the pair
            if x[ii]==x[ii+1]:
                z[ii+1]=1
            else:
                z[ii+1]=x[ii]
                z[ii]=x[ii+1]
            break
    #
    # finally x==00 and y doesn't
    for i in range(0,np.size(x)>>1):
        ii=2*i
        if ((x[ii]+x[ii+1]) == 0) and ((y[ii]+y[ii+1]) != 0):  # found the pair
            if y[ii]==y[ii+1]:
                z[ii+1]=1
            else:
                z[ii+1]=y[ii]
                z[ii]=y[ii+1]
            break
    output[0]=(x+z)%2
    output[1]=(y+z)%2
    return output
    ###################### end findtransvction

################################################################################
def symplectic(i,n): # output symplectic canonical matrix i of size 2nX2n
    ################################################################################
    #    Note, compared to the text the transpose of the symplectic matrix
    #    is returned.  This is not particularly important since
    #                 Transpose(g in Sp(2n)) is in Sp(2n)
    #    but it means the program doesn't quite agree with the algorithm in the
    #    text. In python, row ordering of matrices is convenient, so it is used
    #    internally, but for column ordering is used in the text so that matrix
    #    multiplication of symplectics will correspond to conjugation by
    #    unitaries as conventionally defined Eq. (2).  We can't just return the
    #    transpose every time as this would alternate doing the incorrect thing
    #    as the algorithm recurses.
    #
    nn=2*n   # this is convenient to have
    # step 1
    s=((1<<nn)-1)
    k=(i%s)+1
    i//=s
    #
    # step 2
    f1=int2bits(k,nn)
    #
    # step 3
    e1=np.zeros(nn,dtype=np.int8) # define first basis vectors
    e1[0]=1
    T=findtransvection(e1,f1) # use Lemma 2 to compute T
    #
    # step 4
    # b[0]=b in the text, b[1]...b[2n-2] are b_3...b_2n in the text
    bits=int2bits(i%(1<<(nn-1)),nn-1)
    #
    # step 5
    eprime=np.copy(e1)
    for j in range(2,nn):
        eprime[j]=bits[j-1]
    h0=transvection(T[0],eprime)
    h0=transvection(T[1],h0)
    #
    # step 6
    if bits[0]==1:
        f1*=0
     # T' from the text will be Z_f1 Z_h0.  If f1 has been set to zero
     #                                       it doesn't do anything
     # We could now compute f2 as said in the text but step 7 is slightly
     # changed and will recompute f1,f2 for us anyway
    #
   # step 7
    #  define the 2x2 identity matrix
    id2=np.zeros((2,2),dtype=np.int8)
    id2[0,0]=1
    id2[1,1]=1
    #
    if n!=1:
        g=directsum(id2,symplectic(i>>(nn-1),n-1))
    else:
        g=id2
    #
    for j in range(0,nn):
        g[j]=transvection(T[0],g[j])
        g[j]=transvection(T[1],g[j])
        g[j]=transvection(h0,g[j])
        g[j]=transvection(f1,g[j])
    #
    return g
    ############# end symplectic

def bits2int(b,nn):  
    # converts an nn-bit string b to an integer between 0 and 2^n-1
    # Assumes right most bit is most significant
    output=0
    tmp=1
    for j in range(0,nn):
        if b[j]==1:
           output=output+tmp
        tmp=tmp*2
    return output

def numberofcosets(n): # returns  the number of different cosets
    x=int(np.power(2,2*n-1))*int(np.power(2,2*n)-1)
    return x;

################################################################################
def symplecticinverse(n,gn): # produce an index associated with group element gn
    ################################################################################

    nn=2*n   # this is convenient to have

    # step 1
    v=gn[0];
    w=gn[1];

    # step 2
    e1=np.zeros(nn,dtype=np.int8) # define first basis vectors
    e1[0]=1
    T=findtransvection(v,e1); # use Lemma 2 to compute T

    # step 3
    tw=np.copy(w)
    tw=transvection(T[0],tw)
    tw=transvection(T[1],tw)
    b=tw[0];
    h0=np.zeros(nn,dtype=np.int8)
    h0[0]=1
    h0[1]=0
    for j in range(2,nn):
        h0[j]=tw[j]


    # step 4
    bb=np.zeros(nn-1,dtype=np.int8)
    bb[0]=b;
    for j in range(2,nn):
        bb[j-1] =tw[j];
    zv=bits2int(v,nn)-1; # number between 0...2^(2n)-2
                                    # indexing non-zero bitstring v of length 2n

    zw=bits2int(bb,nn-1);  # number between 0..2^(2n-1)-1
                                         #indexing w (such that v,w is symplectic pair)
    cvw=int(zw*(np.power(2, 2*n)-1)+zv);
    # cvw is a number indexing the unique coset specified by (v,w)
    # it is between 0...2^(2n-1)*(2^(2n)-1)-1=numberofcosets(n)-1


    #step 5
    if n==1:
        return cvw

    #step 6
    gprime=np.copy(gn);
    if b==0:
        for j in range(0,nn):
            gprime[j]=transvection(T[1],transvection(T[0],gn[j]))
            gprime[j]=transvection(h0,gprime[j])
            gprime[j]=transvection(e1,gprime[j])
    else:
        for j in range(0,nn):
            gprime[j]=transvection(T[1],transvection(T[0],gn[j]))
            gprime[j]=transvection(h0,gprime[j])

    # step 7
    gnew=gprime[2:nn,2:nn]; # take submatrix
    gnidx=symplecticinverse(n-1,gnew)*numberofcosets(n)+cvw;
    return gnidx
    ####### end symplecticinverse

########### Smolin Code End #############

def smolin2natural(basestr):
    """
    In the base string representation, there are two possible conventions: smolin and natural

    SMOLIN CONVENTION:
    00 -- 0 -- I
    01 -- 1 -- Z
    10 -- 2 -- X
    11 -- 3 -- Y

    NATURAL CONVENTION:
    00 -- 0 -- I
    01 -- 1 -- X
    10 -- 2 -- Y
    11 -- 3 -- Z
    """

    base_copy = np.array(basestr)
    out = copy.deepcopy(base_copy)
    out[base_copy == 1] = 3
    out[base_copy == 3] = 2
    out[base_copy == 2] = 1

    return out.tolist()

def int2basestr(n, b, l=0):
    # The leftmost digit is the most significant digit

    if n < b:
        a = [n]
    else:
        a = int2basestr(n // b, b)
        a.append(n % b)
    if l > len(a):
        a = [0] * (l - len(a)) + a
    return a

def basestr2int(st, b):
    # Assumes the left most digit is the most significant 
    return np.sum([st[-i-1] * b**i for i in range(len(st))])

def smolinstr2basestr(st):
    return smolin2natural(int2basestr(basestr2int(st, 2), 4, len(st)//2))

def pystr2smolinstr(pystr):
    # Mapping of Pauli operators to their Smolin bit representations
    translation = {'I': [0, 0], 'Z': [0, 1], 'X': [1, 0], 'Y': [1, 1]}
    
    # Remove leading '+' or '-' if present
    if pystr[0] == '+' or pystr[0] == '-':
        pystr = pystr[1:]
    
    smolin_bits = []
    for pauli in pystr:
        # Append the corresponding bits for each Pauli operator
        smolin_bits.extend(translation[pauli])
    
    return smolin_bits

def basestr2pystr(pauli):
    """
    Convert array of pauli indices (possibly complex) to python string format

    E.g.,   [0,-1]      --> -IX
            [1,2,-3,-3] -->  IYZZ
    """

    n = len(pauli)

    pauli_dict = ["I", "X", "Y", "Z"]
    out_str = ""

    for i in range(n):
        out_str += pauli_dict[pauli[i]]

    return out_str

def expToPauli(exp1, exp2):
    """
    Converts exp format UXU^dag = X^alpha Z^beta to pauli str index

    INPUTS:
    exp1 -- int that is either 0 or 1
    exp2 -- int that is either 0 or 1

    OUTPUT:
    int giving index of pauli -- 0, 1,2, or 3

    """

    exps = [exp1, exp2]

    if exps == [0,0]:
        return 0
    elif exps == [0, 1]:
        return 3
    elif exps == [1,0]:
        return 1
    elif exps == [1,1]:
        return 2
    else:
        print("invalid input")

def num_clifford_circuits(n):
    # Number of clifford circuits on n qubits modulo pauli gates

    prod = 1

    for i in range(1, n+1):
        prod *= 4**i - 1

    return 2**(n**2) * prod

####### Simulator Code #######

def row_sum(row1, row2):
    """
    Sums two rows in aaronson format. Note this is equivalent to multiplying to Pauli strings
    so additional logic must be done to make sure the right sign is computed

    If there are n qubits, there should be 2*n + 1 elements in each row (last element is for the sign)

    """

    G = 0
    n = len(row1)//2
    out = [0]*len(row1)
    for i in range(n):
        x1 = row1[i]
        x2 = row2[i]
        z1 = row1[i+n]
        z2 = row2[i+n]
        G += Pauli.g(x1,z1,x2,z2)
        out[i] = (x1+x2) % 2
        out[i+n] = (z1+z2) % 2

    G = (G + 2*row1[-1] + 2*row2[-1]) % 4
    out[-1] = G//2

    return out

def inverse(A):
    """
    Finds inverse of symplectic matrix A with binary finite field

    """
    n = len(A)//2

    # Create matrix [A | I]
    augmented_A = np.zeros([2*n, 4*n+1], dtype=int)
    augmented_A[:,:2*n] = A
    for i in range(2*n): augmented_A[i, 2*n+i] = 1
    augmented_rref = rref(augmented_A)

    # Use Gaussian elimination on reduced rows
    for i in range(2*n):
        for j in range(i):
            if augmented_rref[j,i] == 1:
                augmented_rref[j,:] = (augmented_rref[j,:] + augmented_rref[i,:]) % 2

    return augmented_rref[:,2*n:-1]

def rref(A):
    """
    Computes rref of an aaronson stabilizer tableau
    """

    pivot_count = 0
    for j in range(A.shape[1]):
        rows = np.where(A[:,j]==1)[0].tolist()

        # If column is all zero, just continue; nothing to reduce
        if len(rows) == 0: continue

        # we do not reduce rows that have already been pivoted
        while len(rows) > 0 and rows[0] < pivot_count: rows.pop(0)

        # we must set first unpivoted row to have pivot column
        # if this row does not have a 1, put a 1 there by switching rows
        if len(rows) >0 and rows[0] != pivot_count: 
            temp = copy.copy(A[pivot_count]) 
            A[pivot_count] = A[rows[0]]
            A[rows[0]] = temp
            rows[0] = pivot_count

        if len(rows) > 0:
            for c in rows[1:]:
                A[c] = row_sum(A[c], A[rows[0]])  #add zero'th to c'th row
            pivot_count += 1

    return A

class Clifford:
    """
    As the name suggests, objects of this class represent a single clifford in "compiled" form,
    where compiled form means a single matrix that directly states how basis vectors in symplectic
    space get mapped to cliffords. These objects are supposed to be representing a single smolin 
    matrix. If you want to compose many clifford together, use the CliffordCircuit class

    """

    def __init__(self, cliff_idx : int, sign_idx : int, n : int, s = None):
        """
        Creates object representing Clifford circuit for simulation with symplectic matrices

        INPUTS:
        cliff_idx -- integer index of desired symplectic matrix under smolin convention
        sign_idx -- integer between 0 and 2**(2*n) (see notes below)
        n -- number of qubits
        s -- symplectic matrix if it has already been computed

        NOTES:
        - sign_idx gives the sign factor on each column in the symplectic matrix
        - the symplectic matrix gives the transformation of X_1 and Z_1 in the first columns
        - the first two bits (from the left) of sign_idx tell whether to add a negative sign for these two paulis
        e.g., if sign_idx = [0,1,1,0] for a 2 qubit clifford this means Z_1 and X_2 pick up a negative sign under conjugation
        """

        if s is None: s = np.transpose(symplectic(cliff_idx, n)).astype(int)
        r = np.array(int2basestr(sign_idx, 2, 2*n))

        self.cliff_idx = cliff_idx
        self.sign_idx = sign_idx
        self.s = s
        self.r = r
        self.num_qubits = n

    def __repr__(self):
        """
        Returns a string representation of the Clifford object,
        displaying the symplectic matrix 's' and the sign array 'r',
        with each element of 'r' aligned under the corresponding column of 's'.
        A delimiter line separates the symplectic matrix and the sign array.
        """
        s = self.s
        r = self.r
        n = self.num_qubits

        # Define the width for each element to ensure alignment
        col_width = 2  # Adjust as needed for larger integers

        # Build formatted rows for the symplectic matrix 's'
        formatted_s_rows = []
        for row in s:
            formatted_row = ' '.join(f"{int(elem):>{col_width}}" for elem in row)
            formatted_s_rows.append(formatted_row)

        # Build formatted string for the sign array 'r'
        formatted_r = ' '.join(f"{int(elem):>{col_width}}" for elem in r)

        # Create delimiter line
        total_width = len(formatted_s_rows[0])  # Total width of a matrix row
        delimiter = '-' * total_width

        # Combine the formatted symplectic matrix, delimiter, and sign array
        matrix_str = '\n'.join(formatted_s_rows)
        result = f"{matrix_str}\n{delimiter}\n{formatted_r}"

        return result

    @classmethod
    def random_clifford(cls, n : int):
        """
        Produces random Clifford object sampled uniformly from all n qubit clifford circuits
        """

        num_circuits = num_clifford_circuits(n)
        cliff_idx = random.randint(0, num_circuits-1)
        sign_idx = random.randint(0, 4**n -1)

        return Clifford(cliff_idx, sign_idx, n)

    def inverse(self):
        """
        Outputs ccl Clifford object that is the inverse unitary of self
        """

        s = inverse(self.s)
        cliff_idx = symplecticinverse(self.num_qubits, np.transpose(s))
        sign_idx = self.sign_idx

        return Clifford(cliff_idx, sign_idx, self.num_qubits, s)

class Pauli:

    def __init__(self, smolin_vec, exp = 0):

        self.smolin_vec = np.array(smolin_vec) % 2
        self.exp = exp % 4
        self.num_qubits = len(smolin_vec) // 2

    def __repr__(self):

        return self.exp2str(self.exp) + basestr2pystr(smolinstr2basestr(self.smolin_vec))

    @staticmethod
    def g(x1, z1, x2, z2):
        """
        g function from aaronson gottesman paper

        This function returns the exponent of i that will remain
        after multiplying two paulis P_1 * P_2

        INPUTS:
        x1 -- 1 if Pauli 1 has an X in it (X or Y) else 0
        z1 -- 1 if Pauli 1 has an Z in it (Z or Y) else 0
        x2 -- 1 if Pauli 2 has an X in it (X or Y) else 0
        z2 -- 1 if Pauli 2 has an Z in it (Z or Y) else 0
        """

        if x1 == 0 and z1 == 0: return 0
        if x1 == 1 and z1 == 1: return z2 - x2
        if x1 == 1 and z1 == 0: return z2 * (2*x2 - 1)
        else: return x2 * (1 - 2*z2)

    @staticmethod
    def g_factor(pauli1, pauli2):
        """
        g factor from multiplying pauli1 * pauli2
        
        INPUTS:
        pauli1 -- binary array representing pauli in smolin format 
        pauli2 -- binary array representing pauli in smolin format 

        OUTPUTS:  
        exp -- resulting exp on imaginary i
        """

        G = 0
        n = len(pauli1)//2
        for i in range(n):
            x1, z1 = pauli1[2*i], pauli1[2*i+1]
            x2, z2 = pauli2[2*i], pauli2[2*i+1]

            G += Pauli.g(x1, z1, x2, z2)

        return G

    @staticmethod
    def exp2str(G):
        """
        converts exponent to string suitable for printing
        """

        if G % 4 == 0: return "+"
        if G % 4 == 1: return "i"
        if G % 4 == 2: return "-"
        if G % 4 == 3: return "-i"

    @classmethod
    def from_pystr(cls, pystr):

        if pystr[0] == "-":
            exp = 2
            return Pauli(pystr2smolinstr(pystr), exp)
        else:
            return Pauli(pystr2smolinstr(pystr), 0)

    @classmethod    
    def random_pauli(cls, n):
        """
        Creates random Pauli object on n qubits
        """

        smolin_vec = np.random.randint(0,2,2*n)
        exp = np.random.randint(0,2)*2

        return Pauli(smolin_vec, exp)

    def multiply(self, pauli):
        """
        Multiply object by pauli on the right P_1 P_2, where P_1 is self

        INPUTS:
        pauli -- pauli object for right multiplication

        OUTPUTS:
        resulting Pauli as Pauli object
        """ 

        smolin_vec = (self.smolin_vec + pauli.smolin_vec ) % 2
        exp = (self.exp + pauli.exp + self.g_factor(self.smolin_vec, pauli.smolin_vec) ) % 4

        return Pauli(smolin_vec, exp)

    def conjugate(self, circuit):
        """
        Conjugate Pauli object with  circuit given

        INPUTS:
        circuit -- Pauli or Clifford object

        Returns new Pauli object with correct coefficient
        """
        
        if isinstance(circuit, Clifford):
            P = Pauli([0]*2*self.num_qubits, self.exp)

            # if Pauli is just identity, conjugation leaves it unchanged
            if not 1 in self.smolin_vec: return P

            #increment exp for each Y in smolin_vec
            for i in range(self.num_qubits):
                if self.smolin_vec[2*i] * self.smolin_vec[2*i + 1] == 1: P.exp += 1

            # loop through paulis to multiply together
            ones = np.where(self.smolin_vec == 1)[0] # indices of ones in smolin_vec
            pauli_arr = [Pauli(circuit.s[:,ones[i]], 2*circuit.r[ones[i]]) for i in range(len(ones))]
            
            while len(pauli_arr) > 0:
                P = P.multiply(pauli_arr[0])
                pauli_arr.pop(0)

            return P

        elif isinstance(circuit, Pauli):
            
            if self.commutes(circuit):
                return self
            else:
                P = Pauli(self.smolin_vec, self.exp)
                P.exp = (P.exp + 2) % 4
                return P

    def commutes(self, pauli):
        """
        Returns True if self commutes with Pauli object pauli, else returns False
        """

        if self.num_qubits != pauli.num_qubits:
            print("Paulis differ in length")
            return False

        s = 0
        for i in range(self.num_qubits):    
            s += self.smolin_vec[2*i] * pauli.smolin_vec[2*i + 1]
            s += self.smolin_vec[2*i+1] * pauli.smolin_vec[2*i]

        return not bool(s % 2)

class StabilizerState:

    def __init__(self, pauli_arr, smolin_arr = None, exps=None, generators = None, destabilizers = None):
        """
        Takes array of Paulis specified as strings to create StabilizerState

        INPUTS:
        pauli_arr -- python array of strings, e.g., pauli_arr = ["+IX", "-ZI"]
        smolin_arr -- array where each column is a pauli in smolin format
        exps -- coefficients associated with each column in smolin_arr
        generators -- array of Pauli objects where each is a generator for the stabilizer state
        destabilizers -- python array of Pauli objects where each Pauli anticommutes with exactly one generator (the one with the same index)

        NOTE: 
        - Paulis in pauli arr must have +/- signs; no imaginary i's
        - Each element of exps represents the exponents on imaginary i
          e.g., [0, 0, 2] --> i**0, i**0, i**2, --> +, +, -
        - Primary method for initialization is via array of python strings
        - If you only have smolin_arr, you can still initialize by just passing None for other arguments
        """

        r_func = lambda x : 2 if x == '-' else 0
        if smolin_arr is not None and not isinstance(smolin_arr, np.ndarray): smolin_arr = np.array(smolin_arr)

        if not generators is None: 
            self.generators = generators
        
        elif smolin_arr is None:
            self.generators = [Pauli(pystr2smolinstr(pauli_arr[i]), r_func(pauli_arr[i][0]) ) for i in range(len(pauli_arr))]

        elif exps is None:
            self.generators = [Pauli(smolin_arr[:,i], 0) for i in range(smolin_arr.shape[1])]

        else:
            self.generators = [Pauli(smolin_arr[:,i], exps[i]) for i in range(smolin_arr.shape[1])]
        

        if not destabilizers is None:
            destabs = copy.deepcopy(destabilizers)
            self.destabilizers = destabs
        else:
            print("Destabilizers must be provided. Code for finding destabilizers is broken")
            self.destabilizers = self.get_destabilizers(self.generators)


        self.num_qubits = self.generators[0].num_qubits

    def __repr__(self):
        destabs_str = ', '.join([str(paul) for paul in self.destabilizers])
        gens_str = ', '.join([str(paul) for paul in self.generators])
        return f"Destabilizers: {destabs_str}\nGenerators: {gens_str}"


    @classmethod
    def zero_state(cls, n):
        # Prepares StabilizerState object of zero state on n qubits

        smolin_arr = np.zeros([2*n, n], dtype=int)
        destabs_arr = np.zeros([2*n, n], dtype=int)
        for i in range(n): 
            smolin_arr[2*i + 1, i] = 1
            destabs_arr[2*i, i] = 1

        generators = [Pauli(smolin_arr[:,i]) for i in range(n)]
        destabilizers = [Pauli(destabs_arr[:,i]) for i in range(n)]
        return StabilizerState(None, generators=generators, destabilizers=destabilizers)

    @classmethod
    def get_destabilizers(cls, generators):
        """
        TODO: Fix this code. It produces destabilizers that destabilize more than one generator and are not mutually commuting

        Finds destabilizers for given generators

        INPUTS:
        generators -- array of mutually commuting Pauli objects

        OUTPUTS:
        array of Pauli objects that are detabilizers for generators
        """

        n = generators[0].num_qubits
        smolin_arr = np.zeros([n, 2*n], dtype=int)
        destabilizers = []

        for i in range(n):
            # Generate random Pauli
            h = Pauli(np.random.randint(0,2, 2*n))

            # Keep generating Paulis until we find one that
            # anticommutes with exactly one
            while generators[i].commutes(h) or not 1 in h.smolin_vec: 
                h = Pauli(np.random.randint(0,2, 2*n))

                # make sure it commutes with previous destabs and generators
                for j in range(len(destabilizers)):
                    if not destabilizers[j].commutes(h): h = h.multiply(generators[i])
                    if not generators[j].commutes(h): h = h.multiply(destabilizers[j])

            destabilizers.append(h)

        # # Check that destabilizers commute with exactly one generator
        for i in range(n-2,-1, -1):
            # for j in range(i+1,n):
                # if not destabilizers[i].commutes(generators[j]): destabilizers[i] = destabilizers[i].multiply(destabilizers[j])
            destabilizers[i].exp = 0

        return destabilizers

    def get_tableau(self):
        """
        Returns numpy array with stabilizers formated in gottesman aaronson format
        """

        stab_arr = np.array([g.smolin_vec for g in self.generators])
        signs = np.array([g.exp//2 for g in self.generators]).reshape(-1,1)

        tableau = np.zeros([self.num_qubits, 2*self.num_qubits], dtype=int)
        tableau[:,:self.num_qubits] = stab_arr[:, ::2]
        tableau[:,self.num_qubits:] = stab_arr[:, 1::2]

        tableau = np.hstack((tableau, signs))

        return tableau

    def evolve(self, circuit):
        """
        Evolves stabilizer state under conjugation of circuit 
        
        INPUT:
        clifford -- ccl Clifford object specifying clifford circuit or Pauli object specifying Pauli for conjugation

        OUTPUT:
        Evolved StabilizerState
        """
        generators = [g.conjugate(circuit) for g in self.generators]
        destabilizers = [g.conjugate(circuit) for g in self.destabilizers]
        stab = StabilizerState(None, generators = generators, destabilizers = destabilizers)

        return stab

    def get_probability(self, bitstring : str):
        """
        Compute the probability of outcome given in bitstring

        INPUTS:
        bitstring -- python string specifying outcome of Z basis measurement, e.g., "0010"

        OUTPUT:
        float between 0 and 1 giving probability
        """

        # Create clifford for converting bitstring to all zeros
        bitflips = np.array(list(bitstring)).astype(int)
        smolin_vec = np.zeros(2* self.num_qubits, dtype= int)
        smolin_vec[::2] = bitflips

        # These are the bitflips we need to make to convert this problem
        # to an inner product with the zero state
        P = Pauli(smolin_vec)
        stab = self.evolve(P)

        # Now do gaussian elimination to find inner product with zero state
        rref_tab = rref(stab.get_tableau())
        x_tab = rref_tab[:,:stab.num_qubits]
        z_tab_signs = rref_tab[:, stab.num_qubits:]
        
        x_rank = 0
        while x_rank < stab.num_qubits and 1 in x_tab[x_rank]: x_rank+=1

        if x_rank == stab.num_qubits: return 2**(-1*stab.num_qubits)
        
        # if the z generators have a negative sign, probability is zero
        if 1 in z_tab_signs[x_rank:,-1]: return 0.0

        else: return 2 ** (-1 * x_rank)

    def measure(self, qubits = None):
        """
        By default measures all qubits and returns the resulting bit string and resulting StabilizerState object

        INPUTS:
        qubits -- array of indices to measure (if None all qubits are measured)

        OUPUTS:
        bitstring -- python string of bits giving measurement outcome
        result -- StabilizerState object representing measured state
        """

        if qubits is None: index = range(self.num_qubits)
        else: index = qubits

        psi = copy.deepcopy(self)

        bitstring = ""

        for i in index:

            p = -1
            for j in range(psi.num_qubits): 
                if psi.generators[j].smolin_vec[2*i] == 1: 
                    p = j
                    break

            # Case 1: Z_i is not in stabilizer group
            if p >= 0:
                for j in range(psi.num_qubits):
                    if psi.generators[j].smolin_vec[2*i] == 1 and p != j: psi.generators[j] = psi.generators[j].multiply(psi.generators[p])
                    if psi.destabilizers[j].smolin_vec[2*i] == 1: psi.destabilizers[j] = psi.destabilizers[j].multiply(psi.generators[p])
                    psi.destabilizers[j].exp = 0

                # replace destabilizer
                psi.destabilizers[p] = copy.copy(psi.generators[p])
                psi.destabilizers[p].exp = 0
                
                # prepare new generator
                Zp = [0]*2*psi.num_qubits
                Zp[2*i + 1] = 1
                outcome = np.random.randint(0,2)
                psi.generators[p] = Pauli(Zp, 2*outcome)

                # record measurement outcome
                bitstring += str(outcome)
            
            # Case 2: Z_i is in stabilizer group
            else:
                scratch = Pauli([0]*2*psi.num_qubits)
                
                for j in range(psi.num_qubits):
                    if psi.destabilizers[j].smolin_vec[2*i] == 1: scratch = scratch.multiply(psi.generators[j])

                outcome = scratch.exp//2
                bitstring += str(outcome)

        return bitstring, psi

    def to_qiskit(self):
        # outputs qiskit object version of StabilizerState
        
        return qisq.StabilizerState.from_stabilizer_list([str(g) for g in self.generators])

class CliffordCircuit:
    """
    A class to represent a Clifford circuit composed of individual Clifford gates.
    This class allows building and manipulating Clifford circuits by applying basic
    Clifford gates and composite Clifford operations to qubits and subsets of qubits.
    """

    def __init__(self, num_qubits):
        """
        Initialize a CliffordCircuit with a specified number of qubits.

        INPUTS:
        num_qubits -- int
            The number of qubits in the circuit.

        OUTPUTS:
        None
        """
        self.num_qubits = num_qubits

        self.stab = StabilizerState.zero_state(num_qubits)

    def h(self, qubit):
        """
        Apply a Hadamard gate to the specified qubit.

        INPUTS:
        qubit -- int
            The index of the qubit to which the Hadamard gate is applied.

        OUTPUTS:
        None
        """
        n = self.num_qubits

        # Symplectic matrix for H gate
        S = np.identity(2 * n, dtype=int)
        # Swap X and Z for qubit 'qubit'
        S[2*qubit, 2*qubit], S[2*qubit+1, 2*qubit+1] = 0,0
        S[2*qubit, 2*qubit+1] , S[2*qubit+1, 2*qubit] = 1,1

        # Create Clifford object
        cliff = Clifford(0, 0, n, s=S)

        # Evolve stabilizer state
        self.stab = self.stab.evolve(cliff)

    def s(self, qubit):
        """
        Apply an S gate (phase gate) to the specified qubit.

        INPUTS:
        qubit -- int
            The index of the qubit to which the S gate is applied.

        OUTPUTS:
        None
        """
        n = self.num_qubits

        # Symplectic matrix for S gate
        S = np.identity(2 * n, dtype=int)

        # Update X column
        S[2 * qubit, 2*qubit], S[2 * qubit+1, 2*qubit] = 1,1 

        # Create Clifford object
        cliff = Clifford(0, 0, n, s=S)

        # Evolve stabilizer state
        self.stab = self.stab.evolve(cliff)

    def s_dagger(self, qubit):
        """
        Apply an S† (S-dagger) gate to the specified qubit.

        INPUTS:
        qubit -- int
            The index of the qubit to which the S† gate is applied.

        OUTPUTS:
        None
        """
        n = self.num_qubits

        # Symplectic matrix for S† gate (same as S gate)
        S = np.identity(2 * n, dtype=int)

        # Update X column
        S[2 * qubit, 2*qubit], S[2 * qubit+1, 2*qubit] = 1,1 

        # X --> -Y, so we need to put a sign on X_i
        r = np.zeros(2 * n, dtype=int)
        r[2*qubit] = 1

        # Create Clifford object
        cliff = Clifford(0, 0, n, s=S)
        cliff.r = r

        # Evolve stabilizer state
        self.stab = self.stab.evolve(cliff)

    def x(self, qubit):
        """
        Apply an X gate (Pauli-X gate) to the specified qubit.

        INPUTS:
        qubit -- int
            The index of the qubit to which the X gate is applied.

        OUTPUTS:
        None
        """
        n = self.num_qubits

        # Smolin for X gate (identity matrix)
        S = np.identity(2 * n, dtype=int)

        # Sign array: Applying X gate flips the sign of Z operators
        r = np.zeros(2 * n, dtype=int)
        r[2 * qubit + 1] = 1 # Flip sign of Z

        # Create Clifford object
        cliff = Clifford(0, 0, n, s=S)
        cliff.r = r

        # Evolve stabilizer state
        self.stab = self.stab.evolve(cliff)

    def y(self, qubit):
        """
        Apply a Y gate (Pauli-Y gate) to the specified qubit.

        INPUTS:
        qubit -- int
            The index of the qubit to which the Y gate is applied.

        OUTPUTS:
        None
        """
        n = self.num_qubits

        # Symplectic matrix for Y gate (identity matrix)
        S = np.identity(2 * n, dtype=int)

        # Sign array: Applying Y gate flips the signs of X and Z operators
        r = np.zeros(2 * n, dtype=int)
        r[2 * qubit], r[2 * qubit + 1] = 1, 1     # Flip sign of X and Z

        # Create Clifford object
        cliff = Clifford(0, 0, n, s=S)
        cliff.r = r

        # Evolve stabilizer state
        self.stab = self.stab.evolve(cliff)

    def z(self, qubit):
        """
        Apply a Z gate (Pauli-Z gate) to the specified qubit.

        INPUTS:
        qubit -- int
            The index of the qubit to which the Z gate is applied.

        OUTPUTS:
        None
        """
        n = self.num_qubits

        # Symplectic matrix for Z gate (identity matrix)
        S = np.identity(2 * n, dtype=int)

        # Sign array: Applying Z gate flips the sign of X operators
        r = np.zeros(2 * n, dtype=int)
        r[2 * qubit] =  1 # Flip sign of X

        # Create Clifford object
        cliff = Clifford(0, 0, n, s=S)
        cliff.r = r

        # Evolve stabilizer state
        self.stab = self.stab.evolve(cliff)

    def cx(self, c, t):
        """
        Apply a Controlled-NOT (CX) gate from the control qubit to the target qubit.

        INPUTS:
        c -- int
            The index of the control qubit.
        t -- int
            The index of the target qubit.

        OUTPUTS:
        None
        """
        n = self.num_qubits

        # Symplectic matrix for CX gate
        S = np.identity(2 * n, dtype=int)

        # CX gate action:
        # X_c -> X_c X_t
        # Z_t -> Z_c Z_t
        # Other operators remain the same

        # Update X column
        S[2*t, 2*c] = 1
        # Update Z column
        S[2*c+1, 2*t+1] = 1

        # Sign array (no change in signs)
        r = np.zeros(2 * n, dtype=int)

        # Create Clifford object
        cliff = Clifford(0, 0, n, s=S)
        cliff.r = r

        # Evolve stabilizer state
        self.stab = self.stab.evolve(cliff)

    def cz(self, c, t):
        """
        Apply a Controlled-Z (CZ) gate between the control and target qubits.

        INPUTS:
        c -- int
            The index of the control qubit.
        t -- int
            The index of the target qubit.

        OUTPUTS:
        None
        """
        
        self.h(t)
        self.cx(c,t)
        self.h(t)

    def apply_clifford(self, clifford, qubits):
        """
        Apply a Clifford operation to a subset of qubits.

        INPUTS:
        clifford -- Clifford
            The Clifford object representing the operation to apply.
        qubits -- list of int
            The indices of the qubits to which the Clifford operation is applied.

        OUTPUTS:
        None
        """
        n = self.num_qubits
        n_sub = len(qubits)

        # Initialize the full symplectic matrix as identity
        S_full = np.identity(2 * n, dtype=int)
        r_full = np.zeros(2 * n, dtype=int)

        # Build index arrays for mapping
        full_indices = []
        for q in sorted(qubits):
            full_indices.extend([2 * q, 2 * q + 1])

        # Indices for the small Clifford
        small_indices = list(range(2 * n_sub))

        # Use np.ix_ to select the subsets of the symplectic matrix
        ix_full = np.ix_(full_indices, full_indices)

        # Map the small symplectic matrix into the full symplectic matrix
        S_full[ix_full] = clifford.s

        # Map the small sign array into the full sign array
        r_full[full_indices] = clifford.r

        # Create a new Clifford object for the full system
        full_clifford = Clifford(0, 0, n, s=S_full)
        full_clifford.r = r_full

        # Evolve the stabilizer state
        self.stab = self.stab.evolve(full_clifford)

    def compile(self):
        """
        Convert the CliffordCircuit into a single Clifford object representing the cumulative effect.

        INPUTS:
        None

        OUTPUTS:
        Clifford
            A Clifford object representing the cumulative effect of the circuit.
        """
        n = self.num_qubits
        s_columns = []
        r_bits = []

        # For each qubit, extract the destabilizer (X mapping) and stabilizer (Z mapping)
        for i in range(n):
            # Destabilizer corresponds to how X_i transforms
            destab = self.stab.destabilizers[i]
            destab_vec = destab.smolin_vec  # Symplectic vector of length 2n
            destab_exp = destab.exp % 4

            # Stabilizer corresponds to how Z_i transforms
            stab = self.stab.generators[i]
            stab_vec = stab.smolin_vec  # Symplectic vector of length 2n
            stab_exp = stab.exp % 4

            # Compute sign bits: 0 for positive (+1 or +i), 1 for negative (-1 or -i)
            destab_sign_bit = (destab_exp // 2) % 2  # 1 if exponent corresponds to negative sign
            stab_sign_bit = (stab_exp // 2) % 2      # 1 if exponent corresponds to negative sign

            # Append the symplectic vectors to the columns list
            s_columns.append(destab_vec)
            s_columns.append(stab_vec)

            # Append the sign bits to the r_bits list
            r_bits.append(destab_sign_bit)
            r_bits.append(stab_sign_bit)

        # Construct the symplectic matrix 's' by stacking the columns
        s = np.column_stack(s_columns)  # Shape: (2n, 2n)

        # Convert the list of sign bits to a NumPy array
        r = np.array(r_bits, dtype=int)

        # Create the Clifford object with the compiled symplectic matrix and sign array
        # We can set 'cliff_idx' and 'sign_idx' to 0 since we're providing 's' and 'r' directly
        compiled_clifford = Clifford(cliff_idx=0, sign_idx=0, n=n, s=s)
        compiled_clifford.r = r

        return compiled_clifford