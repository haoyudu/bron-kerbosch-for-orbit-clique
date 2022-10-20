# We use the packages numpy, random, time, and galois.
# numpy lets us use basic math commands in python.
# random lets us pick random items in a set.
# time just lets us see how long things take.
# galois lets us do arithmetic over finite fields.
import numpy as np
import random
import time
import galois

# The following function checks if an integer is prime
from sympy import isprime

def LineSet(q,x,y):
    """
    Summary: Define list of lines on the Fermat surface
    
    Parameters:
    q (int): determined by user's choice of characteristic in setup
    x (galois field element): (q+1)-st root of -1
    y (galois field element): (q^-1)-st root of -1
    
    Returns:
    list of lines of the form [FirstLinearForm,SecondLinearForm], in which each
    linear form is defined by four coefficients
    
    """
    Lines = []
    for i in range(0,q+1):
        for j in range(0,q+1):
            Lines.insert(i*(q+1)+j,[[y**0,x**(2*i+1),x-x,x-x],[x-x,x-x,y**0,x**(2*j+1)],[-1,i,j]]);
            Lines.insert(2*(i*(q+1)+j)+1,[[y**0,x-x,x**(2*i+1),x-x],[x-x,y**0,x-x,x**(2*j+1)],[-2,i,j]]);
            Lines.insert(3*(i*(q+1)+j)+2,[[y**0,x-x,x-x,x**(2*i+1)],[x-x,y**0,x**(2*j+1),x-x],[-3,i,j]]);
    if q>2:
        for i in range(0,q**2-1):
            if not((y**i)**(q+1)== -y**0):
                z = [];
                for j in range(0,q**2):
                    if -y**0-(y**i)**(q+1) == (y**j)**(q+1):
                        z.append(y**j);
                Z=[];
                for j in range(0,len(z)):
                    test = True;
                    for k in range(0,j):
                        if z[k] == z[j]:
                            test = False;
                    if test:
                        Z.append(z[j]);
                for j in range(0,q+1):
                    for k in range(0,q+1):
                        Lines.append([[-y**0,y**i,x-x,Z[j]],[-(y**i)**q,-y**0,Z[k],x-x],[i,j,k]]);
    return Lines


def setup():
    """
    Summary: Generates the base field f of order q^2 where q = p^e and
    p is the characteristic of the field, e is a prime power
    
    p,e are user inputs
    
    Parameters:
    None
    
    Returns:
    None
    
    Modifies:
    global variables p, e, q, power, F, y, x, Lines
    - F is the galois field of characteristic p and order q^2
    - Lines is the set of lines on F
    
    Outputs:
    Print number of lines on Fermat surface
    
    """
    global p, e, q, power, F, y, x, Lines
    print("What characteristic are you working in?")
    p = input()
    while not(p.isdigit()) or not(isprime(int(p))):
        print("Please enter a prime number.")
        p = input()
    print("Which power of it do you wish to do use?")
    e = input()
    while not(e.isdigit()) or not(int(e)>0):
        print("Please input a positive integer.")
        e = input()
    p = int(p)
    e = int(e)
    q = p**e
    # F is the field F_{q^2} & y is a generator of it,
    # i.e. a primitive (q^2-1)-st root of unity
    # uses package galois
    power = e*2;
    F = galois.GF(p**power);
    y = F.primitive_element;
    # x is a primitive (q+1)-st root of -1
    # if p is odd that means it is a primitive 2(q+1)-st root of 1
    if not(p==2):
        x = y**int((q-1)/2);
    else:
        x = y**(q-1);
    # outputs the set of lines and at least checks that it is the right size.
    Lines = LineSet(q,x,y);
    if len(Lines) == (q**3+1)*(q+1):
        print("There are "+str(len(Lines))+" lines on the smooth extremal surface of degree "+str(q+1)+"!");
    else:
        print("Error! The correct number of lines is ",(q**3+1)*(q+1),", but the program only constructed ",len(Lines));


def areSkew_Fq2(i,j):
    """
    Summary: Determines if two lines are skew in field F
    
    Parameters:
    i,j (list): lines of the form [FirstLinearForm,SecondLinearForm]
    
    Returns:
    True if i,j are skew in field F; False otherwise
    
    """
    
    # Create matrix (2d field arrays) of coefficients of two lines
    # First two rows are L1, Last two rows are L2
    L1 = F([Lines[i][0],Lines[i][1]])
    L2 = F([Lines[j][0],Lines[j][1]])
    Mat = np.append(L1,L2,axis=0)
    zero = F(0)
    # If the determinant is nonzero, then we would have a nontrivial null space,
    # which means that there are common zeros between the four polynomials defining both lines
    if(np.linalg.det(Mat) == zero):
        return False
    else:
        return True

def listOfSkewMatrix(N):
    """
    Summary: Returns a matrix with 1's to indicate pairwise skewness of lines i,j and 0 otherwise
    
    Parameters: 
    arg1 (int): Number of lines on Fermat surface over given field F
    
    Returns: 
    2D np array: Adjacency Matrix with 1 in (i,j)th position if line i and line j are skew
    
    """
    skew_mat = np.zeros((N,N))
    for i in range(0,N):
        for j in range(i+1,N):
            if areSkew_Fq2(i,j):
                skew_mat[i][j] = 1
                skew_mat[j][i] = 1
            else:
                skew_mat[i][j] = 0
                skew_mat[j][i] = 0
                
    
    return skew_mat

        
def CommonNeighbours(vertices,graph):
    """
    Summary: Function that takes in line numbers and returns the line numbers of lines that are skew to 'all' of the input lines.
    i.e. returns lines that are skew to all of the input lines.
    
    Parameters:
    arg1 (list): Vertex numbers (Line numbers in the case of lines on surface)
    arg2 (2D np array): Adjacency matrix of graph (0's on leading diagonal)
    
    Returns:
    list: list of vertex(line) numbers adjacent(skew) to all the vertices(lines) in arg1
    
    """

    l = [x for x in Neighbours[vertices[0]]]
    for i in range(1,len(vertices)):
        for j in l[:]:
            if graph[vertices[i]][j] == 0:
                l.remove(j)
               
    return l  

def rowreduce(i,j,k):
    """
    Summary: Using the row_reduce function in numpy to row reduce coefficient matrix on field f.
    Here, the inputs are the a in the linear combination phi(x+a^(2i+1)z) = A(x+a^(2i+1)z)+B(y+z^(2i+1)w) where phi is the stabilizer

    Parameters: 
    i,j,k (int): line m_0 sent to m_i, line m_1 sent to m_j, line m_2 sent to m_j
    
    Returns: 
    2D np field array: Row reduced form of the matrix 
    
    """
    
    matrix = F([[1, x, -((x)**((2*i) + 1)), -((x)**((2*i) + 2))],
                [1, x**3, -((x)**((2*j) + 1)), -((x)**((2*j) + 4))],
                [1, x**5, -((x)**((2*k) + 1)), -((x)**((2*k) + 6))]])
    
    return matrix.row_reduce()


def nullspace(i,j,k):
    """
    Summary: Using the null_space function in numpy to row reduce coefficient matrix (described in summary of rowreduce()) on field f

    Parameters: 
    i,j,k (int): line m_0 sent to m_i, line m_1 sent to m_j, line m_2 sent to m_j
    
    Returns: 
    2D np field array: null space of the matrix discussed in class
    
    """
    matrix = F([[1, x, -((x)**((2*i) + 1)), -((x)**((2*i) + 2))],
                [1, x**3, -((x)**((2*j) + 1)), -((x)**((2*j) + 4))],
                [1, x**5, -((x)**((2*k) + 1)), -((x)**((2*k) + 6))]])
    
    return matrix.null_space()


def stab_auto(i,j,k):
    """
    Summary: Enumerating the stabilizer

    Parameters: 
    i,j,k (int): line m_0 sent to m_i, line m_1 sent to m_j, line m_2 sent to m_j
    
    Returns: 
    2D np field array:  returns  4x4 matrix that stabilizes l_0 l_1 and l_2 and 
                        takes line m_0 to m_i, line m_1 to m_j and line m_2 to m_j
    
    """
    vec = nullspace(i,j,k)
    
    return F([[vec[0][2], 0, vec[0][0], 0],
              [0, vec[0][2], 0, vec[0][0]],
              [vec[0][3], 0, vec[0][1], 0],
              [0, vec[0][3], 0, vec[0][1]]])
    

def automorphism(i, n, matrix):
    """
    Summary: Determining where a line is sent to under automorphism (matrix)

    Parameters: 
    i (int): line number (0 <= i < n) we are applying the automorphism to
    n (int): Total number of lines 
    matrix (2D np field array): Automorphism being applied
    
    Returns: 
    int: returns line number that line i is sent to
    
    """
    matrix = matrix.transpose() #First transpose the automorphism matrix, so we can
                                # apply it straight to the coefficients defining the line
        
    L1 = F([Lines[i][0],Lines[i][1]]).transpose() #to get coefficients of both polys as cols
    
    updatedL1 = np.matmul(matrix, L1).transpose() #2 rows of coefficients of new polys
    
    #Now we check with every line
    
    for j in range(0,n):
        L2 = F([Lines[j][0],Lines[j][1]])
        Mat = np.append(updatedL1,L2,axis=0)
        
        if(np.linalg.matrix_rank(Mat)==2):
            return j
        
    return "No corresponding line found"



def BuildStab():
    """
    Summary: Building the stabilizer of the standard (ordered) three lines.
    
    Parameters: none
    
    Returns: generator function
    
    """
    for i in range(q+1):
        for j in range(q+1):
            if j != i:
                for k in range(q+1):
                    if k != i and k != j:
                        yield stab_auto(i,j,k)


def Stabilizer(NewLine,SetOfMatrices):
    """
    Summary: Finding the matrices from a set which stabilize a particular line.
    
    Parameters:
    NewLine (int): The index of the line we want to stabilize
    SetOfMatrices (list of indices of matrice): This is the set of indices of the matrices we test for stabilizes the line.
    
    Returns: generator function
    
    """
    for matrix in SetOfMatrices:
        if LineAuto[matrix][NewLine] == NewLine:
            yield matrix

def BuildAuto():
    """
    Summary: Building the automorphisms in the stabilizer.
    
    Parameters: none
    
    Returns:
    list: matrix (automorphism)
    
    """
    LineAuto  = [];
    for matrix in range(len(ThreeLineStab)):
        matrixRow = [];
        for line in range(len(Lines)):
            matrixRow.append(automorphism(line, len(Lines), ThreeLineStab[matrix]))
        LineAuto.append(matrixRow)
    return LineAuto



def FindingAnOrbit(Line,SetOfMatrices):
    """
    Summary: This function finds the orbit of a line under the action of a set of matrices
   
    Parameters:
    Line (int): The index of the line we are considering
    SetOfMatrices (list of int): This is the set of indices of matrices acting on the line.
    
    Returns: generator function
    
    """
    for matrix in SetOfMatrices:
        yield LineAuto[matrix][Line]

def FindingOrbitRepresentatives(SetOfLines,SetOfMatrices):
    """
    Summary: This function yields one representative of each orbit of the action of a set of matrices on a set of lines.
    
    Parameters:
    SetOfLines (list of int): This is the set of lines which we want orbit representatives from.
    SetOfMatrices (list of int): This is the indices of the set of matrices acting on the lines.
    
    Returns: generator function
    """
    while len(SetOfLines)>0:
        yield SetOfLines[0]
        Orbit = list(FindingAnOrbit(SetOfLines[0],SetOfMatrices))
        for Line in Orbit:
            if Line in SetOfLines:
                SetOfLines.remove(Line)

def BkMaximalCliquesPivot2Orbit(r,p,x,graph,rStab,isTrivial):
    """
    Summary: Modification of BkMaximalCliquesPivot2.
    If action of the stabilizer of the skew set at each step is not trivial to
    the previous step, finds one representative of the orbits of p under the
    action of Stab(r) and reduces number of calculations by using only that
    representative
    
    Parameters:
    arg1 (list): R (from algorithm defn)
    arg2 (list): P (from algorithm defn)
    arg3 (list): X (from algorithm defn)
    arg4 (2D np array): Adjacency matrix of graph (0's on leading diagonal)
    
    Returns:
    generator object
    
    """
    p_copy = list(p)
    if len(p) == 0 and len(x) == 0:
        yield r
    else:
        for l in Neighbours[(p_copy+x)[0]]:
            if l in p_copy:
                p_copy.remove(l)
        if not(isTrivial):
            OrbitRepresentatives = list(set(FindingOrbitRepresentatives(p_copy[:],rStab)))
            if OrbitRepresentatives == p_copy:
                isTrivial = True
        else:
            OrbitRepresentatives = p_copy;
        for vertex in OrbitRepresentatives:
            r_new = r + [vertex]
            p_new = [val for val in p if val in Neighbours[vertex]] # p intersects N(vertex)
            x_new = [val for val in x if val in Neighbours[vertex]] # x intersects N(vertex)
            for j in BkMaximalCliquesPivot2Orbit(r_new,p_new,x_new,graph,list(Stabilizer(vertex,rStab)),isTrivial):
                yield j
            p.remove(vertex)
            x.append(vertex)


def ListCliqueCounter(CliqueList):
    """
    Summary: Counts number of maximal skew sets of a certain size
    
    Parameters:
    CliqueList (list): list of cliques
    
    Returns:
    None
    
    """
    NumberOfMaximalCliquesOfEachSize = []
    for i in range(len(Lines)):
        Counter = 0
        for clique in CliqueList:
            if len(clique) == i:
                Counter += 1;
        NumberOfMaximalCliquesOfEachSize.append(Counter)
    return NumberOfMaximalCliquesOfEachSize

def CliqueCounter(BaseClique,NumberOfMaximalCliquesOfEachSize):
	"""
    Summary: Printing out how many orbits of maximal skew sets of a certain size containing a base clique
    
    Parameters:
    BaseClique (list): clique to check for
    NumberOfMaximalCliquesOfEachSize (list): output list of ListCliqueCounter
    
    Returns:
    None
    
    """
    for i in range(len(NumberOfMaximalCliquesOfEachSize)):
        if not(NumberOfMaximalCliquesOfEachSize[i] ==0):
            print("There are",NumberOfMaximalCliquesOfEachSize[i],"orbits of maximal skew sets of size",i,"containing the set "+str(BaseClique)+" up to the action of the stablizer of that set.")



def BKCliquesFixedSizeOrbit(k,r,p,x,graph,rStab,isTrivial):
    """
    Summary: Outputs all relevant skew sets of size k in the orbit method.
    
    Parameters:
    arg1 (int): k (relevant length)
    arg2 (list): R (from algorithm defn)
    arg3 (list): P (from algorithm defn)
    arg4 (list): X (from algorithm defn)
    arg5 (2D np array): Adjacency matrix of graph (0's on leading diagonal)
    
    Returns:
    generator object
    
    """
    p_copy = list(p)
    if len(r) == k or (len(p) ==0 and len(x)==0):
        yield [r,p,x,rStab,isTrivial]
            
    else:
        for l in Neighbours[(p_copy+x)[0]]:
            if l in p_copy:
                try:
                    p_copy.remove(l)
                except:
                    pass
                
        if not(isTrivial):
            OrbitRepresentatives = list(set(FindingOrbitRepresentatives(p_copy[:],rStab)))
            if OrbitRepresentatives == p_copy:
                isTrivial = True
        else:
            OrbitRepresentatives = p_copy;
            
        
        for vertex in OrbitRepresentatives:
            r_new = r + [vertex]
            p_new = [val for val in p if val in Neighbours[vertex]] # p intersects N(vertex)
            x_new = [val for val in x if val in Neighbours[vertex]] # x intersects N(vertex)
            for j in BKCliquesFixedSizeOrbit(k,r_new,p_new,x_new,graph, 
                                                    list(Stabilizer(vertex,rStab)),isTrivial):
                yield j
                
            p.remove(vertex)
            x.append(vertex)



def BuildNeighbours():
    """
    Summary: Building the set of neighbors for the overall graph to reduce run time, determining the neighbors of every vertex
    Uses global variable graph
	Credit: https://stackoverflow.com/questions/13904636/implementing-bron-kerbosch-algorithm-in-python for a starting point.

    Returns:
    list of lists: in row i this is a list of vertex numbers adjacent to i
    
    """
    for vertex in range(len(Lines)):
        c = 0
        l = []
        for i in graph[vertex]:
            if i == 1 :
                l.append(c)
            c+=1
        yield l


def main():
    # create the field and the set of lines
    global graph, ThreeLineStab, LineAuto, Neighbours
    setup()
    if q==2 or q==3 or q==4:
        graph = listOfSkewMatrix(len(Lines))
        ThreeLineStab = list(BuildStab())
        LineAuto = list(BuildAuto())
        Neighbours = list(BuildNeighbours())
        r = [0,q+2,2*q+4]
        p = CommonNeighbours(r,graph)
        rStab = range(len(ThreeLineStab))
        isTrivial = False
    if q==2 or q==3:
        CliqueCounter(r,ListCliqueCounter(list(BkMaximalCliquesPivot2Orbit(r,p,[],graph,rStab,isTrivial))))
    elif q==4:
        k=8;
        CliqueSizeList = [];
        for index in range(len(Lines)):
            CliqueSizeList.append(0)
        OrbitList = list(BKCliquesFixedSizeOrbit(k,r,p,[],graph,range(len(ThreeLineStab)),False))
        for index in range(16):
            SkewSetList = [];
            for j in range(5000*index,min(5000*(index+1),len(OrbitList))):
                TempOrbits = list(BkMaximalCliquesPivot2Orbit(OrbitList[j][0],OrbitList[j][1],OrbitList[j][2],graph,OrbitList[j][3],OrbitList[j][4]))
                for SkewSet in TempOrbits:
                    SkewSetList.append(SkewSet);
            NewCliquesSizes = ListCliqueCounter(SkewSetList)
            for i in range(len(Lines)):
                CliqueSizeList[i] = CliqueSizeList[i]+NewCliquesSizes[i]
        CliqueCounter(r,CliqueSizeList)
    else:
        print("This program currently can only run without modification for q=2,3,4. Please input one of those integers or modify the code.")
        main()
    
    
if __name__ == "__main__":
    main()
