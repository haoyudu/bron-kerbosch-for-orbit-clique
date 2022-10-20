# We use the packages numpy, random, time, and galois.
# numpy lets us use basic math commands in python.
# random lets us pick random items in a set.
# time just lets us see how long things take.
# galois lets us do arithmetic over finite fields.
import numpy as np
import random
import time
import galois
import graph4
import LineAuto4
import multiprocessing

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


def listOfSkew(N):
    """
    Summary: Produces list of indices of skew lines in F
    
    Parameters:
    N (int): number of lines in field
    
    Returns:
    list of pairs of indices (i,j) of the list Lines
    
    """
    skew_mat = [list(range(1 + N * i, 1 + N * (i + 1)))
                            for i in range(N)]
    for i in range(0,N):
        for j in range(0,N):
            if areSkew_Fq2(i,j):
                skew_mat[i][j] = 1
            else:
                skew_mat[i][j] = 0
    skew_mat_triu = np.triu(skew_mat)
    list_of_pairs = []
    for i in range(0,N):
        for j in range(0,N):
            if skew_mat_triu[i][j] == 1:
                list_of_pairs.append((i,j))
    return list_of_pairs


def listOfSkewMatrix(N):
    """
    Summary: Variation of the listOfSkew function that returns a matrix
    with 1's to indicate skew and 0 to indicate not skew
    
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



#the Bron-Kerbosch recursive algorithm
def BkMaximalCliques(r,p,x,graph):
    """
    Summary: Uses the Bron-Kerbosch without pivoting
    
    Parameters:
    arg1 (list): R (from algorithm defn)
    arg2 (list): P (from algorithm defn)
    arg3 (list): X (from algorithm defn)
    arg4 (2D np array): Adjacency matrix of graph (0's on leading diagonal)
    
    Returns:
    generator object
    
    """
    if len(p) == 0 and len(x) == 0:
        yield r
    for vertex in p[:]:
        r_new = r + [vertex]
        p_new = [val for val in p if val in Neighbours[vertex]] # p intersects N(vertex)
        x_new = [val for val in x if val in Neighbours[vertex]] # x intersects N(vertex)
        for j in  BkMaximalCliques(r_new,p_new,x_new,graph):
            yield j
            
        p.remove(vertex)
        x.append(vertex)

def BkMaximalCliquesPivot(r,p,x,graph):
    """
    Summary: Uses the algorithm with pivoting.
    
    Parameters:
    arg1 (list): R (from algorithm defn)
    arg2 (list): P (from algorithm defn)
    arg3 (list): X (from algorithm defn)
    arg4 (2D np array): Adjacency matrix of graph (0's on leading diagonal)
    
    Returns:
    generator object
    
    """
    if len(p) == 0 and len(x) == 0:
        yield r
    
    else:
        u = Neighbours[(p+x)[0]]
        
    for vertex in p[:]:
        flag = True
        for k in u:
            if k == vertex:
                flag = False
        
        if flag:
            r_new = r + [vertex]
            p_new = [val for val in p if val in Neighbours[vertex]] # p intersects N(vertex)
            x_new = [val for val in x if val in Neighbours[vertex]] # x intersects N(vertex)
            for j in BkMaximalCliquesPivot(r_new,p_new,x_new,graph):
                yield j
            
            p.remove(vertex)
            x.append(vertex)
            
def BkMaximalCliquesPivot2(r,p,x,graph):
    """
    Summary: Uses the algorithm with pivoting - quicker version than above function. Halves the time.
    
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
            try:
                p_copy.remove(l)
            except:
                pass
        
    for vertex in p_copy:
        
        r_new = r + [vertex]
        p_new = [val for val in p if val in Neighbours[vertex]] # p intersects N(vertex)
        x_new = [val for val in x if val in Neighbours[vertex]] # x intersects N(vertex)
        for j in BkMaximalCliquesPivot2(r_new,p_new,x_new,graph):
            yield j
            
        p.remove(vertex)
        x.append(vertex)


def BkMaximalCliquesPivot2Test(r,p,x,graph,n):
    """
    Summary: Uses the algorithm with pivoting with additional parameter n to only consider the n x n submatrix of 
    the graph matrix. This function is essentially to test the speed of the above function for matrices that the algorithm is
    taking too long for, and in order to just apply the algorithm on a submatrix to see how long it is taking.
    
    Parameters:
    arg1 (list): R (from algorithm defn)
    arg2 (list): P (from algorithm defn)
    arg3 (list): X (from algorithm defn)
    arg4 (2D np array): Adjacency matrix of graph (0's on leading diagonal)
    arg5 (int): the size of the submatrix you want to run the algorithm on
    
    Returns:
    generator object
    
    """
    p_copy = list(p)
    graph = [graph[i][0:n] for i in range(0,n)]
    
    if len(p) == 0 and len(x) == 0:
        yield r
    
    else:
        for l in Neighbours[(p_copy+x)[0]]:
            p_copy.remove(l)
    
    for vertex in p_copy:
        
        r_new = r + [vertex]
        p_new = [val for val in p if val in Neighbours[vertex]] # p intersects N(vertex)
        x_new = [val for val in x if val in Neighbours[vertex]] # x intersects N(vertex)
        for j in BkMaximalCliquesPivot2Test(r_new,p_new,x_new,graph):
            yield j
            
        p.remove(vertex)
        x.append(vertex)
        
def CommonNeighbours(vertices,graph):
    """
    Summary: Function that takes in line numbers and returns the line numbers of lines that are skew to 'all' of the input lines.
    i.e. returns lines that are skew to all of the input lines.
    
    Parameters:
    arg1 (list): Vertex Numbers (Line number in the case of lines on surface)
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

def AllNeighbours(vertices, graph):
    """
    Summary: Function that takes in line numbers and returns the line numbers of lines that are skew to 'any' of the input lines.
    i.e. returns lines that skew to atleast one of the input lines.
   
    Parameters:
    vertices (list): Vertex Numbers (Line number in the case of lines on surface)
    graph (2D np array): Adjacency matrix of graph (0's on leading diagonal)
    
    Returns:
    list: list of vertex(line) numbers adjacent(skew) to atleast one of the vertices(lines) in arg1
    
    """
    l = []
    for i in range(0,len(graph)):
        for j in range(0, len(vertices)):
            if graph[vertices[j]][i] == 1:
                l.append(i)
                break
    
    return l
 

def CliqueCounter(BaseClique,CliqueList):
    """
    Summary: Counts number of maximal skew sets of a certain size that contains
    a base clique
    
    Parameters:
    BaseClique (list): clique to check for
    CliqueList (list): list of cliques
    
    Returns:
    None
    
    """
    for i in range(len(Lines)):
        Counter = 0
        for clique in CliqueList:
            if len(clique) == i:
                Counter += 1;
        if not (Counter ==0):
            print("There are",Counter,"maximal skew sets of size",i,"containing the set "+str(BaseClique)+".")


def rowreduce(i,j,k):
    """
    Parameters: 
    i,j,k (int): line m_0 sent to m_i, line m_1 sent to m_j, line m_2 sent to m_j
    
    Returns: 
    2D np field array: Row reduced form of the matrix discussed in class
    
    """
    
    matrix = F([[1, x, -((x)**((2*i) + 1)), -((x)**((2*i) + 2))],
                [1, x**3, -((x)**((2*j) + 1)), -((x)**((2*j) + 4))],
                [1, x**5, -((x)**((2*k) + 1)), -((x)**((2*k) + 6))]])
    
    return matrix.row_reduce()


def nullspace(i,j,k):
    """
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
    Summary: This function builds the stabilizer of the standard three lines.
    
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
    Summary: This function finds the matrices from a set which stabilize a particular line.
    
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
    Summary: This function builds the automorphisms in the stabilizer.
    
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
    global counter
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

def BkMaximalCliquesFixedSizeOrbit(k,r,p,x,graph,rStab,isTrivial):
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
    if len(r) == k or (len(p) == 0 and len(x)==0):
        if (len(p) ==0 and len(x)==0):
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
            for j in BkMaximalCliquesFixedSizeOrbit(k,r_new,p_new,x_new,graph,
                                                   list(Stabilizer(vertex,rStab)),isTrivial):
                yield j
                
            p.remove(vertex)
            x.append(vertex)



def BKPool(orbit):
    #for orbit in Orbits:
    for SkewSet in BkMaximalCliquesPivot2Orbit(orbit[0],orbit[1],orbit[2],graph,orbit[3],orbit[4]):
        yield SkewSet


def BKHybrid(k,r,p,x,graph,rStab,isTrivial):
    if k>3:
        Orbits = BKCliquesFixedSizeOrbit(k,r,p,x,graph,rStab,isTrivial)
    else:
        for SkewSet in BkMaximalCliquesPivot2(r,p,x,graph):
            yield SkewSet

#Graph theory functions:
#Credit: https://stackoverflow.com/questions/13904636/implementing-bron-kerbosch-algorithm-in-python for a starting point.
#Function determines the neighbors of every vertex
def BuildNeighbours():
    """
    Uses global variable graph

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


def RandomSkewSet(Option,Print="False"):
    while not(Option=="minimized") and not(Option=="maximized") and not(Option=="random"):
        print("Please input one of the three valid options: minimized, maximized, or random.")
        Option = input();
    if Option == "random":
        start_time = time.time()
        SkewSet = [0,q+2,2*q+4];
        RemainingPossibleLines = CommonNeighbours(SkewSet,graph);
        while len(RemainingPossibleLines)>0:
            counter = len(RemainingPossibleLines)
            Index = random.randint(0,len(RemainingPossibleLines)-1);
            SkewSet.append(RemainingPossibleLines[Index]);
            RemainingPossibleLines.remove(RemainingPossibleLines[Index])
            RemainingPossibleLines = [x for x in RemainingPossibleLines if x in Neighbours[SkewSet[-1]]];
    elif Option == "maximized" or Option == "minimized":
        start_time = time.time()
        Index =random.randint(0,len(Lines)-1);
        SkewSet = [0,q+2,2*q+4];
        RemainingPossibleLines = CommonNeighbours(SkewSet,graph);
        '''
        for i in range(2):
            print(len(SkewSet),len(RemainingPossibleLines),round(time.time() - start_time,2));
            start_time = time.time();
            Index =random.randint(0,len(RemainingPossibleLines)-1);
            SkewSet.append(RemainingPossibleLines[Index]);
            RemainingPossibleLines.remove(RemainingPossibleLines[Index])
            NewRemainingPossibleLines = [line for line in RemainingPossibleLines if line in Neighbours(Index,graph)];
            RemainingPossibleLines = NewRemainingPossibleLines;
        '''
        while len(RemainingPossibleLines)>0:
            counter = len(RemainingPossibleLines)
            #print(len(SkewSet),len(RemainingPossibleLines),round(time.time() - start_time,2));
            start_time = time.time();
            Chance = random.randint(0,5);
            if Chance<4:
                Min = len(RemainingPossibleLines);
                Max = 0;
                PossibleNextSkewLines = [];
                for i in range(0,len(RemainingPossibleLines)):
                    Temp = len([line for line in RemainingPossibleLines if line in Neighbours[RemainingPossibleLines[i]]]);
                    if (Temp < Min and Option =="minimized") or (Temp > Max and Option =="maximized"):
                        Min = Temp;
                        PossibleNextSkewLines = [RemainingPossibleLines[i]];
                    elif (Temp == Min and Option =="minimized") or (Temp == Max and Option =="maximized"):
                        PossibleNextSkewLines.append(RemainingPossibleLines[i]);
                Index = random.randint(0,len(PossibleNextSkewLines)-1);
                SkewSet.append(PossibleNextSkewLines[Index]);
                RemainingPossibleLines.remove(PossibleNextSkewLines[Index])
                RemainingPossibleLines = [line for line in RemainingPossibleLines if line in Neighbours[PossibleNextSkewLines[Index]]];
            else:
                counter = len(RemainingPossibleLines)
                Index = random.randint(0,len(RemainingPossibleLines)-1);
                SkewSet.append(RemainingPossibleLines[Index]);
                RemainingPossibleLines.remove(RemainingPossibleLines[Index])
                RemainingPossibleLines = [x for x in RemainingPossibleLines if x in Neighbours[SkewSet[-1]]];        
    if q ==4:
        if len(SkewSet)<13  or len(SkewSet)>25:
            print(len(SkewSet),SkewSet)
    elif Print == True or len(SkewSet) >(q+1)**2 or len(SkewSet)<q**2+1 or (p==2 and len(SkewSet) >q**2+q):
        print(SkewSet)
    return len(SkewSet);

# This finds a fixed number of random skew sets and outputs how many of each size occur.
# number is a positive integer
def RandomSkewSetLengths(number,Option="random"):
    start_time = time.time()
    Lengths = [];
    NumberOfLengths = 0;
    for i in range(0,number):
        Lengths.append(RandomSkewSet(Option));
        if i%10000 == 9999:
            print(str(i+1)+" of "+str(number)+" has taken "+str(round(time.time() - start_time,2))+" seconds.")
    for k in range(1,max(Lengths)+1):
        x = 0;
        for i in range(0,number):
            if Lengths[i]==k:
                x = x+1;
        if not(x==0):
            NumberOfLengths = NumberOfLengths+x;
            print("There are "+str(x)+" skew sets of size "+str(k)+" in this collection.")
    if not(NumberOfLengths == number):
        print("An error occured");

def StarsOnALine(Line):
    StarList = [];
    IntersectingLines = [x for x in range(len(Lines)) if not(x in Neighbours[Line])]
    IntersectingLines.remove(Line)
    while len(IntersectingLines)>0:
        IntersectingLine = IntersectingLines[0]
        Star = [x for x in IntersectingLines if not(x in Neighbours[IntersectingLine])]
        Star.insert(0,Line)
        StarList.append(Star)
        for NewLine in Star:
            if not(NewLine ==Line):
                IntersectingLines.remove(NewLine)
    return StarList



def GhostLinePairs(ListOfStars,LineList):
    ListOfGhostLinePairs = []
    while len(ListOfStars)>0:
        Star = ListOfStars[0]
        GhostLinePair = [Star]
        ListOfStars.remove(Star)
        Test = True;
        while Test:
            NewStar = GhostLinePair[-1]
            Counter = len(ListOfStars);
            for Line in [x for x in NewStar if not(x in LineList)]:
                for TestStar in ListOfStars:
                    if Line in TestStar:
                        GhostLinePair.append(TestStar)
                        ListOfStars.remove(TestStar)
            if len(ListOfStars) == Counter:
                Test = False
        ListOfGhostLinePairs.append(GhostLinePair)
    return ListOfGhostLinePairs

'''
for i in range(1):
    LineList = [25,31,37,43,49];
    StarList = [];
    for Line in LineList:
        for Star in StarsOnALine(Line):
            StarList.append(Star)
    GhostLinePairList = GhostLinePairs(StarList,LineList)
    for GhostLinePair in GhostLinePairList:
        print(GhostLinePair)
'''


def SkewSetTest(skew_set, graph):
    """
    Summary: Checks if an input set in a skew set or not.
    
    Parameters:
    arg1 (list): skew_set (set to be tested)
    arg2 (2D np array): Adjacency matrix of graph (0's on leading diagonal)
    
    Returns:
    Boolean based on if input list is a skew set or not
    """
    for i in range(len(skew_set)):
        for j in range(i+1, len(skew_set)):
            if graph[skew_set[i]][skew_set[j]] == 0:
                return False
    
    return True


def MaximalSkewSetTest(skew_set, graph):
    """
    Summary: Checks if an input set in a maximal skew set or not.
    
    Parameters:
    arg1 (list): skew_set (set to be tested)
    arg2 (2D np array): Adjacency matrix of graph (0's on leading diagonal)
    
    Returns:
    Boolean based on if input list is a maximal skew set or not
    """
    for i in range(len(skew_set)):
        for j in range(i+1, len(skew_set)):
            if graph[skew_set[i]][skew_set[j]] == 0:
                return False
    
    if (len(CommonNeighbours(skew_set, graph)) > 0):
        return False
    else:
        return True

def LargestSkewSet(X):
    indexes = [];
    size = len(X[0]);
    for i in range(len(X)):
        if len(X[i]) > size:
            size = len(X[i]);
    for i in range(len(X)):
        if len(X[i]) == size:
            indexes.append(i)
    return size, indexes
    
def SmallestSkewSet(X):
    indexes = [];
    size = len(X[0]);
    for i in range(len(X)):
        if len(X[i]) < size:
            size = len(X[i]);
            index = i;
    for i in range(len(X)):
        if len(X[i]) == size:
            indexes.append(i)
    return size, indexes

def CountCliquesOfSizeN(X,n):
    counter = 0;
    indexlist = [];
    for i in range(len(X)):
        if len(X[i]) == n:
            counter += 1;
            indexlist.append(i)
    return counter, indexlist



def main():
    # create the field and the set of lines
    global graph, ThreeLineStab, counter, LineAuto, Neighbours
    setup()
    if q==4:
        graph = graph4.graph
    else: 
        graph = listOfSkewMatrix(len(Lines))
    ThreeLineStab = list(BuildStab())
    Neighbours = list(BuildNeighbours())
    if q==4:
        LineAuto = LineAuto4.LineAuto
    else:
        LineAuto = list(BuildAuto())
    r = [0,q+2,2*q+2]
    p = CommonNeighbours(r,graph)
    rStab = range(len(ThreeLineStab))
    print("Do you want to run BkMaximalCliquesPivot2Orbit (y/n)? (Output will be called X)")
    answer = input()
    if answer == "y":
        start = time.time()
        X = BkMaximalCliquesPivot2Orbit(r,p,[],graph,rStab,False)
        print(time.time()-start)
    
    #start = time.time()
    #CliqueCounter([0,4,8],list(BkMaximalCliquesPivot2(r,p,[],graph)))
    #print(time.time()-start)
    
    #start = time.time()
    #counter = 0
    #CliqueCounter([0,q+2,2*q+2],list(BkMaximalCliquesPivot2Orbit(r,p,[],graph,range(len(ThreeLineStab)),False)))
    #print(time.time()-start)
    
#if __name__ == "__main__":
    #setup()
    #k = 8
    #r = [0,q+2,2*q+4]
    #p = CommonNeighbours(r,graph)
    #x = []
    #rStab = range(len(ThreeLineStab))
    #isTrivial =  False
    #start = time.time()
    #counter = 0
    #Orbits = list(BKCliquesFixedSizeOrbit(k,r,p,x,graph,rStab,isTrivial))
    #pool = multiprocessing.Pool(2)
    #pool.map(BKPool,Orbits)
    #pool.close()
    #print(time.time()-start,len(Orbits))
