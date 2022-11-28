"""
This program is useful for studying the set of lines on the Fermat surface S of
degree q+1 over the field with q^2 elements where q = p^e, p is a prime integer,
and e is a positive integer. This is a smooth extremal surface.

The code is roughly broken up into six sections.
The first section constructs the skewness matrix of the lines on the Fermat
surface of degree q+1 over F_q^2 (this is the incidence matrix of the
complement of the incidence graph of the lines).
The second section constructs the matrix of neighbours for each vertex in a
graph.
The third section is the BK algorithm with pivoting and a related function.
The fourth section constructs the stabilizer of three (fixed) lines and the
matrix that records where each automorphism sends each line.
The fifth section is the BK algorithm with pivoting for orbits and the
necessary functions.
The sixth section is the setup and main functions that run the program.

Sections 2 and 3 should work for an arbitrary graph.
Section 5 should work for an arbitrary graph given the matrix that records
where each of its automorphisms sends each vertex.
Sections 1, 4, and 6 are very specific to the case of lines on S.
"""

import numpy as np
import galois
"""
numpy : allows basic math commands in python.
galois : allows arithmetic over finite fields.
"""

from sympy import isprime
"""
Summary: Determines if an integer is prime.

Parameters:
n (int): The integer to test

Returns:
True if n is prime. False if it is not.
"""



#Section 1: Constructing the graph
def LineSet(q,x,y):
    """
    Summary: Define list of lines on the Fermat surface
    
    Parameters:
    q (int): determined by user's choice of characteristic in setup
    x (galois field element): primitive(q+1)-st root of -1
    y (galois field element): primitive(q^2-1)-st root of 1
    
    Returns:
    The list of lines on S given in the form [FirstLinearForm,SecondLinearForm],
    in which each linear form is given by a list of its four coefficients
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
                SolutionList = [];
                for j in range(0,q**2):
                    if -y**0-(y**i)**(q+1) == (y**j)**(q+1):
                        SolutionList.append(y**j);
                SolutionSet=[];
                for j in range(0,len(SolutionList)):
                    isNew = True;
                    for k in range(0,j):
                        if SolutionList[k] == SolutionList[j]:
                            isNew = False;
                    if isNew:
                        SolutionSet.append(SolutionList[j]);
                for j in range(0,q+1):
                    for k in range(0,q+1):
                        Lines.append([[-y**0,y**i,x-x,SolutionSet[j]],[-(y**i)**q,-y**0,SolutionSet[k],x-x],[i,j,k]]);
    return Lines


def AreSkew(i,j):
    """
    Summary: Determines if two lines indexed by i and j are skew over a field F
    by checking the determinant of the matrix of their coefficients.
    
    Parameters:
    i,j (integers): indices of lines of the form
    [FirstLinearForm,SecondLinearForm]
    
    Returns:
    True if i,j are skew in field F; False otherwise.
    """
    Matrix = np.append(F([Lines[i][0],Lines[i][1]]),F([Lines[j][0],Lines[j][1]]),axis=0)
    if(np.linalg.det(Matrix) == F(0)):
        return False
    else:
        return True


def SkewnessMatrix(q):
    """
    Summary: Returns a matrix with 1's to indicate pairwise skewness of lines
    i,j and 0 otherwise
    
    Parameters: 
    q (int): q where the field F has q^2 elements 
    
    Returns: 
    (2D np array): Matrix with 1 in (i,j)th position if line i and line j
    are skew and 0 otherwise
    """
    N = (q**3+1)*(q+1)
    skew_mat = np.zeros((N,N))
    for i in range(0,N):
        for j in range(i+1,N):
            if AreSkew(i,j):
                skew_mat[i][j] = 1
                skew_mat[j][i] = 1
    return skew_mat



#Section 2: Creating the neighbours matrix
def CommonNeighbours(VertexList,graph):
    """
    Summary: This function takes in a list of vertices and returns the indices
    of all vertices that are connected to all of those vertices by an edge.
        In our setting, this function that takes in a list of indices of lines
        and returns the indices of lines that are skew to 'all' of the input
        lines, i.e., returns lines that are skew to all of the input lines.
    
    Parameters:
    VertexList (list): A list of vertices in the graph.
        In our setting, a list of indices of lines on S.
    graph (2D np array): A graph containing the vertices.
        In our setting, the skewness matrix of the lines on S
    
    Returns:
    CommonSkewLines(list): The list of indices which are connected to every
    vertex in VertexList.
        In our setting, the list of line indices skew to all the lines in
        VertexList.
    """

    CommonSkewLines = [x for x in Neighbours[VertexList[0]]]
    for i in range(1,len(VertexList)):
        for j in CommonSkewLines[:]:
            if graph[VertexList[i]][j] == 0:
                CommonSkewLines.remove(j)         
    return CommonSkewLines  


def BuildNeighbours(graph):
    """
    Summary: Building the set of neighbors for each vertex in a graph to
    reduce run time.
    	Credit: This is originally based on the code at
            https://stackoverflow.com/questions/13904636/
                implementing-bron-kerbosch-algorithm-in-python
        In our setting, this builds the list of lines skew to each lines on S.

    Parameters:
    graph(2D np array): The graph for which we are building the set of
    neighbors for each vertex.
        In our setting, this is the skewness graph for the lines on S.

    Returns:
    NeighboursOfVertex (list of lists): row i is the list of neighbours of
    the vertex with index i
        In our setting, row i is the list of lines skew to the line with
        index i.
    """
    for Vertex in range(len(graph)):
        TestVertex = 0
        NeighboursOfVertex = []
        for i in graph[Vertex]:
            if i == 1 :
                NeighboursOfVertex.append(TestVertex)
            TestVertex+=1
        yield NeighboursOfVertex 



#Section 3: BK algorithm with pivoting
def ListSizeCounter(ListOfLists):
    """
    Summary: Counts number of lists of each size in a list of lists.
    
    Parameters:
    ListOfLists (list): list of lists whose size is to be counted
    
    Returns:
    NumberOfListsOfEachSize(list): A list whose whose i-th entry is the number
    of lists of size i bigger than the smallest list in ListOfLists.
    """
    MaxLength = len(ListOfLists[0]);
    MinLength = len(ListOfLists[0]);
    for TestList in ListOfLists:
        if len(TestList)<MinLength:
            MinLength = len(TestList);
        elif len(TestList)>MaxLength:
            MaxLength = len(TestList);
    NumberOfListsOfEachSize = []
    for Length in range(MinLength,MaxLength+1):
        Counter = 0
        for List in ListOfLists:
            if len(List) == Length:
                Counter += 1;
        if not(Counter ==0):
            NumberOfListsOfEachSize.append([Length,Counter])
    return NumberOfListsOfEachSize


def BKWithPivoting(R,P,X,graph):
    """
    Summary: This is an implementation of the Bron-Kerbosch algorithm with
    pivoting.
    
    Parameters:
    R (list): R (from algorithm defn) The vertices already in the clique.
        In our setting, the lines already in your skew set.
    P (list): P (from algorithm defn) The vertices that can be added to the
    clique.
        In our setting, this is the lines that can be added to the skew set.
    X (list): X (from algorithm defn) The vertices that are connected to every
    vertex in the clique but have already been considered.
        In our setting, this is the lines that have already been considered
    graph (2D np array): Adjacency matrix of graph
        In our setting, this is the skewness matrix of the lines on S
    
    Returns:
    (generator object): All maximal cliques in graph.
        In our setting, this is all maximal skew sets on S.
    """
    P_copy = list(P)
    if len(P) == 0 and len(X) == 0:
        yield R
    else:
        for vertex in Neighbours[(P_copy+X)[0]]:
            if vertex in P_copy:
                P_copy.remove(vertex)
        for vertex in P_copy:
            R_new = R + [vertex]
            P_new = [val for val in P if val in Neighbours[vertex]] 
            X_new = [val for val in X if val in Neighbours[vertex]]
            for j in BKWithPivoting(R_new,P_new,X_new,graph):
                yield j
            P.remove(vertex)
            X.append(vertex)



#Section 4: Constructing the action of the automorphisms
def Automorphism(i,j,k):
    """
    Summary: Gives the matrix of the linear transformation that takes the lines
    with indices 0,q+2,2q+4,(q+1)^2,(q+1)^2+q+2,(q+1)^2+2q+4 to the lines with
    indices 0,q+2,2q+4,(q+1)^2+i(q+2),(q+1)^2+j(q+2),(q+1)^2+k(q+2),
    respectively.

    Parameters: 
    i,j,k (int): 0<=i,j,k<q+1, none equal. These determine the output matrix,
    see summary.
    
    Returns: 
    (2D np field array): returns 4x4 matrix of the linear transformation in
    the summary.

    Note: The transformation fixing the lines with indices 0,q+2, and 2q+4
    forces it to be of the form [[a,0,b,0],[0,e,0,f],[c,0,d,0],[0,g,0,h]].
    The line with index (q+1)^2 being sent to a line with index (q+1)^2+i(q+2)
    where 0<=i<q+1 forces e=a,f=b,g=c, and h=d. Writing out what it means to
    take the lines with indices (q+1)^2,(q+1)^2+q+2,(q+1)^2+2q+4 to the lines
    with indices (q+1)^2+i(q+2), (q+1)^2+j(q+2),(q+1)^2+k(q+2) is the same as
    saying that [a,b,c,d]^T is in the kernel of the matrix given below as
    "matrix".
    """
    matrix = F([[1, x, -((x)**((2*i) + 1)), -((x)**((2*i) + 2))],
                [1, x**3, -((x)**((2*j) + 1)), -((x)**((2*j) + 4))],
                [1, x**5, -((x)**((2*k) + 1)), -((x)**((2*k) + 6))]])
    vec = matrix.null_space()
    
    return F([[vec[0][2], 0, vec[0][0], 0],
              [0, vec[0][2], 0, vec[0][0]],
              [vec[0][3], 0, vec[0][1], 0],
              [0, vec[0][3], 0, vec[0][1]]]).transpose()


def ImageLine(i, matrix):
    """
    Summary: Determines where a line is sent by an automorphism
    
    Parameters: 
    i (int): Index of the line we are applying the automorphism to
    matrix (2D np field array): The automorphism being applied
    
    Returns: 
    (int): returns index of the line that line i is sent to by matrix
    """
        
    ImageLine = np.matmul(matrix, F([Lines[i][0],Lines[i][1]]).transpose()).transpose() 
    
    for j in range(0,len(Lines)):
        TestLine = F([Lines[j][0],Lines[j][1]])
        Mat = np.append(ImageLine,TestLine,axis=0)
        
        if(np.linalg.matrix_rank(Mat)==2):
            return j
        
    return "No corresponding line found"


def BuildStab():
    """
    Summary: Builds the ordered stabilizer of the three lines with indices
    0,q+2,2q+4.
    
    Parameters: None
    
    Returns:
    (generator function): Of the matrices which stabilizes lines 0,q+2,2q+4.
    """
    for i in range(q+1):
        for j in range(q+1):
            if j != i:
                for k in range(q+1):
                    if k != i and k != j:
                        yield Automorphism(i,j,k)


def Stabilizer(Line,SetOfMatrices):
    """
    Summary: Finding the matrices from a set which stabilize a particular line.
    
    Parameters:
    Line (int): The index of the line we want to stabilize
    SetOfMatrices (list of indices of matrices): This is the set of indices of the
    matrices we test for stabilizing the line
    
    Returns:
    (generator function): Of the indices of the matrices stabilizing the line
    """
    for matrix in SetOfMatrices:
        if LineAuto[matrix][Line] == Line:
            yield matrix


def BuildAuto():
    """
    Summary: Building the matrix that keeps track of where each automorphism in
    the stabilizer of the ordered list of lines with indices 0,q+2,2q+4 sends
    each line.
    
    Parameters: None

    Returns:
    LineAuto (matrix): LineAuto is matrix whose (i,j)-th entry is the index of
    the line which the i-th matrix in the stabilizer of the three lines with
    indices 0,q+2,2q+4 sends the line with index j to.
    """
    LineAuto  = [];
    for matrix in range(len(ThreeLineStab)):
        matrixRow = [];
        for line in range(len(Lines)):
            matrixRow.append(ImageLine(line, ThreeLineStab[matrix]))
        LineAuto.append(matrixRow)
    return LineAuto



#Section 5:BK algorithm with pivots for orbits
def FindingAnOrbit(Vertex,SetOfAutomorphisms,LineAuto):
    """
    Summary: This function finds the orbit of a vertex under the action of a
    set of automorphisms of the graph those vertices are in.
        In our setting, it finds the orbit of a line under the action of a
    set of matrices.
   
    Parameters:
    Line (int): The index of the vertex we are considering.
        In our setting, this is the index of a line on S
    SetOfMatrices (list of int): This is the set of indices of automorphisms
    acting on the vertex.
        In our setting, this is a set of indices of matrices which are
        automorphisms of S
    LineAuto (matrix): This is the matrix whose (i,j)-th entry is the image of
    the vertex with index j under the image of the automorphism with index i.
        In our setting, this is the matrix whose (i,j)-th entry is the image of
        the line with index j under the image of the matrix with index i.
    
    Returns:
    (generator function): Of the indices of the images of the vertex under each
    automorphism.
        In our setting, it is the generator function of the indices of the
        images of the vertex under each automorphism.
    """
    for aut in SetOfAutomorphisms:
        yield LineAuto[aut][Vertex]


def FindingOrbitRepresentatives(SetOfVertices,SetOfAutomorphisms):
    """
    Summary: In this setting, this function yields one representative of each
    orbit of the action of a set of automorphisms on a set of vertices.
        In this setting, this function yields one representative of each orbit
        of the action of a set of matrices on a set of lines.
    
    Parameters:
    SetOfVertices (list of int): this is the set of vertices which we want orbit
    representatives from.
        In our setting, this is the set of lines which we want orbit
        representatives from.
    SetOfAutomorphisms (list of int): This is the indices of the set of
    automorphisms acting on the vertices.
        In our setting, this is the indices of the set of matrices acting on
        the lines.
    
    Returns:
    (generator function): Of the indices of the orbit representatives.
    """
    while len(SetOfVertices)>0:
        yield SetOfVertices[0]
        Orbit = list(FindingAnOrbit(SetOfVertices[0],SetOfAutomorphisms,LineAuto))
        for Vertex in Orbit:
            if Vertex in SetOfVertices:
                SetOfVertices.remove(Vertex)


def BKWithPivotingOrbit(R,P,X,graph,StabR,isTrivial):
    """
    Summary: Modification of BKWithPivoting to account for orbits.
    If the action of the stabilizer of the skew set at each step is not
    trivial to the previous step, finds one representative of the orbits of P
    under the action of Stab(R) and reduces number of calculations by using
    only that representative.
    
    Parameters:
    R (list): R (from algorithm defn) The vertices already in the clique.
        In our setting, the lines already in your skew set.
    P (list): P (from algorithm defn) The vertices that can be added to the
    clique.
        In our setting, this is the lines that can be added to the skew set.
    X (list): X (from algorithm defn) The vertices that are connected to every
    vertex in the clique but have already been considered.
        In our setting, this is the lines that have already been considered
    graph (2D np array): Adjacency matrix of graph.
        In our setting, this is the skewness matrix of the lines on S.
    StabR (list): The indices of the automorphisms that stabilize R.
    isTrivial (boolean): Records if StabR acts trivially on P.

    Returns:
    (generator object): Representative of all orbits of maximal cliques in
    the graph.
        In our setting, this is representative of all orbits of maximal skew
        sets on S.
    """
    P_copy = list(P)
    if len(P) == 0 and len(X) == 0:
        yield R
    else:
        for vertex in Neighbours[(P_copy+X)[0]]:
            if vertex in P_copy:
                P_copy.remove(vertex)
        if not(isTrivial):
            OrbitRepresentatives = list(set(FindingOrbitRepresentatives(P_copy[:],StabR)))
            if OrbitRepresentatives == P_copy:
                isTrivial = True
        else:
            OrbitRepresentatives = P_copy;
        for vertex in OrbitRepresentatives:
            R_new = R + [vertex]
            P_new = [val for val in P if val in Neighbours[vertex]] 
            X_new = [val for val in X if val in Neighbours[vertex]]
            for j in BKWithPivotingOrbit(R_new,P_new,X_new,graph,list(Stabilizer(vertex,StabR)),isTrivial):
                yield j
            P.remove(vertex)
            X.append(vertex)


def BKCliquesFixedSizeOrbit(Size,R,P,X,graph,StabR,isTrivial):
    """
    Summary: Outputs a representative of every orbit of cliques of size Size.
    
    Parameters:
    Size (int): The size of the cliques we want to list
    R (list): R (from algorithm defn) The vertices already in the clique.
        In our setting, the lines already in your skew set.
    P (list): P (from algorithm defn) The vertices that can be added to the
    clique.
        In our setting, this is the lines that can be added to the skew set.
    X (list): X (from algorithm defn) The vertices that are connected to every
    vertex in the clique but have already been considered.
        In our setting, this is the lines that have already been considered
    graph (2D np array): Adjacency matrix of graph.
        In our setting, this is the skewness matrix of the lines on S
    StabR (list): The indices of the automorphisms that stabilize R.
    isTrivial (boolean): Records if StabR acts trivially on P.
    
    Returns:
    (generator object): Of represenativers of the orbits of the cliques of
    size Size in graph.
        In our setting, this is represenatives of the orbits of skew sets of
        size Size of lines on S.
    """
    P_copy = list(P)
    if len(R) == Size or (len(P) ==0 and len(X)==0):
        yield [R,P,X,StabR,isTrivial]
            
    else:
        for vertex in Neighbours[(P_copy+X)[0]]:
            if vertex in P_copy:
                try:
                    P_copy.remove(vertex)
                except:
                    pass
                
        if not(isTrivial):
            OrbitRepresentatives = list(set(FindingOrbitRepresentatives(P_copy[:],StabR)))
            if OrbitRepresentatives == P_copy:
                isTrivial = True
        else:
            OrbitRepresentatives = P_copy;
            
        for vertex in OrbitRepresentatives:
            R_new = R + [vertex]
            P_new = [val for val in P if val in Neighbours[vertex]]
            X_new = [val for val in X if val in Neighbours[vertex]] 
            for clique in BKCliquesFixedSizeOrbit(Size,R_new,P_new,X_new,graph, 
                                             list(Stabilizer(vertex,StabR)),isTrivial):
                yield clique
                
            P.remove(vertex)
            X.append(vertex)



#Section 6: Main running functions
def setup():
    """
    Summary: Generates the base field f of order q^2 where q = p^e,
    p is the characteristic of the field, and e is any integer power.
    In addition, generates the set of lines on S, their skewness matrix, and
    the matrix listing the lines skew to each line.
    
    User inputs:
    p(int): the characteristic of the field
    e(int): the power determining the size of the field
    
    Parameters: None
    
    Returns: None
    
    Modifies:
    global variables p, e, q, F, y, x, Lines, graph, Neighbours, ThreeLineStab,
    LineAuto
    - F is the galois field of characteristic p and order q^2
    - Lines is the set of lines on F
    - graph is the skewness matrix of the lines in Lines
    - Neighbours is the matrix whose ith row is the indices of lines skew to
        line i
    - ThreeLineStab is the ordered stabilizer of the list of lines with indices
        0,q+2,2q+4
    - LineAuto is the matrix whose (i,j)-th entry is the index of the line
        which the i-th element of ThreeLineStab sends the line with index j to.
    
    Outputs:
    Print number of lines on Fermat surface
    """
    global p, e, q, F, y, x, Lines, graph, Neighbours, ThreeLineStab, LineAuto
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
    power = e*2;
    F = galois.GF(p**power);
    y = F.primitive_element;
    if not(p==2):
        x = y**int((q-1)/2);
    else:
        x = y**(q-1);
    Lines = LineSet(q,x,y);
    graph = SkewnessMatrix(q);
    Neighbours = list(BuildNeighbours(graph));
    ThreeLineStab = list(BuildStab())
    LineAuto = list(BuildAuto())  

def main():
    global graph, ThreeLineStab, LineAuto, Neighbours
    setup()
    if q==2 or q==3 or q==4:
        R = [0,q+2,2*q+4]
        P = CommonNeighbours(R,graph)
        StabR = range(len(ThreeLineStab))
        isTrivial = False
    if q==2 or q==3:
        CliqueSizeList = ListSizeCounter(list(BKWithPivotingOrbit(R,P,[],graph,StabR,isTrivial)))
        for i in range(len(CliqueSizeList)):
            print("There are at most",CliqueSizeList[i][1],"orbits of maximal skew sets of size",CliqueSizeList[i][0],"containing the set",R,"up to the action of the stablizer of that skew set.")
        R = [0,q+2,2*q+4]
        P = CommonNeighbours(R,graph)
        CliqueSizeList =ListSizeCounter(list(BKWithPivoting(R,P,[],graph)))
        for i in range(len(CliqueSizeList)):
            print("There are",CliqueSizeList[i][1],"maximal skew sets of size",CliqueSizeList[i][0],"containing the set",R)
    elif q==4:
        Size=8;
        OrbitList = list(BKCliquesFixedSizeOrbit(Size,R,P,[],graph,range(len(ThreeLineStab)),False))
        for index in range(16):
            SkewSetList = [];
            for j in range(5000*index,min(5000*(index+1),len(OrbitList))):
                TempOrbits = list(BKWithPivotingOrbit(OrbitList[j][0],OrbitList[j][1],OrbitList[j][2],graph,OrbitList[j][3],OrbitList[j][4]))
                for SkewSet in TempOrbits:
                    SkewSetList.append(SkewSet);
            print("The count for this batch is as follows.")
            NewCliquesSizes = ListSizeCounter(SkewSetList)
            for i in range(len(NewCliquesSizes)):
                print("There are",NewCliquesSizes[i][1],"maximal skew sets of size",NewCliquesSizes[i][0],"containing the set",R)
    else:
        print("This program currently can only run without modification for q=2,3,4. Please input one of those integers or modify the code.")
        main()
    
    
if __name__ == "__main__":
    main()

