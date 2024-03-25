import numpy as np
import scipy as scipy
import math as math
def read_matrices(filename: str,delimiter: str = '\t'):
    matrices = []
    current_matrix = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()  
            if not line or line.startswith("#"):
                if current_matrix: 
                    matrices.append(current_matrix)
                    current_matrix = []  
                continue
            
            try:
                row = [float(num) for num in line.split(delimiter)]
                current_matrix.append(row)
            except ValueError:
                # print("Skipping non-numeric line:", line)
                pass
        if current_matrix:
            matrices.append(current_matrix)
    return matrices
def print_matrix(matrix):

    for row in matrix:
        for element in row:
            print("\t",element, end="\t")
        print()

def check_symmetry(A):

    for i in range(len(A)):
        for j in range(len(A[i])):
            if i != j:
                if A[i][j] != A[j][i]:
                    return False
    return True
def get_transpose(A):

    for i in range(len(A)):
        for j in range(i+1,len(A[i])):
            A[i][j],A[j][i]=A[j][i],A[i][j]
    return A
def Cholesky_decompose(A):

    if check_symmetry(A) == False:
        raise ValueError('The matrix is not symmetric')
      
    from math import sqrt
    if len(A) != len(A[0]):
        return None
    n = len(A)

    for row in range(n):
        for col in range(row,n): 
            if row == col:
                sum = 0
                for i in range(row):
                    sum += (A[row][i] ** 2)
                A[row][row] = sqrt(A[row][row] - sum)
            else:
                sum = 0
                for i in range(row):
                    sum += (A[row][i] * A[i][col])
                A[row][col] = (A[row][col] - sum) / A[row][row]
                A[col][row] = A[row][col]   
    for row in range(n):
        for col in range(row + 1,n):
            A[row][col] = 0   
    return A

def sup_chol_dec(A):
    L=Cholesky_decompose(A)
    LT=get_transpose(np.copy(L).tolist())
    LP=[[0 for i in range(len(L))] for i in range(len(L))]
    for i in range(len(L)):
        for j in range(len(L)):
            LP[i][j]=L[i][j]+LT[i][j]
    for i in range(len(L)):
        LP[i][i]=LP[i][i]/2
    return LP   
def Cholesky_solve(A,B):
    LP=sup_chol_dec(np.copy(A))
    del A
    A=LP
    Y = []
    n = len(B)
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += (A[i][j] * Y[j])
        Y.append((B[i][0]-sum)/A[i][i])

    X = Y
    for i in range(n-1,-1,-1):
        sum = 0
        for j in range(i + 1, n):
            sum += (A[i][j] * X[j])
        X[i] = (Y[i] - sum) / A[i][i]

    for i in range(n):
        X[i] = [X[i]]

    return Y

def Check_diagonal_dominance(m) :
    n=len(m)
    for i in range(0, n) :        
        sum = 0
        for j in range(0, n) :
            sum = sum + abs(m[i][j])    
        sum = sum - abs(m[i][i])
        if (abs(m[i][i]) < sum) :
            return False 
    return True
def gauss_seidel(A,B,tolerance):

    n = len(A)
    X = list([0] for i in range(n))
    for steps in range(200):
        flag = 1
        for i in range(n):
            sum = 0
            for j in range(i):
               sum += (A[i][j] * X[j][0])
            for j in range(i+1,n):
                sum += (A[i][j] * X[j][0])
            temp = (B[i][0] - sum) / (A[i][i])
            if (abs(temp) - abs(X[i][0])) > tolerance:
               flag = 0
            X[i][0] = temp
        if flag == 1:
            return X,steps + 1
    print('Eqn not solved after 200 steps')
    return None,200

def Gauss_seidel_solve(A,B,T):

    if Check_diagonal_dominance(A)==True:
        return gauss_seidel(A,B,T)
    else:
        raise ValueError('The matrix is not diagonally dominant')  
    
def gauss_jordan_solve(A,B):
    augmat = A #constructing augmented matrix
    for row in range(len(augmat)):
        augmat[row].append(B[row][0])
    for curr in range(len(augmat)): #curr takes the index of each column we are processing
        if augmat[curr][curr] == 0: #row swap if zero
            max_row = curr
            for row in range(curr + 1,len(augmat)):
                if abs(augmat[row][curr]) > abs(augmat[max_row][curr]):
                    max_row = row
            if max_row == curr: #if max elemnt is still zero, max_row is not changed; no unique solution
                return None
            augmat[curr],augmat[max_row] = augmat[max_row], augmat[curr]
        #making the pivot element 1
        if augmat[curr][curr] != 1:
            pivot_term = augmat[curr][curr]
            for i in range(len(augmat[curr])):
                augmat[curr][i] = augmat[curr][i] / pivot_term
        #making others zero
        for i in range(0,len(augmat)):
            if i == curr: #skipping the pivot point
                continue
            if augmat[i][curr] != 0:
                lead_term = augmat[i][curr]
                for j in range(curr,len(augmat[i])): #elements before the curr column are zero in curr row, so no need to calculate
                    augmat[i][j] = augmat[i][j] - (augmat[curr][j] * lead_term)
    solution = []
    for i in range(len(augmat)):
        solution.append([augmat[i][-1]]) #Taking last elements into a list to form column matrix
    return solution 

def LU_decompose(A):
    '''
    # LU Decomposition of a Matrix
    '''
    n = len(A) 
    if n != len(A[0]): 
        print('Not square!')
        return None

    for j in range(n):
        for i in range(1,n):
            if i <= j:
                sum = 0
                for k in range(0,i):
                    sum += (A[i][k] * A[k][j])
                A[i][j] = A[i][j] - sum
            else:
                sum = 0
                for k in range(0,j):
                    sum += (A[i][k] * A[k][j])
                A[i][j] = (A[i][j] - sum) / A[j][j]
    return A

def for_back_LU(A,B):
    '''
    # Forward and Backward Substitution for LU Decomposition
    '''
    Y = []
    n = len(B)
    for i in range(n):
        s = 0
        for j in range(i):
            s += (A[i][j] * Y[j])
        Y.append(B[i][0]-s)
    X = Y[:]
    for i in range(n-1,-1,-1):
        s = 0
        for j in range(i + 1, n):
            s+= (A[i][j] * X[j])
        X[i] = (Y[i] - s) / A[i][i]

    for i in range(n):
        X[i] = [X[i]]

    return X
def Make_diag_dominant(A,B):
    '''
    # Making the matrix A Diagonally Dominant
    Returns the Diaonally Dominant matrix A and the corresponding B
    '''
    n = len(A)
    for i in range(n):
        sum=0
        for j in range(n):
            if j != i:
                sum += abs(A[i][j])
        if abs(A[i][i])>sum:
            continue
        else:
            c = i + 1
            flag = 0
            while c<n:
                sum = 0
                for j in range(n):
                    if j != i:
                        sum += abs(A[c][j])
                if abs(A[c][i])>sum:
                    flag = 1
                    break
                else:
                    c+=1
            if flag==0:
                return None,None
            else:
                A[i],A[c]=A[c],A[i]
                B[i],B[c]=B[c],B[i]
    return A,B


def LU_Solve_eqn(A,B):
    '''
    # LU Decomposition Method for solving the equation A.X = B
    '''
    # A,B = make_diag_non_0(A,B)
    A,B = Make_diag_dominant(A,B)
    L=LU_decompose(A)
    L1=for_back_LU(L,B)
    return L1

def conjugate_gradient_solve(A: list,B: list,guess: list,T: float):
    x0=guess
    r0 = np.add(B, -1 * np.matmul(A, x0))
    d0 = np.copy(r0)
    i=1
    while True:
        alpha1 = np.matmul(np.transpose(r0), r0) / np.matmul(np.transpose(d0), np.matmul(A, d0))
        x1 = np.add(x0, alpha1[0][0]*d0)
        r1 = np.add(r0, -1 * alpha1[0][0] * np.matmul(A, d0))
        if np.linalg.norm(r1) < T and i<=len(A):
            return x1.tolist(),i
        
        elif i>len(A):
            print("Maybe the matrix A dosent satisfy the conditions for the Conjugate Gradient Method")
            return None
        else: 
            beta1 = np.matmul(np.transpose(r1), r1) / np.matmul(np.transpose(r0), r0)
            d1 = np.add(r1, beta1[0][0] * d0)
            x0 = x1
            del x1
            r0 = r1
            del r1
            d0 = d1
            del d1
            i+=1

