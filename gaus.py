import numpy as np

def gaussElimin(a, b):
    n = len(b)
    
    # Elimination phase
    for k in range(0, n-1):
        for i in range(k+1, n):
            if a[i, k] != 0.0:
                lam = a[i, k] / a[k, k]
                a[i, k+1:n] = a[i, k+1:n] - lam * a[k, k+1:n]
                b[i] = b[i] - lam * b[k]
    
    # Back substitution
    for k in range(n-1, -1, -1):
        b[k] = (b[k] - np.dot(a[k, k+1:n], b[k+1:n])) / a[k, k]
    
    return b

# Example usage:
np.set_printoptions(suppress=True) #to convert very small values like e-16 to 0
a = np.array([[2.0, 1.0, -4.0],
              [1.0, -1.0, 1.0],
              [-1.0, 3.0, -2.0]])

b = np.array([[1.0, 0, 0],
              [0, 1.0,0],
              [0, 0, 1.0]])

# X can be found through X = B * A^(-1) 
ainv = np.linalg.inv(a) # = A^(-1)
x = np.matmul(b, ainv) #matrix multiplication
print("Matrix X is\n", x)

# check by substituting
xcheck = np.matmul(a, x)
print("We check if matrix X gives us matrix B result\n", xcheck)