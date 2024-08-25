import numpy as np

def brent(f, a, b, tol=1.0e-9):
    x1 = a
    x2 = b
    f1 = f(x1)
    
    if f1 == 0.0: 
        return x1
    
    f2 = f(x2)
    if f2 == 0.0: 
        return x2
    
    if f1 * f2 > 0.0: 
        raise ValueError('Root is not bracketed')
    
    x3 = 0.5 * (a + b)
    
    for i in range(30):
        f3 = f(x3)
        
        if abs(f3) < tol: 
            return x3
        
        # Tighten the brackets on the root
        if f1 * f3 < 0.0: 
            b = x3
        else: 
            a = x3
        
        if (b - a) < tol * max(abs(b), 1.0): 
            return 0.5 * (a + b)
        
        # Try quadratic interpolation
        denom = (f2 - f1) * (f3 - f1) * (f2 - f3)
        
        # Handle potential division by zero
        if denom != 0:
            numer = x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 * (f2 - f3) + f1 * x2 * (f3 - f1)
            dx = f3 * numer / denom
        else:
            dx = b - a
        
        x = x3 + dx
        
        # If interpolation goes out of bounds, use bisection
        if (b - x) * (x - a) < 0.0:
            dx = 0.5 * (b - a)
            x = a + dx
        
        # Let x3 <-- x & choose new x1 and x2 so that x1 < x3 < x2
        if x < x3:
            x2 = x3
            f2 = f3
        else:
            x1 = x3
            f1 = f3
        
        x3 = x
    
    print('Too many iterations in brent')
    return None

def rootsearch(f, a, b, dx):
    x1 = a
    f1 = f(x1)
    
    x2 = a + dx
    f2 = f(x2)
    
    while x2 <= b:
        if f1 * f2 < 0.0:
            return x1, x2
        x1 = x2
        f1 = f2
        x2 = x1 + dx
        f2 = f(x2)
    
    return None, None

# Define the function
def f(x):
    return x * np.sin(x) + 3 * np.cos(x) - x

# Interval and step size
a, b = -6, 6
dx = 0.1

# Finding all roots in the interval
roots = []
while True:
    x1, x2 = rootsearch(f, a, b, dx)
    if x1 is not None and x2 is not None:
        root = brent(f, x1, x2)
        if root is not None:
            roots.append(root)
        a = x2  # Move to the next subinterval
    else:
        break

# Display the roots
print("Roots found:")
for root in roots:
    print(root)
