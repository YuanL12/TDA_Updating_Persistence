import bats

def freudenthal_grid(m, n):
    """
    Freudenthal triangulation of a m x n grid
    """
    
    def _get_idx(i, j):
        """
        get index of grid in row-major order
        """
        return j + n * i;
    
    
    X = bats.SimplicialComplex()
    
    for i in range(m-1):
        for j in range(n-1):
            k1 = _get_idx(i,j)
            k2 = _get_idx(i+1,j)
            k3 = _get_idx(i,j+1)
            k4 = _get_idx(i+1,j+1)
            X.add_recursive([k1,k2,k4])
            X.add_recursive([k1,k3,k4])

    return X

def freudenthal_grid_light(m, n):
    """
    Freudenthal triangulation of a m x n grid
    """
    
    def _get_idx(i, j):
        """
        get index of grid in row-major order
        """
        return j + n * i;
    
    
    X = bats.LightSimplicialComplex(m*n, 2)
    
    for i in range(m-1):
        for j in range(n-1):
            k1 = _get_idx(i,j)
            k2 = _get_idx(i+1,j)
            k3 = _get_idx(i,j+1)
            k4 = _get_idx(i+1,j+1)
            X.add_recursive([k1,k2,k4])
            X.add_recursive([k1,k3,k4])

    return X

def levelset_cube(n1, n2, n3):
    """
    Freudenthal triangulation of a n1 x n2 x n3 grid
    """
    
    def _get_idx(i, j, k):
        """
        get index of grid in row-major order
        """
        return k + n2 * (j + n1 * i);
    
    
    X = bats.SimplicialComplex()
    
    for i in range(n1-1):
        for j in range(n2-1):
            for k in range(n3-1):
                s1 = _get_idx(i,j,k)
                s2 = _get_idx(i,j,k+1)
                s3 = _get_idx(i,j+1,k)
                s4 = _get_idx(i+1,j,k)
                s5 = _get_idx(i,j+1,k+1)
                s6 = _get_idx(i+1,j,k+1)
                s7 = _get_idx(i+1,j+1,k)
                s8 = _get_idx(i+1,j+1,k+1)
                X.add_recursive([s1,s2,s5,s8]) #k, j, i
                X.add_recursive([s1,s2,s6,s8]) #k, i, j
                X.add_recursive([s1,s4,s7,s8]) #i, j, k
                X.add_recursive([s1,s4,s6,s8]) #i, k, j
                X.add_recursive([s1,s3,s7,s8]) #j, i, k
                X.add_recursive([s1,s3,s5,s8]) #j, k, i

    return X

def levelset_cube_light(n1, n2, n3):
    """
    Freudenthal triangulation of a n1 x n2 x n3 grid
    """
    
    def _get_idx(i, j, k):
        """
        get index of grid in row-major order
        """
        return k + n2 * (j + n1 * i);
    
    
    X = bats.LightSimplicialComplex(n1*n2*n3, 3) # 3 is the max complex dimension
    
    for i in range(n1-1):
        for j in range(n2-1):
            for k in range(n3-1):
                s1 = _get_idx(i,j,k)
                s2 = _get_idx(i,j,k+1)
                s3 = _get_idx(i,j+1,k)
                s4 = _get_idx(i+1,j,k)
                s5 = _get_idx(i,j+1,k+1)
                s6 = _get_idx(i+1,j,k+1)
                s7 = _get_idx(i+1,j+1,k)
                s8 = _get_idx(i+1,j+1,k+1)
                X.add_recursive([s1,s2,s5,s8]) #k, j, i
                X.add_recursive([s1,s2,s6,s8]) #k, i, j
                X.add_recursive([s1,s4,s7,s8]) #i, j, k
                X.add_recursive([s1,s4,s6,s8]) #i, k, j
                X.add_recursive([s1,s3,s7,s8]) #j, i, k
                X.add_recursive([s1,s3,s5,s8]) #j, k, i

    return X


def cubical_grid(m,n):
    X = bats.CubicalComplex(2)
    for i in range(m-1):
        for j in range(n-1):
            X.add_recursive([i, i+1, j, j+1])
            
    return X


def cubical_cube(n1,n2,n3):
    X = bats.CubicalComplex(3)
    for i in range(n1-1):
        for j in range(n2-1):
            for k in range(n3-1):
                X.add_recursive([i, i+1, j, j+1, k, k+1])

    return X