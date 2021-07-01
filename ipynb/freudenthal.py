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