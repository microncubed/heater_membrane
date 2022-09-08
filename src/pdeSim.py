import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix

class PoissonFlow2D:
    """
    Contains a sparse solve method for the steady state of the advection-diffusion equation.
    Dirichlet (constant temperature) boundary conditions are used all around.
    """	

    def __init__(self, Q, kappa, v, dx,dy,c,rho,freq):
        '''
        Parameters
        ----------
        Q np.ndarray: heat source
        kappa np.ndarray: thermal conductivity (W/m/K)
        v np.ndarray: velocity field (m/s)
        dx float: step-size in x (m)
        dy float: step-size in y (m)
        c np.ndarray: heat capacity (J/K/m3)
        rho np.ndarray: density (kg/m3)
        freq float: non-zero in case a frequency domain simulation is performed (Hz)
        '''        
        self.Q=Q
        self.kappa=kappa
        self.v=v
        self.dx=dx
        self.dy=dy
        self.c=c
        self.rho=rho        
        self.freq=freq
        
    def sparseSolve(self):
        '''
        Solves the Poisson equation with the given parameters
        
        Returns
        ----------
        unwrapped np.ndarray(dtype = complex): the result of the simulation 
        '''
        row = []
        col = []
        data = []
        nx = self.kappa.shape[0]
        ny = self.kappa.shape[1]
          
        # build up the COO sparse matrix
        for j in range(ny):
            for i in range(nx):
                Adiag = j*nx+i
                if j==0:
                    self.Q[i,j] = 0
                    row.append(Adiag)
                    col.append(Adiag)
                    data.append(1)
                elif j ==(ny-1):
                    self.Q[i,j] = 0
                    row.append(Adiag)
                    col.append(Adiag)
                    data.append(1)
                elif i ==0:
                    self.Q[i,j] = 0
                    row.append(Adiag)
                    col.append(Adiag)
                    data.append(1)
                elif i==(nx-1):
                    self.Q[i,j] = 0
                    row.append(Adiag)
                    col.append(Adiag)
                    data.append(1)
                else:
                    data.append(-1/2*(self.kappa[i-1,j]+self.kappa[i-1,j-1])*self.dy/self.dx-\
                    1/2*self.c[i,j]*self.rho[i,j]*self.v[i,j]*self.dy-\
                    1/2*(self.kappa[i,j]+self.kappa[i,j-1])*self.dy/self.dx+\
                    1/2*self.c[i,j]*self.rho[i,j]*self.v[i,j]*self.dy-\
                    1/2*(self.kappa[i-1,j]+self.kappa[i,j])*self.dx/self.dy-\
                    1/2*(self.kappa[i,j-1]+self.kappa[i-1,j-1])*self.dx/self.dy-\
                    self.c[i,j]*self.rho[i,j]*self.freq*np.pi*2*self.dx*self.dy*1j)
                        
                    row.append(Adiag)
                    col.append(Adiag)
                        
                    data.append(1/2*(self.kappa[i-1,j]+self.kappa[i-1,j-1])*self.dy/self.dx-\
                    1/2*self.c[i,j]*self.rho[i,j]*self.v[i,j]*self.dy)
                    row.append(Adiag)
                    col.append(Adiag-1)
                        
                    data.append(1/2*(self.kappa[i-1,j]+self.kappa[i,j])*self.dx/self.dy)
                    row.append(Adiag)
                    col.append(Adiag+nx)
                        
                    data.append(1/2*(self.kappa[i,j]+self.kappa[i,j-1])*self.dy/self.dx+\
                    1/2*self.c[i,j]*self.rho[i,j]*self.v[i,j]*self.dy)
                    row.append(Adiag)
                    col.append(Adiag+1)
                        
                    data.append(1/2*(self.kappa[i,j-1]+self.kappa[i-1,j-1])*self.dx/self.dy)
                    row.append(Adiag)
                    col.append(Adiag-nx)
                                          
        As=coo_matrix((data, (row, col)), shape=(nx*ny, nx*ny))
        As = As.tocsr()
                          
        # flatten the input heat source array    
        Q_in_flat = np.zeros(nx*ny)
        for j in range(ny):
            for i in range(nx):
                Q_in_flat[nx*j+i] = self.Q[i,j]
            

        #solve and wrap up in a 2-d array
        n = spsolve(As,-Q_in_flat)   
        wrapped = np.zeros((nx,ny),dtype=complex)
        for j in range(ny):
            for i in range(nx):
                wrapped[i,j]= n[nx*j+i]
        
        return wrapped
            
