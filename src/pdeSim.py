import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix

class PoissonFlow2D:
    """
    Contains a sparse solve method for the steady state of the advection-diffusion equation.

    """	

    def __init__(self, Q, kappa, v, dx,dy,c,rho,freq):        
        self.Q=Q
        self.kappa=kappa
        self.v=v
        self.dx=dx
        self.dy=dy
        self.c=c
        self.rho=rho        
        self.freq=freq
        
    def sparseSolve(self):
            row = []
            col = []
            data = []
            nx = self.kappa.shape[0]
            ny = self.kappa.shape[1]
          

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
                          
            Q_in_flat = np.zeros(nx*ny)
            for j in range(ny):
                for i in range(nx):
                    Q_in_flat[nx*j+i] = self.Q[i,j]
            

            n = np.zeros(nx*ny,dtype=complex)
            n = spsolve(As,-Q_in_flat)
            
            unwrapped = np.zeros((nx,ny),dtype=complex)
            for j in range(ny):
                for i in range(nx):
                    unwrapped[i,j]= n[nx*j+i]
        
            return unwrapped
            
