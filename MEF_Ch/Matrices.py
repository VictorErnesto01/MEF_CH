import numpy as np
from shapely import Polygon
from hat_utils import GradienteGorro
def MatrizMasa2D(P:np.array,T:np.array):
    """
    Construye la matriz de masa M de una malla 2D de elementos finitos lineales
    Lagrange P1

    Parameters
    ----------
    P : np.array
        Matriz de Nodos.
    T : np.array
        Matriz de conectividad.

    Returns
    -------
    M : np.array
        Matriz de Masa.

    """
    
    n = np.shape(P)[1]
    M = np.zeros((n,n))
    
    for i in range(np.shape(T)[1]):
        Tri = T[:,i]
        Nodos = []
        for e in range(len(Tri)):
            Nodos.append(P[:,Tri[e]])
        Nodos = np.array(Nodos)
        K = Polygon(Nodos).area
        Mk = K/12 * np.array([[2,1,1],[1,2,1],[1,1,2]])
        form = np.shape(Mk)
        for m in range(form[0]):
            for j in range(form[1]):
                M[Tri[m],Tri[j]] += Mk[m,j]
    return M
def MatrizMasa1D(x:np.array):
    """
    Construye la matriz de Masa para una malla de elemento finitos 1D.

    Parameters
    ----------
    x : np.array
        Malla de elementos finitos 1D.

    Returns
    -------
    M : np.array
        Matriz de Masa.

    """
    
    n = len(x)-1
    M = np.zeros((n+1,n+1))
    for i in range(0, n):
        h = x[i + 1] - x[i]
        M[i, i] += h / 3
        M[i, i + 1] += h / 6
        M[i + 1, i] += h / 6
        M[i + 1, i + 1] += h / 3
    return M
def MatrizRigidez1D(x:np.array,f=lambda x:1,kappa=[0,0]):
    """
    Calcula la matriz de Rigidez para un esquema de elementos finitos 1D

    Parameters
    ----------
    x : np.array
        Arreglo de elementos finitos.
    f : lamda, optional
        Función en caso de requerirse.
    kappa : lista de 2 elementos, optional
        Parametros por condiciones de Robin.

    Returns
    -------
    M : np.array
        Matriz de Rigidez.

    """
    
    n = len(x)-1
    M = np.zeros((n+1,n+1))
    for i in range(0, n):
        fa = f((x[i + 1] + x[i]) / 2)
        h = x[i + 1] - x[i]
        M[i, i] += fa / h
        M[i, i + 1] += -fa / h
        M[i + 1, i] += -fa/h
        M[i + 1, i + 1] += fa/h
    #Aplicar Condiciones de Frontera
    M[0,0] = kappa[0]
    M[-1,-1] = kappa[1]
    return M
def MatrizRigidez2D(P,T,f = lambda x,y:1):
    n = np.shape(P)[1]
    M = np.zeros((n,n))
    for i in range(np.shape(T)[1]):
        Tri = T[:,i]
        Nodos = []
        for e in range(len(Tri)):
            Nodos.append(P[:,Tri[e]])
        Nodos = np.array(Nodos)
        K = Polygon(Nodos).area
        b = []
        x = np.sum(Nodos[:,0])/3
        y = np.sum(Nodos[:,1])/3
        a= f(x,y)
        for i in range(len(Nodos)):
            b.append(GradienteGorro(Nodos, i))
        Mk = np.array([[b[0][0]**2+b[0][1]**2,b[0][0]*b[1][0]+b[0][1]*b[1][1],
                        b[0][0]*b[2][0]+b[0][1]*b[2][1]],
                       [b[1][0]*b[0][0]+b[0][1]*b[1][1]
                                            ,b[1][0]**2+b[1][1]**2,b[1][0]*b[2]
                                            [0]+b[1][1]*b[2][1]],
                       [b[0][0]*b[2][0]+b[0][1]*b[2][1],b[2][0]*b[1][0]
                        +b[2][1]*b[1][1],b[2][0]**2+b[2][1]**2,]])
        Mk = a * K * Mk
        form = np.shape(Mk)
        for m in range(form[0]):
            for j in range(form[1]):
                M[Tri[m],Tri[j]] += Mk[m,j]
    return M
    