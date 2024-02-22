import numpy as np
def PuntosMedios(x):
    """
    Calcula los puntos medios de cada arista de un triangulo dadas sus coordenadas x,y
    """
    xm0 = x[0] + x[1]
    xm1 = x[0] + x[2]
    xm2 = x[1] + x[2]
    xm = np.array([xm0, xm1, xm2]) / 2
    return xm
def FuncionGorro(Nodos,i):
    M = np.ones((3,3))
    M[:,1:] = Nodos
    b = np.zeros(3)
    b[i] = 1
    u = np.linalg.solve(M,b)
    return lambda x,y: u[0]+u[1]*x+u[2]*y
def GradienteGorro(Nodos,i):
    M = np.ones((3,3))
    M[:,1:] = Nodos
    b = np.zeros(3)
    b[i] = 1
    u = np.linalg.solve(M,b)
    return [u[1],u[2]]