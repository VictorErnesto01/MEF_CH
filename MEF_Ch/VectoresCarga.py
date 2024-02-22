import numpy as np
import traceback
from hat_utils import PuntosMedios,FuncionGorro
from shapely import Polygon
def VectorCarga1D(x:np.array,f,Integrador="Trapecio"):
    """
    Calcula el Vector de Carga de unm sistema de Elemento Finito unidimensional
    de acuerdo a la Proyección L2 en una base nodal.

    Parameters
    ----------
    x : Numpy Array.
        Arreglo de Elementos finitos.
    f : Lambda.
        Función.
    Integrador : String, opcional.
        El integrador a utilizar, por defecto se utiliza la regla del trapecio.
        "Trapecio" para regla del trapecio
        "Medio" para Punto Medio
        "Simpson" para Simpson 1/3

    Returns
    -------
    Vector de Carga correspondiente como un numpy array.

    """
    def Trapecio(x,f):
        n = len(x)-1
        b=np.zeros((n+1))
        for i in range(0, n):
            h = x[i + 1] - x[i]
            b[i] += f(x[i]) * h / 2
            b[i + 1] += f(x[i + 1]) * h / 2
        return b
    def Medio(x,f):
        n = len(x)-1
        b=np.zeros((n+1))
        for i in range(0, n):
            h = (x[i + 1] - x[i])
            m = (x[i + 1] + x[i])/2
            b[i] += f(m) * h/2
            b[i+1] += f(m) * h/2
        return b
    def Simpson(x,f):
        n = len(x)-1
        b=np.zeros((n+1))
        for i in range(0,n):
            h = x[i + 1] - x[i]
            m = (x[i + 1] + x[i])/2
            b[i] += h/6 * (f(x[i])+2*f(m))
            b[i+1] += h/6 * (f(x[i])+2*f(m))
        return b
    try:
        Cuadratura = {
            "Trapecio": Trapecio,
            "Medio": Medio,
            "Simpson": Simpson
            }[Integrador](x,f)
        return Cuadratura
    except KeyError:
        traceback.print_exc()
    except Exception:
        traceback.print_exc()
def VectorCarga2D(P:np.array,T:np.array,f,Integrador = "Simpson"):
    """
    Calcula el vector de carga para una malla lineal de Polinomios de Lagrange
    P1.

    Parameters
    ----------
    P : np.array
        Malla de Nodos.
    T : np.array
        Malla de Conectividad.
    f : Lambda
        Función de x,y.
    Integrador : String, opcional.
        El integrador a utilizar, por defecto se utiliza la regla del trapecio.
        "Trapecio" para regla del trapecio
        "Medio" para Punto Medio
        "Centro" para Centro de Masa

    Returns
    -------
    Vector de Carga correspondiente como un numpy array.

    """
    def Trapecio(P,T,f = lambda x,y:0):
        n = np.shape(P)[1]
        b = np.zeros(n)
        bk = 1/3 * np.array([1,1,1])
        for i in range(np.shape(T)[1]):
            Tri = T[:,i]
            Nodos = []
            for j in range(3):
                Nodos.append(P[:,Tri[j]])
            Nodos = np.array(Nodos)
            K = Polygon(Nodos).area
            bk2 = bk * K
            for m in range(len(bk2)):
                b[Tri[m]] += f(Nodos[m,0],Nodos[m,1])*bk2[m]
        return b
    def Centro(P,T,f = lambda x,y:0):
        n = np.shape(P)[1]
        b = np.zeros(n)
        bk = np.array([1,1,1])
        for i in range(np.shape(T)[1]):
            Tri = T[:,i]
            Nodos = []
            for j in range(3):
                Nodos.append(P[:,Tri[j]])
            Nodos = np.array(Nodos)
            K = Polygon(Nodos).area
            bk2 = bk * K
            for m in range(len(bk2)):
                hat = FuncionGorro(Nodos, m)
                x = np.sum(Nodos[:,0])/3
                y = np.sum(Nodos[:,1])/3
                b[Tri[m]] += f(x,y)*hat(x,y)*bk2[m]
        return b
    def Medio(P,T,f = lambda x,y:0):
        n = np.shape(P)[1]
        b = np.zeros(n)
        bk = 1/3*np.array([1,1,1])
        for i in range(np.shape(T)[1]):
            Tri = T[:,i]
            Nodos = []
            for j in range(3):
                Nodos.append(P[:,Tri[j]])
            Nodos = np.array(Nodos)
            K = Polygon(Nodos).area
            bk2 = bk * K
            xm = PuntosMedios(Nodos[:,0])
            ym = PuntosMedios(Nodos[:,1])
            for m in range(len(bk2)):
                hat = FuncionGorro(Nodos, m)
                b[Tri[m]] += np.sum(f(xm,ym)*hat(xm,ym))*bk2[m]
        return b
    try:
        Cuadratura = {
            "Trapecio": Trapecio,
            "Medio": Medio,
            "Centro": Centro
            }[Integrador](P,T,f)
        return Cuadratura
    except KeyError:
        traceback.print_exc()
    except Exception:
        traceback.print_exc()
    
    