import gmsh
import meshio
import os
import datetime
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
class Malla:
    def __init__(self,fecha,ruta,nombre,user):
        self.__fecha = fecha
        self.ruta = ruta
        self.__nombre = nombre
        self.__malla = meshio.read(ruta)
        self.P = self.__malla.points[:, :2].T
        self.T = self.__malla.cells_dict["triangle"].T
        self.__t = self.__malla.cells_dict["triangle"].T 
        self.__frontera = self.__malla.cells_dict["line"]
        self.__nodos = np.shape(self.P)[1]
        self.__elementos = np.shape(self.T)[1]
        self.info = {
            "Creador" : "Esta malla fue diseñada con el código de prueba",
            "Nombre de la malla" : nombre,
            "Fecha de Creación" : fecha,
            "Ruta de acceso": ruta,
            "Elementos de la malla": self.__elementos,
            "Nodos de la malla": self.__nodos,
            "Creador de la malla: ": user
            }
        
    def DibujarMalla(self, label_triangles=False):
        # Load the mesh
        mesh = meshio.read(self.ruta)

        # Extract points and triangle connectivity
        points = mesh.points
        triangles = mesh.cells_dict["triangle"]

        # Create a matplotlib Triangulation object
        triangulation = tri.Triangulation(points[:, 0], points[:, 1], triangles)

        # Plot the triangulation
        plt.figure()
        plt.triplot(triangulation, "ko-")

        # Optionally label triangles
        if label_triangles:
            for i, triangle in enumerate(triangles):
                x = points[triangle, 0].mean()
                y = points[triangle, 1].mean()
                plt.text(x, y, str(i), ha="center", va="center", color="red")

        plt.xlabel("X")
        plt.ylabel("Y")
        plt.axis("equal")
        plt.title(f"Mesh Plot: {self.__nombre}")
        plt.show()
    def boundary_nodes(self):
        nodes = np.unique(self.__frontera.flatten())
        return nodes
def MeshCreator2D(x,y, lcs,archivo="malla.msh"):
    if os.path.exists(archivo):
        entr = input("El archivo ya existe, desea reemplazarlo: [S/N]: ")
        if entr == "S" or entr == "s":
            archivo = archivo
        else:
            archivo = str(archivo.split(".")[0] + "_nuevo.msh")
            print("Se guardo como un nuevo archivo, para evitar eliminarlo")
    #Inicializa el mallador
    gmsh.initialize()
    #Carga un modelo vacío
    # Ensure the input lists are of the same length
    assert len(x) == len(y) == len(lcs), "x, y, and lcs must be the same length."

    # Create points in the model
    for xi, yi, lci in zip(x, y, lcs):
        gmsh.model.geo.addPoint(xi, yi, 0, lci)

    # Define the lines connecting the points
    num_points = len(x)
    for i in range(num_points):
        gmsh.model.geo.addLine(i + 1, (i + 1) % num_points + 1)

    # Create a curve loop and a plane surface
    curve_loop = gmsh.model.geo.addCurveLoop(list(range(1, num_points + 1)))
    gmsh.model.geo.addPlaneSurface([curve_loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Generate the mesh
    gmsh.model.mesh.generate(2)

    # Save the mesh to a file
    #Guardar la malla a un archivo
    gmsh.write(archivo)
    
    #Extraer Información
    gmsh.finalize()
    username = os.getlogin()
    fecha = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    ruta = os.path.abspath(archivo)
    Mesh = Malla(fecha, ruta, archivo, username)
    
    return Mesh
def MeshLoader2D(archivo:str):
    """
    Carga un archivo .msh a Python

    Parameters
    ----------
    archivo : str
        Ruta de acceso

    Returns
    -------
    Mesh.

    """
    return Malla("Proximamente", os.path.abspath(archivo), archivo, "Proximamente")
#%%
# x = [0,1,1,0]
# y = [0,0,1,1]
# lc = (0.02,0.2,0.2,0.2)
# Mesh = MeshCreator2D(x,y,lc,"prueba.msh")
# #%%
# P = Mesh.P
# T = Mesh.T
# Mesh.DibujarMalla()
# #%%
# Mesh = MeshLoader2D("prueba.msh")
# P = Mesh.P
# Frontera = Mesh.boundary_nodes()

    
