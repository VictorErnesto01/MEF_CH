from setuptools import setup, find_packages

setup(
    name='MEF_Ch',
    version='0.1',
    packages=find_packages(),
    license='MIT',
    description='Varias utilidades para resolver Ecuaciones Diferenciales con Elemento Finito',
    install_requires=["numpy",
                      "matplotlib.pyplopt",
                      "gmsh",
                      "shapely",
                      "meshio"],
    tests_require=['pytest'],
    setup_requires=['pytest-runner'],
    test_suite='tests',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
