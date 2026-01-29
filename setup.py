from setuptools import setup, Extension
import pybind11
import sys
import os

# Detecção de Sistema Operacional
if sys.platform == 'win32':
    # Configurações para Visual Studio (Windows)
    cpp_args = ['/O2', '/openmp', '/EHsc', '/std:c++14']
    link_args = []
else:
    # Configurações para GCC (Linux/AWS)
    cpp_args = ['-O3', '-fopenmp', '-std=c++11']
    link_args = ['-fopenmp']

ext_modules = [
    Extension(
        'imm_module',
        # Liste aqui todos os seus .cpp, EXCETO o imm.cpp (pois ele já foi incluído no bindings.cpp)
        # Se models.cpp e tools.cpp já estão inclusos no imm.cpp via #include, não coloque aqui.
        # Se eles são compilados separadamente, adicione-os na lista abaixo.
        # Assumindo que você usa #include "models.cpp" dentro do imm.cpp (comum em scripts unicos):
        ['bindings.cpp', 'models.cpp', 'tools.cpp', 'mt19937ar.cpp'],
        
        include_dirs=[pybind11.get_include()],
        language='c++',
        extra_compile_args=cpp_args,
        extra_link_args=link_args,
    ),
]

setup(
    name='imm_module',
    version='1.0',
    ext_modules=ext_modules,
)