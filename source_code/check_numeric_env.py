import sys
import platform
import numpy as np
import scipy
import pandas as pd

print("Python executable:", sys.executable)
print("Python version:", sys.version)
print("Platform:", platform.platform())

print("\nPackage versions")
print("NumPy:", np.__version__)
print("SciPy:", scipy.__version__)
print("pandas:", pd.__version__)

print("\nDefault NumPy dtypes")
print("np.array([1.0]).dtype:", np.array([1.0]).dtype)
print("np.zeros(1).dtype:", np.zeros(1).dtype)
print("np.ones(1).dtype:", np.ones(1).dtype)
print("np.empty(1).dtype:", np.empty(1).dtype)
print("np.linspace(0.0, 1.0, 3).dtype:", np.linspace(0.0, 1.0, 3).dtype)

print("\nFloating point info")
print("np.finfo(float).bits:", np.finfo(float).bits)
print("np.finfo(np.float64).bits:", np.finfo(np.float64).bits)
print("np.finfo(np.float32).bits:", np.finfo(np.float32).bits)

print("\nBLAS/LAPACK configuration")
np.show_config()