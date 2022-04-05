 # Modern Computing Technologies, WUT, 2022
 # 
 # Author: Jakub Zieliński
 # Date: 30.03.2022
 #
 # Calculate covariance using precompiled
 # library in c
 #
from numpy.ctypeslib import ndpointer
from ctypes import *
import numpy as np
import time

NELEMENTS = 134217728
libcov = CDLL("./libcov.so")

# check connection
libcov.connect()

# Describe function
calculate_covariance_1Darray = libcov.calculate_covariance_1Darray
calculate_covariance_1Darray.restype = None
calculate_covariance_1Darray.argtypes = [c_long, ndpointer(dtype=np.double, ndim=1, shape=(10*NELEMENTS), flags='C_CONTIGUOUS'), ndpointer(dtype=c_double, shape=(10*10), flags='C_CONTIGUOUS')]
                                      
# read_data_test = libcov.read_data_test
# read_data_test.restype = None
# read_data_test.argtypes = [c_long,c_int, ndpointer(dtype=np.double, ndim=1, shape=(10*NELEMENTS))]

# create covariance array as 1D array
cov = np.zeros(shape=(10*10),dtype=np.double)
print("cov shape= {0}, is contiguous = {1}".format(cov.shape, cov.flags['C_CONTIGUOUS']))
#avg = np.zeros(10, dtype = np.double)

# create variables array as 1D array
var = np.empty(shape=(10*NELEMENTS), dtype=np.double) 
fileName = "/home2/archive/MCT-2022/lab1/var{0}.dat"
# read data
for i in range(5):
    var[i*NELEMENTS:(i+1)*NELEMENTS] = np.fromfile(fileName.format(i+1), dtype=np.double)
    print("File {0} read: [10] = {1}".format(i, var[i*NELEMENTS+10]))

print("var shape= {0}, is contiguous = {1}".format(var.shape, var.flags['C_CONTIGUOUS']))
# print("Inside c function 10th elements are:\n")
# read_data_test(NELEMENTS, 0, var);
# read_data_test(NELEMENTS, 1, var);
# read_data_test(NELEMENTS, 2, var);
# read_data_test(NELEMENTS, 3, var);
# read_data_test(NELEMENTS, 4, var);
# read_data_test(NELEMENTS, 5, var);
# read_data_test(NELEMENTS, 6, var);
# read_data_test(NELEMENTS, 7, var);
# read_data_test(NELEMENTS, 8, var);
# read_data_test(NELEMENTS, 9, var);
# calculate covariance
bt = time.time()
calculate_covariance_1Darray(NELEMENTS,var,cov)
et = time.time()
print("# COMPUTATION TIME: {0} sec".format(et-bt))

# changing data to 2D array in the end
var = var.reshape(10,NELEMENTS)
print("Check 10th elements ", var.shape)
for row in var:
    print(row[10], end=",")
print()

cov = cov.reshape((10,10))
print(cov.shape)
# print results
for i in range(10):
    for j in range(i+1):
        print("cov[{0:2d}][{1:2d}]={2:16.8f}".format(i+1, j+1, cov[i][j]))
