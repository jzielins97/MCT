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
calculate_covariance_1Darray.restype = ndpointer(c_double)
calculate_covariance_1Darray.argtypes = [c_long, ndpointer(dtype=np.double, ndim=1, shape=(10*NELEMENTS), flags='C_CONTIGUOUS'), ndpointer(dtype=c_double, shape=(10*10), flags='C_CONTIGUOUS')]

calculate_other_variables = libcov.calculate_other_variables
calculate_other_variables.restype = ndpointer(c_double)
calculate_other_variables.argtypes = [c_long,
                                      ndpointer(c_double), #var[0]
                                      ndpointer(c_double), #var[1]
                                      ndpointer(c_double), #var[2]
                                      ndpointer(c_double), #var[3]
                                      ndpointer(c_double), #var[4]
                                      ndpointer(c_double), #var[5]
                                      ndpointer(c_double), #var[6]
                                      ndpointer(c_double), #var[7]
                                      ndpointer(c_double), #var[8]
                                      ndpointer(c_double)] #var[9]

calculate_covariance = libcov.calculate_covariance
calculate_covariance.restype = ndpointer(c_double)
calculate_covariance.argtypes = [c_long, ndpointer(c_double), ndpointer(c_double), ndpointer(c_double)]
                                      

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
    print("File {0} read {1} elements: [10] = {2}".format(i,var.shape, var[i*NELEMENTS+10]))

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
# calculate_covariance(NELEMENTS,var[0], var[1], var[2], var[3], var[4], var[5], var[6], var[7], var[8], var[9])
calculate_covariance_1Darray(NELEMENTS,var, cov)
et = time.time()
print("# COMPUTATION TIME: {0} sec".format(bt-et))

var.reshape((10,NELEMENTS))
print("Check 10th elements ", var.shape)
print(var[0][10], var[1][10], var[2][10], var[3][10], var[4][10])
cov.reshape((10,10))
print(cov.shape)
# print results
for i in range(10):
    for j in range(i):
        print("cov[{0:2d}][{1:2d}]={3:16.8f}".format(i, j, cov[i][j]))
