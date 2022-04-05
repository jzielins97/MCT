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
calculate_covariance = libcov.calculate_covariance
calculate_covariance.restype = None
calculate_covariance.argtypes = [c_long, ndpointer(dtype=np.double, ndim=2, shape=(10,NELEMENTS)), ndpointer(dtype=c_double, ndim=2, shape=(10,10), flags='C_CONTIGUOUS')]

read_data_test = libcov.read_data_test
read_data_test.restype = None
read_data_test.argtypes = [c_long,ndpointer(dtype=np.double, ndim=1, shape=(NELEMENTS))]

cov = np.zeros(shape=(10,10),dtype=np.double)
print("cov shape= {0}, is contiguous = {1}".format(cov.shape, cov.flags['C_CONTIGUOUS']))

var = np.empty(shape=(10,NELEMENTS), dtype=np.double)
fileName = "/home2/archive/MCT-2022/lab1/var{0}.dat"
# read data
for i in range(5):
    var[i] = np.fromfile(fileName.format(i+1), dtype=np.double)
    print("File {0} read {1} elements: [10] = {2}".format(i,var[i].shape, var[i][10]))

print("var shape= {0}, is contiguous = {1}".format(var.shape, var.flags['C_CONTIGUOUS']))
read_data_test(NELEMENTS, var[0]);
read_data_test(NELEMENTS, var[1]);
read_data_test(NELEMENTS, var[2]);
read_data_test(NELEMENTS, var[3]);
read_data_test(NELEMENTS, var[4]);
read_data_test(NELEMENTS, var[5]);
read_data_test(NELEMENTS, var[6]);
read_data_test(NELEMENTS, var[7]);
read_data_test(NELEMENTS, var[8]);
read_data_test(NELEMENTS, var[9]);
# calculate covariance
bt = time.time()
calculate_covariance(NELEMENTS, np.ascontiguousarray(var, dtype=np.double), cov)
et = time.time()
print("# COMPUTATION TIME: {0} sec".format(bt-et))

# print results
for i in range(10):
    for j in range(i):
        print("cov[{0:2d}][{1:2d}]={3:16.8f}")

        
        








 
