import matplotlib.pyplot as plt
import csv

X = [] # number of threads
t_read = [] # time in seconds
bw_read = [] # bandwidth [MB/s]
t_write = []
bw_write = []

mpi_read = [] # read time for the send/recv routine

with open('stats_mpiio.txt','r') as file:
    data = csv.reader(file,delimiter='\t') # read data from the file wiht '\t' separating columns
    next(data) # skipping first line in the file
    for ROWS in data:
        X.append(int(ROWS[0])) # add thread value to the X list
        t_read.append(float(ROWS[1]))
        bw_read.append(float(ROWS[2]))
        t_write.append(float(ROWS[4]))
        bw_write.append(float(ROWS[5]))

print(X)
print(t_read)
print(t_write)
with open('stats_mpi.txt','r') as file:
    data = csv.reader(file,delimiter='\t') # read data from the file wiht '\t' separating columns
    next(data) # skipping first line in the file
    for ROWS in data:
        mpi_read.append(float(ROWS[1]))
print(mpi_read)

# sort data in terms of threads in case it wasn't saved sorted
# X, Y = (list(x) for x in zip(*sorted(zip(X,Y), key=lambda pair: pair[0])))


plt.plot(X,t_read,'o-',color='r')
plt.plot(X,t_write,'o-',color='g')
plt.plot(X,mpi_read,'o-', color='b')
plt.title("Computational time as function of number of threads used") # set title for the graph
plt.legend(['MPII/O read','MPII/O write','MPI send/recv'])
plt.title("Comparison of read/write times for the MPII/O and send/recv") # set title for the graph
plt.xlabel("Number of threads [-]") # set x label title
plt.ylabel("Operation time [s]") # set y label title
plt.savefig("times.png", transparent=False, facecolor="white") # save plot as png with white background


