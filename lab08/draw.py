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


plt.plot(X,bw_read,'o-',color='r')
plt.plot(X,bw_write,'o-',color='g')
plt.legend(['MPII/O read','MPII/O write'])
plt.title("Bandwidth for read and write with MPII/O") # set title for the graph
plt.xlabel("Number of threads [-]") # set x label title
plt.ylabel("Bandwidth [MB/s]") # set y label title
plt.savefig("bandwidth.png", transparent=False, facecolor="white") # save plot as png with white background


