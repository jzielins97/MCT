import matplotlib.pyplot as plt
import csv

# number of threads 
X1 = [] # openMP
X2 = [] # MPI
X3 = [] # theory
# time in seconds
Y1 = [] # openMP
Y2 = [] # MPI
Y3 = [] # theory

# get data from openMP
with open('../lab02/time.txt','r') as file:
    data = csv.reader(file,delimiter='\t') # read data from the file wiht '\t' separating columns
    
    for ROWS in data:
        X1.append(int(ROWS[0])) # add thread value to the X list
        Y1.append(float(ROWS[1])) # add time value to the Y list


# get data from MPI
with open('../lab03/time_mpi.txt','r') as file:
    data = csv.reader(file,delimiter='\t') # read data from the file wiht '\t' separating columns
    
    for ROWS in data:
        X2.append(int(ROWS[0])) # add thread value to the X list
        Y2.append(float(ROWS[1])) # add time value to the Y list


NP = 100
T1 = X2[0] # (X1[0] + X2[0])/2
for i in range(NP):
    X3.append(int(i))
    Y3.append(float(T1 / (i+1)))
    
# sort data in terms of threads in case it wasn't saved sorted
X1, Y1 = (list(x) for x in zip(*sorted(zip(X1,Y1), key=lambda pair: pair[0])))
X2, Y2 = (list(x) for x in zip(*sorted(zip(X2,Y2), key=lambda pair: pair[0])))
X3, Y3 = (list(x) for x in zip(*sorted(zip(X3,Y3), key=lambda pair: pair[0])))

plt.plot(X1,Y1,'o-',color='r') # plot data as a graf with cirles as markers connected by a line
plt.plot(X2,Y2, 'o-',color='b')
plt.plot(X3,Y3,'--',color='g')

plt.legend(['openMP','MPI','ideal scaling'])
plt.title("Computational time as function of number of threads used") # set title for the graph
plt.xlabel("Number of threads [-]") # set x label title
plt.ylabel("Computational time [s]") # set y label title
plt.savefig("plot.png", transparent=False, facecolor="white") # save plot as png with white background
