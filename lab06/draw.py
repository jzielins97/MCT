import matplotlib.pyplot as plt
import csv

X = [] # number of threads
Y = [] # time in seconds

with open('time.txt','r') as file:
    data = csv.reader(file,delimiter='\t') # read data from the file wiht '\t' separating columns
    
    for ROWS in data:
        X.append(int(ROWS[0])) # add thread value to the X list
        Y.append(float(ROWS[1])) # add time value to the Y list

# sort data in terms of threads in case it wasn't saved sorted
X, Y = (list(x) for x in zip(*sorted(zip(X,Y), key=lambda pair: pair[0])))

plt.plot(X,Y,'o-') # plot data as a graf with cirles as markers connected by a line
plt.title("Computational time as function of number of threads used") # set title for the graph
plt.xlabel("Number of threads [-]") # set x label title
plt.ylabel("Computational time [s]") # set y label title
plt.savefig("plot.png", transparent=False, facecolor="white") # save plot as png with white background
