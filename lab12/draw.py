import matplotlib.pyplot as plt
import numpy as np
import csv

# number of threads 
procs = [] # no of processes
nodes = [1, 2, 3, 4, 5, 6, 7, 8] # number of nodes (one node is 24 processes)
# time in seconds
invert_t = [] # time of invertion of the matrix
mult_t = [] # multiplication time

# get data from openMP
with open('time.csv','r') as file:
    data = csv.DictReader(file,delimiter='\t') # read data from the file wiht '\t' separating columns
    for ROWS in data:
        procs.append(int(ROWS['No of processes'])) # add number of processes to the procs list
        invert_t.append(float(ROWS['Invertion time'])) # add time value to the invert_t list
        mult_t.append(float(ROWS['Multiplication time'])) # add time value to the mult_t list


fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

# Add some extra space for the second axis at the bottom
fig.subplots_adjust(bottom=0.15)

ax1.plot(procs,invert_t,'o-',color='r') # plot data as a graf with cirles as markers connected by a line
ax1.plot(procs,mult_t, '*-',color='b')

ax1.legend(['Invertion time','Multiplication time'])
plt.title("Computational time as function of number of processes used for a 10000x10000 matrix") # set title for the graph
ax1.set_xlabel("Number of processes [-]") # set x label title
ax1.set_ylabel("Computational time [s]") # set y label title


# Move twinned axis ticks and label from top to bottom
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")

# Offset the twin axis below the host
ax2.spines["bottom"].set_position(("axes", -0.1))

ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(np.linspace(24,192,8))
ax2.set_xticklabels(nodes)
ax2.set_xlabel("Number of nodes")
#ax2.spines["bottom"].set_visible(True)
plt.savefig("plot.png", transparent=False, facecolor="white") # save plot as png with white background
