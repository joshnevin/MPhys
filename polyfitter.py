# readDataBasic.py
# program to read in data from a file
# and simply display it
# This is an edit

import matplotlib.pyplot as plt
import numpy as np 
import math 

# store data in these list variables
# create empty lists (arrays) x and y
x = [];
y = [];

# open the file for reading
readFile = open('villoisgraph.csv','r');

# read in data from the file
# line by line
for line in readFile:
    # variable line now holds one of the lines in the
    # data file

    # split up the string 'line' based on whitespace
    splitUp = line.split(',');
        
    # x = splitUp[0], y = splitUp[1]
    # append these values to the arrays x and y
    x.append(splitUp[0]);
    y.append(splitUp[1]);

del x[0] #delete header for column 
del y[0]
del x[len(x)-1] # delete last point - after the cut-off
del y[len(y)-1]


for i in range(len(y)):
    
    y[i] = y[i].strip('\n')
    y[i] = float(y[i])
    
for j in range(len(x)):
  
    x[j] = float(x[j])    
 
readFile.close();

degree = 10

#fit = np.polyfit(x,y,3,full=True) 
fit = np.polyfit(x, y, degree, full=True) 
z = np.linspace(-1.5,1.5,num=100)
print(fit)
p = np.poly1d(np.polyfit(x, y, degree))

plt.plot(x,y)  # plots the data and the polynomial fit 
plt.plot(z,p(z))
plt.xticks(np.arange(-1.5,2,0.5))
plt.show()




print(fit) 


