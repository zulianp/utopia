import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv
import os



with open(path) as csvfile:
    data = list(csv.reader(csvfile))
   
norm = []
time = []
name = []

for i in range(0, len(data)-1):
	for j in range (0, 3):
		if j == 0:
			name.append(data[i][j])
		elif j == 1:
			norm.append(float(data[i][j]))
		elif j == 2:
			time.append(float(data[i][j]))


x = np.arange(len(name))  # the label locations
x2 = np.arange(len(name)) 
width = 0.35  # the width of the bars

fig, ax = plt.subplots()

fig2, ax2 = plt.subplots()


# Add some text for name, title and custom x-axis tick name, etc.
# ax.set_ylabel('')
ax.set_title('Norm')
ax.set_xticks(x)
ax.set_xticklabels(name)
ax.legend()



fig.tight_layout()
fig2.tight_layout()

plt.show()

# # ------------------------------




