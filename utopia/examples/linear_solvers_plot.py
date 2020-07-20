import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv

with open('linear_solvers.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))

print(data)    

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

# print("name: ", name)
# print("------------------ \n")
# print("norm: ", norm)
# print("------------------ \n")
# print("time: ", time)


x = np.arange(len(name))  # the label locations
x2 = np.arange(len(name)) 
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, norm, width, label='Men')
fig2, ax2 = plt.subplots()
rects2 = ax2.bar(x2 + width/2, time, width, label='Women')

# Add some text for name, title and custom x-axis tick name, etc.
ax.set_ylabel('Scores')
ax.set_title('Scores by group and gender')
ax.set_xticks(x)
ax.set_xticklabels(name)
ax2.set_xticklabels(name)
ax.legend()


# def autolabel(rect):
#     """Attach a text label above each bar in *rects*, displaying its height."""
#     for rect in rects:
#         height = rect.get_height()
#         ax.annotate('{}'.format(height),
#                     xy=(rect.get_x() + rect.get_width() / 2, height),
#                     xytext=(0, 3),  # 3 points vertical offset
#                     textcoords="offset points",
#                     ha='center', va='bottom')

# autolabel(rects1)
# autolabel(rects2)

fig.tight_layout()
fig2.tight_layout()

plt.show()

# ------------------------------





