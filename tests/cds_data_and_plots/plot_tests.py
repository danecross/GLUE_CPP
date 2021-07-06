

import csv
from matplotlib import pyplot as plt
stylesheet = "../../../../default_stylesheet.mplstyle"
plt.style.use(stylesheet)

def plot_coords(file_name):
    x = [] ; y = []
    with open(file_name, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            x += [float(row[0])] ; y += [float(row[1])]



    plt.gcf().set_size_inches(20,10)
    plt.gca().set_aspect('equal', 'box')
    plt.plot(x, y, '.')


###############
## CONST BOX ##
###############

plot_coords('const_box_coords.csv')
plt.savefig("const_box.png")
plt.cla()

#################
## CONST SHELL ##
#################

plot_coords('const_shell_coords.csv')
plt.savefig("const_shell.png")
plt.cla()


