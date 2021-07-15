
from plotting_aid import plot_axes_from_evecs_2D

import os
import csv
from matplotlib import pyplot as plt
stylesheet = "../../../default_stylesheet.mplstyle"
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

def plot_eigs(file_name):
    f = open(file_name, 'r')

    M = []
    M += [float(f.readline().split()[-1])]
    M += [float(f.readline().split()[-1])]
    evecs1 = [float(x) for x in f.readline().split()[1:]]
    evecs2 = [float(x) for x in f.readline()[-19:].split()]

    a = 10 ; b = 5.0
    b = a*(M[1]/M[0])**(1/2)
    axes = plot_axes_from_evecs_2D([evecs1, evecs2], [a, b])
    for ax in axes:
        plt.plot(ax[0], ax[1], color="black")

    plt.gca().set_aspect('equal', 'box')

    f.close()



for i in range(40):
    coords_path = "data/coords_%i.csv"%i
    eigen_path = "data/eigenresults_%i.csv"%i

    if not os.path.exists(coords_path) or not os.path.exists(eigen_path):
        continue

    plot_coords(coords_path)
    plot_eigs(eigen_path)
    plt.savefig("images/iteration_%i.png"%i)
    plt.cla()


