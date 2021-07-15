
import ellipse_functions as ef
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

def plot_ellipses(file_name):
    f = open(file_name, 'r')
    
    lower_ellipse = [float(x) for x in f.readline().split()]
    upper_ellipse = [float(x) for x in f.readline().split()]

    th = ef.FULL_THETA ; alpha = upper_ellipse[2]
    lower_r = ef.r(th, lower_ellipse[0], lower_ellipse[1])
    upper_r = ef.r(th, upper_ellipse[0], upper_ellipse[1], alpha=alpha)

    x_low, y_low = ef.polar_to_cartesian(lower_r, th)
    x_up, y_up = ef.polar_to_cartesian(upper_r, th)

    plt.plot(x_low, y_low, label="lower ellipse")
    plt.plot(x_up, y_up, label="upper fit ellipse")

    f.close()

###############
## CONST BOX ##
###############

plot_coords('data/const_box_coords.csv')
plt.savefig("images/const_box.png")
plt.cla()

#################
## CONST SHELL ##
#################

plot_coords('data/const_shell_coords.csv')
plt.savefig("images/const_shell.png")
plt.cla()

#############
## CDS FIT ##
#############

plot_coords("data/coords_final.csv")
plot_ellipses("data/ellipses_final.csv")
plt.savefig("images/cds_final_results.png")
plt.cla()




