import matplotlib.pyplot as plt
import numpy as np

def get_averages(type):
    data = []
    for i in range(1,6):
        file = open("Data/"+type+str(i)+".csv", "r")
        file.readline()
        line = file.readline()
        j = 0
        print(len(data))
        while line != "":
            line = line.split(',')
            #print(line)
            if i == 1:
                data.append(float(line[1]))
            else:
                data[j] += float(line[1])
            j += 1
            line = file.readline()

    for i in range(len(data)):
        data[i] /= 5
    return data

def plot_data(points, title):
    plt.plot(points, linewidth=.5)
    plt.title('Opt Cost Over Iterations of ' + title+ ' Algorithm')
    plt.ylabel('Cost')
    plt.xlabel('Generations')
    plt.legend(loc="best")
    plt.savefig(title+"NLGraph")
    plt.show()

if __name__ == "__main__":
    data = get_averages("local_custom_data")
    print("Average final cost is {}".format(data[-1]))
    plot_data(data, "Localized Custom")

