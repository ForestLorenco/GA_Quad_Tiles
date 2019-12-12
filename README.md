# Quadrilateral Tile Image Matching with Genetic Algorithms

This is the code for our final project in nonlinear optimization. The base of the code and original idea was contributed by our professor **Bob Bosch** (http://www.dominoartwork.com/). We have modified it and expanded the crossover methods for 
testing.

This code will attempt to create a image similar to a given image in data file format out of black and white quadrilateral tiles.

## Crossover Method
We developed a crossover method that employs three major parts. 

### Intermediate Values
First we employ an intermediate value crossover tecnique to get the values for the children. Given two parents x and y where (x < y), the childs value at i would be:

c_i = x_i + a(y_i-x_i)

Where a is in [0,1]. This is similar to a convex combination.

### Extinction
Next we employed an extinction threshold. We noticed that the fitness of the children would tend to stagnate. If it stagnated for over 1/4 of total training time, it didn't find new minimums. We added an extinction threshold of 1/4 of total runtime. When this threshold is hit, the entire population is randomly reset, but the old opt value is kept. This has been shown to help find lower mins in practice.


### Localization
The last piece we added was a localization component. Instead of crossing over a large piece from the entire original picture, we seperate the picture into 4 quadrants and crossover 4 small pieces local to those quadrants. While there is a risk of having the crossover be to localized, this has shown better results, especially wheb considerung images with more tiles than the baseline of 22x15.

## How to Run
If you want to try it on your own I made a makefile. Use cd to navigate to the folder. Just go to terminal, cd into the folder by typing:

cd <name of folder>

into terminal and then type

make

into terminal. You will get a lot of warnings from the rans.c file but dont worry. Then simply type

./ga_quads

and it will run for 100000 iterations. This will take about 20 to 30 minutes. Once done, go to file explorer and double click on ga_quads.eps in the folder to see final image. 

## Image Examples
Before and After of custom algorithm


<img src="images/frankenstein300x440.jpg" alt="Kitten"
	title="Before" width="200" height="300" />

<img src="images/ga_quads_local_costum_1M_Ext.jpg" alt="Kitten"
	title="Before" width="200" height="300" />
