/*
   This program reads the files

      picture.dat

   and produces the file

      ga_quads.eps
 */

/*
Struct definition for dark and light structures
*/
struct Point{
  int x, y;
};
//Declare point as struct p1 = {x,y};

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ROWS 60
#define MAX_COLS 60

#define LOW 0.25 //ask about these
#define HIGH 0.75

#define ACCENTUATE 1  /* If this = 1, the b[i][j]'s are accentuated: */
#define WEIGHT 3      /*   b[i][j] = 0.5 + WEIGHT*(b[i][j] - 0.5)    */
#define ENCOURAGE 1   /* If this = 1, certain mutations are encouraged */

/* GA parameters */

#define POP_SIZE 100
#define MAX_IT 100000

#define TRUE 1
#define FALSE 0

#define PROB_CROSSOVER 0.00
#define PROB_MUTATION 0.01
#define COMPRESS 0.5

double h[POP_SIZE][MAX_ROWS+1][MAX_COLS+1];
double v[POP_SIZE][MAX_ROWS+1][MAX_COLS+1];

double h_temp[POP_SIZE][MAX_ROWS+1][MAX_COLS+1];
double v_temp[POP_SIZE][MAX_ROWS+1][MAX_COLS+1];

double best_h[MAX_ROWS+1][MAX_COLS+1];
double best_v[MAX_ROWS+1][MAX_COLS+1];

double opt_h[MAX_ROWS+1][MAX_COLS+1];
double opt_v[MAX_ROWS+1][MAX_COLS+1];

double cost_h[MAX_ROWS+1][MAX_COLS+1];
double cost_v[MAX_ROWS+1][MAX_COLS+1];

double fitness[POP_SIZE];

double cdf[POP_SIZE];

//Define variables for custom crossover
int d_len;
struct Point dark[MAX_COLS*MAX_ROWS];
int l_len;
struct Point light[MAX_ROWS*MAX_COLS];

/* Random number generation parameters */

double ran_result, rans_(); // rans returns number in  [0,1)

void random_();

/* Postscript parameters */

#define LINEWIDTH 0.05
#define BW 10.0 /* border width */ 
#define PSH 500.0 /* postscript height */
#define PSW 500.0 /* postscript width */

#define ZERO 0.1

#define RGRID 0.5
#define GGRID 0.5
#define BGRID 0.5

#define BLACK 0.0
#define WHITE 1.0

#define JUST_BORDER 0

/* Global variables */

int r, c, ch, best_ch, worst_ch, opt_ch, row, col, i, j, t;
int rows, cols, low_r, high_r, low_c, high_c, pair;
int it, counter, parent1, parent2, temp_int;

double b[MAX_ROWS][MAX_COLS];  

double north, south, east, west;
double tri_area, cost, best_fitness, worst_fitness, range, sum_of_fitnesses;
double opt_fitness;
double sum_of_costs, average_of_costs;
double prob;

double old_partial_cost, new_partial_cost, best_partial_cost;

double max_vert, max_horiz;
double min_vert, min_horiz;

double height, width;
double vert_scale, horiz_scale, scale;

double shade, rshade, gshade, bshade;
double error;

float result;

char trash_string[25];

FILE *picture_file, *ps;

/*
====================
=  draw postcript  =
====================
*/

void draw_postscript()
{

/*
   Create the postscript file
 */

if ((ps=fopen("ga_quads.eps", "w")) == NULL)
  {
    printf("Couldn't open ga_quads.ps!\n");
    exit(0);
  }

fprintf(ps, "%%!PS-Adobe-3.0 EPSF-3.0\n"); 
fprintf(ps, "%%%%BoundingBox: 0 0 %d %d\n\n", (int) PSW, (int) PSH); 

fprintf(ps, "1 setlinecap\n"); 
fprintf(ps, "1 setlinejoin\n");
fprintf(ps, "%2f setlinewidth\n",scale*LINEWIDTH); 

fprintf(ps, "0 %.2f translate\n\n", (PSH-scale*height-2.0*BW)/2.0);

for (r = 1; r <= rows; r++)
  for (c = 1; c <= cols; c++)
    {
      north = opt_h[r-1][c];
      south = opt_h[r][c];
      east = opt_v[r][c];
      west = opt_v[r][c-1];
      if (((r+c) % 2) == 0)
        {
          if (JUST_BORDER == 0)
            {
              shade = BLACK;
              fprintf(ps, "%.3f %.3f %.3f setrgbcolor\n\n", shade, shade, shade);

              fprintf(ps, "newpath\n");
              fprintf(ps, "%.6f %.6f moveto\n",
                scale*(c-1 + north - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 0.0)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + 1.0 - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + east)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + south - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 1.0)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + 0.0 - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + west)) + BW);
              fprintf(ps, "closepath\n");
              fprintf(ps, "fill\n\n");
            }
          fprintf(ps, "%.3f %.3f %.3f setrgbcolor\n\n", RGRID, GGRID, BGRID);
          fprintf(ps, "newpath\n");
          fprintf(ps, "%.6f %.6f moveto\n",
            scale*(c-1 + north - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 0.0)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + 1.0 - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + east)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + south - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 1.0)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + 0.0 - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + west)) + BW);
          fprintf(ps, "closepath\n");
          fprintf(ps, "stroke\n\n");
        }
      else
        {
          if (JUST_BORDER == 0)
            {
              shade = BLACK;
              fprintf(ps, "%.3f %.3f %.3f setrgbcolor\n\n", shade, shade, shade);

              fprintf(ps, "newpath\n");
              fprintf(ps, "%.6f %.6f moveto\n",
                scale*(c-1 + 0.0 - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 0.0)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + north  - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 0.0)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + 0.0 - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + west)) + BW);
              fprintf(ps, "closepath\n");
              fprintf(ps, "fill\n\n");

              fprintf(ps, "newpath\n");
              fprintf(ps, "%.6f %.6f moveto\n",
                scale*(c-1 + 1.0 - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 0.0)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + 1.0  - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + east)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + north - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 0.0)) + BW);
              fprintf(ps, "closepath\n");
              fprintf(ps, "fill\n\n");

              fprintf(ps, "newpath\n");
              fprintf(ps, "%.6f %.6f moveto\n",
                scale*(c-1 + 1.0 - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 1.0)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + south  - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 1.0)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + 1.0 - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + east)) + BW);
              fprintf(ps, "closepath\n");
              fprintf(ps, "fill\n\n");

              fprintf(ps, "newpath\n");
              fprintf(ps, "%.6f %.6f moveto\n",
                scale*(c-1 + 0.0 - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 1.0)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + 0.0  - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + west)) + BW);
              fprintf(ps, "%.6f %.6f lineto\n",
                scale*(c-1 + south - min_horiz) + (PSW-scale*width)/2.0, 
                scale*(max_vert - (r-1 + 1.0)) + BW);
              fprintf(ps, "closepath\n");
              fprintf(ps, "fill\n\n");
            }
          fprintf(ps, "%.3f %.3f %.3f setrgbcolor\n\n", RGRID, GGRID, BGRID);
          fprintf(ps, "newpath\n");

          fprintf(ps, "%.6f %.6f moveto\n",
            scale*(c-1 + 0.0 - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 0.0)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + north  - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 0.0)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + 0.0 - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + west)) + BW);
          fprintf(ps, "closepath\n");
          fprintf(ps, "stroke\n\n");

          fprintf(ps, "newpath\n");
          fprintf(ps, "%.6f %.6f moveto\n",
            scale*(c-1 + 1.0 - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 0.0)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + 1.0  - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + east)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + north - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 0.0)) + BW);
          fprintf(ps, "closepath\n");
          fprintf(ps, "stroke\n\n");

          fprintf(ps, "newpath\n");
          fprintf(ps, "%.6f %.6f moveto\n",
            scale*(c-1 + 1.0 - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 1.0)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + south  - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 1.0)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + 1.0 - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + east)) + BW);
          fprintf(ps, "closepath\n");
          fprintf(ps, "stroke\n\n");

          fprintf(ps, "newpath\n");
          fprintf(ps, "%.6f %.6f moveto\n",
            scale*(c-1 + 0.0 - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 1.0)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + 0.0  - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + west)) + BW);
          fprintf(ps, "%.6f %.6f lineto\n",
            scale*(c-1 + south - min_horiz) + (PSW-scale*width)/2.0, 
            scale*(max_vert - (r-1 + 1.0)) + BW);
          fprintf(ps, "closepath\n");
          fprintf(ps, "stroke\n\n");
        }
    }
fprintf(ps, "showpage\n");

fclose(ps);
}

/*
=====================
= calculate fitness =
=====================
*/
void calc_fitness(){
  
  for (r = 1; r <= rows; r++)//loop to find the the cost of chromosome ch
    for (c = 1; c <= cols; c++)
      {
        /*
          *Get the four sides of the quadrilateral
        */
        north = h[ch][r-1][c];
        south = h[ch][r][c];
        east = v[ch][r][c];
        west = v[ch][r][c-1];
        tri_area = 0.5*(1.0 + north*west - north*east + south*east - south*west); //Get area of trianges
        if (((r+c)%2) == 0)//we are switching tile colors every other
          cost += fabs(tri_area - b[r-1][c-1]);
        else
          cost += fabs((1.0 - tri_area) - b[r-1][c-1]);
      }
    
}

/*
=====================
=     crossover     =
=====================
*/

//Gets the info from the picture
void get_pic_info(){
  float threshold = 0; //threshold for brightness
  d_len = 0;
  l_len = 0;
  //b has picture data stored in it
  for(r = 0; r < rows; r++)
    for(c = 0; c < cols; c++){
        if(b[r][c] < threshold){
            dark[d_len].x = r+1;
            dark[d_len].y = c+1;
            d_len ++;
        }else{
            light[l_len].x = r+1;
            light[l_len].y = c+1;
            l_len ++;
        }
    }
}
//our custom crossover algorithm
void crossover_custom(){

}


//Bobs original crossover algorithm
void crossover_bob(){
  for (pair = 0; pair < POP_SIZE/2; pair++) //Loop through pairs of ch
      {
        ran_result = rans_(); // get random number
        counter = 0;
        while ((counter < POP_SIZE) && (ran_result > cdf[counter]))//loop through ch until r < cdf
          counter++;
        parent1 = counter; //set parent one to be that ch
        ran_result = rans_(); // Get new random number
        counter = 0;
        while ((counter < POP_SIZE) && (ran_result > cdf[counter]))
          counter++;
        parent2 = counter;//get second parent

        ran_result = rans_(); // get new random number
        if (ran_result < PROB_CROSSOVER) // if we are within probality threshold
          {
            ran_result = rans_();//get ran number
            low_r = 1 + floor(rows * ran_result);//get random low row 		
            ran_result = rans_();
            high_r = 1 + floor(rows * ran_result); //get random high row	
            if (low_r > high_r) // if low > high swap them
              {
                temp_int = low_r;
                low_r = high_r;
                high_r = temp_int;
              }
            //repeat to get low column and high column
            ran_result = rans_(); 
            low_c = 1 + floor(cols * ran_result);		
            ran_result = rans_();
            high_c = 1 + floor(cols * ran_result);		
            if (low_c > high_c)
              {
                temp_int = low_c;
                low_c = high_c;
                high_c = temp_int;
              }
            //copy the parents into the new genome as pairs
            //pair[0] = parent1
            for (r = 0; r <= rows; r++)
              for (c = 1; c <= cols; c++)
                h[2*pair][r][c] = h_temp[parent1][r][c];
            for (r = 1; r <= rows; r++)
              for (c = 0; c <= cols; c++)
                v[2*pair][r][c] = v_temp[parent1][r][c];

            //pair[1] = parent2
            for (r = 0; r <= rows; r++)
              for (c = 1; c <= cols; c++)
                h[2*pair+1][r][c] = h_temp[parent2][r][c];
            for (r = 1; r <= rows; r++)
              for (c = 0; c <= cols; c++)
                v[2*pair+1][r][c] = v_temp[parent2][r][c];
            //loop through rows and columns
            for (r = 1; r <= rows; r++)
              for (c = 1; c <= cols; c++)
              //if lowr <= r <= highr (same for c)
                if ((r >= low_r) && (r <= high_r) && (c >= low_c) && (c <= high_c))
                  {
                    //set pair[0] to have parent2s info
                    h[2*pair][r][c] = h_temp[parent2][r][c];
                    h[2*pair][r-1][c] = h_temp[parent2][r-1][c];
                    v[2*pair][r][c] = h_temp[parent2][r][c];
                    v[2*pair][r][c-1] = h_temp[parent2][r][c-1];

                    //set pair[1] to have parent1s info
                    h[2*pair+1][r][c] = h_temp[parent1][r][c];
                    h[2*pair+1][r-1][c] = h_temp[parent1][r-1][c];
                    v[2*pair+1][r][c] = h_temp[parent1][r][c];
                    v[2*pair+1][r][c-1] = h_temp[parent1][r][c-1];
                  }
          }
        //if no crossover, just copy parents into new generation
        else
          {
            //pair[0] = parent1
            for (r = 0; r <= rows; r++)
              for (c = 1; c <= cols; c++)
                h[2*pair][r][c] = h_temp[parent1][r][c];
            for (r = 1; r <= rows; r++)
              for (c = 0; c <= cols; c++)
                v[2*pair][r][c] = v_temp[parent1][r][c];

            //pair[1] = parent2
            for (r = 0; r <= rows; r++)
              for (c = 1; c <= cols; c++)
                h[2*pair+1][r][c] = h_temp[parent2][r][c];
            for (r = 1; r <= rows; r++)
              for (c = 0; c <= cols; c++)
                v[2*pair+1][r][c] = v_temp[parent2][r][c];
          }
      }
}

/*
=====================
=      mutation     =
=====================
*/
void mutation(){
  //Mutations for h
  for (r = 0; r <= rows; r++)
    for (c = 1; c <= cols; c++)
      {
        ran_result = rans_(); //get random num
        //if encourage, increase prob
        if (ENCOURAGE == 1)
          prob = PROB_MUTATION*cost_h[r][c]/average_of_costs;
        else
          prob = PROB_MUTATION;

        //if we are in proportion
        if (ran_result < prob)
          {
            ran_result = rans_(); //get random num
            ran_result = COMPRESS*(ran_result - 0.5); //compress result
            h[ch][r][c] += ran_result; //add it to h(r,c)
            //scale it within acceptable bounds
            if (h[ch][r][c] < LOW)
              h[ch][r][c] = LOW;
            if (h[ch][r][c] > HIGH)
              h[ch][r][c] = HIGH;
          }
      }
  //Mutations for v, same as above
  for (r = 1; r <= rows; r++)
    for (c = 0; c <= cols; c++)
      {
        ran_result = rans_(); 
        if (ENCOURAGE == 1) 
          prob = PROB_MUTATION*cost_v[r][c]/average_of_costs;
        else
          prob = PROB_MUTATION;
        if (ran_result < prob)
          {
            ran_result = rans_(); 
            ran_result = COMPRESS*(ran_result - 0.5);
            v[ch][r][c] += ran_result;
            if (v[ch][r][c] < LOW)
              v[ch][r][c] = LOW;
            if (v[ch][r][c] > HIGH)
              v[ch][r][c] = HIGH;
          }
      }
}

/*
========================
=        main          =
========================
*/
int main ()
{

/* 
    Read picture.dat
 */

if ((picture_file = fopen("picture_frank22x15.dat", "r")) == NULL)
  {
    printf("Couldn't open picture.dat!\n");
    exit(0);
  }

fscanf(picture_file, "%s", trash_string);//Load the first line into trash string
fscanf(picture_file, "%d", &rows);//Load the # rows
fscanf(picture_file, "%s", trash_string);//load the next line into trash string
fscanf(picture_file, "%d", &cols);//load the # cols

/*
*Loop through the remaining data to get the values for each box in file
*/
for (r = 0; r < rows; r++)
  for (c = 0; c < cols; c++)
    {
      fscanf(picture_file, "%f", &result); //load (r,c) into result
      if (ACCENTUATE == 1)  //if acc modify result
        b[r][c] = 0.5 + 3*(result-0.5);
      else
        b[r][c] = result;
    }

fclose(picture_file);

get_pic_info();//get picture info for custom crossover

/*
 *  Construct the initial population of chromosomes.
 */

random_();//call random to get seed

for (ch = 0; ch < POP_SIZE; ch++)//loop through population with ch (chromosome?)
  {
    for (r = 0; r <= rows; r++)//loop through rows = 1 and cols = 0
      for (c = 1; c <= cols; c++)//construct the h population
        {
          ran_result = rans_(); //get random number
          h[ch][r][c] = LOW + ran_result*(HIGH-LOW); //I think this keeps it in the sliders range
          if (h[ch][r][c] < LOW)//edge case to hande if we go outside the slider range
            h[ch][r][c] = LOW;
          if (h[ch][r][c] > HIGH)
            h[ch][r][c] = HIGH;
        }
    for (r = 1; r <= rows; r++)
      for (c = 0; c <= cols; c++)//same as above for the v population
        {
          ran_result = rans_();
          v[ch][r][c] = LOW + ran_result*(HIGH-LOW);
          if (v[ch][r][c] < LOW)
            v[ch][r][c] = LOW;
          if (v[ch][r][c] > HIGH)
            v[ch][r][c] = HIGH;
        }
    cost = 0.0; //give initial cost
    calc_fitness();
    fitness[ch] = rows*cols - cost;//save in fitness
  }

/*
 *  Find the chromosomes with the best and worst fitness values.
 */

best_fitness = 0.0;
opt_fitness = 0.0;
opt_ch = -1;
best_ch = -1;

for (ch = 0; ch < POP_SIZE; ch++)//loop through chromosomes in the population
  if (fitness[ch] > best_fitness)//find the best fitness
    {
      best_ch = ch;
      best_fitness = fitness[ch];
    }

for (ch = 0; ch < POP_SIZE; ch++)
  if (fitness[ch] > opt_fitness) //find opt fitness (will be best, redundent)
    {
      opt_ch = ch;
      opt_fitness = fitness[ch];
    }

worst_fitness = rows*cols; // Set worst fitness to maximum area
worst_ch = -1;

for (ch = 0; ch < POP_SIZE; ch++)
  if (fitness[ch] < worst_fitness)//find worst fitness in population
    {
      worst_ch = ch;
      worst_fitness = fitness[ch];
    }

for (ch = 0; ch < POP_SIZE; ch++)//set fitness of ch to fitness-worst
  fitness[ch] = fitness[ch] - worst_fitness;

for (r = 0; r <= rows; r++)
  for (c = 1; c <= cols; c++)
    best_h[r][c] = h[best_ch][r][c]; //set the best h(r,c) to h with best fitness

for (r = 1; r <= rows; r++)
  for (c = 0; c <= cols; c++)
    best_v[r][c] = v[best_ch][r][c]; //set the best v(r,c) to h with best fitness

for (r = 0; r <= rows; r++)
  for (c = 1; c <= cols; c++)
    h[worst_ch][r][c] = h[best_ch][r][c]; //replace the worst h with best h

for (r = 1; r <= rows; r++)
  for (c = 0; c <= cols; c++)
    v[worst_ch][r][c] = v[best_ch][r][c]; //replace the worst v with best v

fitness[worst_ch] = fitness[best_ch];//set fitness of worst to best

if (opt_ch != -1)//if we have found an opt, set opt h and v to opt h and v
  {
    for (r = 0; r <= rows; r++)
      for (c = 1; c <= cols; c++)
        opt_h[r][c] = h[opt_ch][r][c];

    for (r = 1; r <= rows; r++)
      for (c = 0; c <= cols; c++)
        opt_v[r][c] = v[opt_ch][r][c];
  }

sum_of_fitnesses = 0.0;
for (ch = 0; ch < POP_SIZE; ch++)
  sum_of_fitnesses += fitness[ch]; // sum all fitnesses

cdf[0] = fitness[0]/sum_of_fitnesses; //this does something by computing the sum of (fitnesses/sum of fitnesses) up to ch
for (ch = 1; ch < POP_SIZE; ch++)
  cdf[ch] = cdf[ch-1] + fitness[ch]/sum_of_fitnesses;
    
/*
    Compute scaling for postscript file.
 */

min_horiz = 0.0;
max_horiz = 1.0*cols;
min_vert = 0.0;
max_vert = 1.0*rows;

height = max_vert - min_vert;
width = max_horiz - min_horiz;

vert_scale = (PSH-2.0*BW)/height;
horiz_scale = (PSW-2.0*BW)/width;
if (vert_scale < horiz_scale)
  scale = vert_scale;
else
  scale = horiz_scale;

printf(" Initially, avg. cost = %.6f\n", 
        (rows*cols - opt_fitness)/(rows*cols)); //print average cost

draw_postscript(); // draw the image

for (it = 1; it <= MAX_IT; it++) // loop through the max iterations
  {
    for (ch = 0; ch < POP_SIZE; ch++) // loop through each chromosome and save the current population
      {
        for (r = 0; r <= rows; r++)
          for (c = 1; c <= cols; c++)
            h_temp[ch][r][c] = h[ch][r][c];
        for (r = 1; r <= rows; r++)
          for (c = 0; c <= cols; c++)
            v_temp[ch][r][c] = v[ch][r][c];
      }

/*
   Natural selection and (possibly) crossover.
 */
    crossover_bob(); // Bobs crossover method

/*  
   Mutate some parts of some of the chromosomes.
 */

    for (ch = 0; ch < POP_SIZE; ch++)//loop through each ch
      {
        //reset cost h
        for (r = 0; r <= rows; r++)
          for (c = 1; c <= cols; c++) 
            cost_h[r][c] = 0.0; 
        //reset cost v
        for (r = 1; r <= rows; r++)
          for (c = 0; c <= cols; c++)
            cost_v[r][c] = 0.0; 

        //recalculate the cost

        
        for (r = 1; r <= rows; r++)
          for (c = 1; c <= cols; c++)
            {
              north = h[ch][r-1][c];
              south = h[ch][r][c];
              east = v[ch][r][c];
              west = v[ch][r][c-1];
              tri_area = 0.5*(1.0 + north*west - north*east + south*east - south*west);
              if (((r+c)%2) == 0)
                cost += fabs(tri_area - b[r-1][c-1]);
              else
                cost += fabs((1.0 - tri_area) - b[r-1][c-1]);
              cost_h[r-1][c] += cost;
              cost_h[r][c] += cost;
              cost_v[r][c] += cost;
              cost_v[r][c-1] += cost;
            }
        sum_of_costs = 0.0;
        //sum all costs together
        for (r = 0; r <= rows; r++)
          for (c = 1; c <= cols; c++)
            sum_of_costs += cost_h[r][c];
        for (r = 1; r <= rows; r++)
          for (c = 0; c <= cols; c++)
            sum_of_costs += cost_v[r][c];
        //average the costs
        average_of_costs = sum_of_costs/((rows+1)*cols + rows*(cols+1));

        //perform the mutations
        mutation();
      }
    //recalculate fitness again after mutations
    for (ch = 0; ch < POP_SIZE; ch++)
      {
        cost = 0.0;
        calc_fitness();
         fitness[ch] = rows*cols - cost;
      }

    //Get best, worst and opt fitness
    best_fitness = 0.0;
    best_ch = -1;
    opt_ch = -1;

    for (ch = 0; ch < POP_SIZE; ch++)
      if (fitness[ch] > best_fitness)
        {
          best_ch = ch;
          best_fitness = fitness[ch];
        }

    for (ch = 0; ch < POP_SIZE; ch++)
      if (fitness[ch] > opt_fitness)
        {
          opt_ch = ch;
          opt_fitness = fitness[ch];
        }

    worst_fitness = rows*cols;
    worst_ch = -1;

    for (ch = 0; ch < POP_SIZE; ch++)
      if (fitness[ch] < worst_fitness)
        {
          worst_ch = ch;
          worst_fitness = fitness[ch];
        }

    for (ch = 0; ch < POP_SIZE; ch++)
      fitness[ch] = fitness[ch] - worst_fitness;

    for (r = 0; r <= rows; r++)
      for (c = 1; c <= cols; c++)
        best_h[r][c] = h[best_ch][r][c];

    for (r = 1; r <= rows; r++)
      for (c = 0; c <= cols; c++)
        best_v[r][c] = v[best_ch][r][c];

    for (r = 0; r <= rows; r++)
      for (c = 1; c <= cols; c++)
        h[worst_ch][r][c] = h[best_ch][r][c];

    for (r = 1; r <= rows; r++)
      for (c = 0; c <= cols; c++)
        v[worst_ch][r][c] = v[best_ch][r][c];

    fitness[worst_ch] = fitness[best_ch];

    if (opt_ch != -1)
      {
        for (r = 0; r <= rows; r++)
          for (c = 1; c <= cols; c++)
            opt_h[r][c] = h[opt_ch][r][c];

        for (r = 1; r <= rows; r++)
          for (c = 0; c <= cols; c++)
            opt_v[r][c] = v[opt_ch][r][c];
      }

    //recalc sum of fitness
    sum_of_fitnesses = 0.0;
    for (ch = 0; ch < POP_SIZE; ch++)
      sum_of_fitnesses += fitness[ch];

    //recalculate cdf
    cdf[0] = fitness[0]/sum_of_fitnesses;
    for (ch = 1; ch < POP_SIZE; ch++)
      cdf[ch] = cdf[ch-1] + fitness[ch]/sum_of_fitnesses;
    
    //every 500 iterations print results
    if (it % 500 == 0)
      {
        printf(" After %d iterations, avg. cost = %.6f\n", it, 
          (rows*cols - opt_fitness)/(rows*cols));
        draw_postscript();//draw postcript
      }
  }

}

