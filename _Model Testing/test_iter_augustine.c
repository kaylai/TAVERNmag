#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>

//All compositions are given in mol%
//The below values are given in the order: EA-1 equilibrium fluid; 97009 equilibrium fluid; DVDP 3-295 equilibrium fluid; Lava Lake euqilibrium fluid; Bas-to-Teph MIs; Tehp-to-Phon MIs; Phon-to-LL MIs
float co2_vals[] = {19.93/100.0, 3.83/100.0, 4.68/100.0, 0.04/100.0, 0.06/100.0, 0.19/100.0};
float h2o_vals[] = {71.41/100.0, 24.43/100.0, 35.03/100.0, 98.36/100.0, 27.43/100.0, 90.97/100.0};
float so2_vals[] = {8.67/100.0, 71.74/100.0, 60.29/100.0, 1.60/100.0, 72.51/100.0, 8.84/100.0};

//The below values are all for the surface gas, given in the order: x; x-1; x+1, where x = measured value
float surface_gas_co2[] = {1.30, 3.30, 5.30};
float surface_gas_h2o[] = {90.38, 92.38, 94.38};
float surface_gas_so2[] = {2.32, 4.32, 6.32};

//Defines an inline function that multiplies the possible percentage values a,b,c,d,e,f, and g by the concentration of some species (CO2, H2O, or Stot) in each contributing fluid
float getVals(int a, int b, int c, int d, int e, int f, float *ptr){
  return (a*ptr[0])+(b*ptr[1])+(c*ptr[2])+(d*ptr[3])+(e*ptr[4])+(f*ptr[5]);
}

//Define the variables to be used
int main(void){
  register int a,b,c,d,e,f;
  double x=0;
  float val;
  int i;
  int start_val, end_val;
  char s[200];
  FILE *ofp;
  
  pid_t pid;
  
	//Generates a list of one-dimensional arrays, each array 7 values long, where the sum of each array is 100
  for (i = 0; i < 4; i++){
    if ((pid = fork()) == 0){
      sprintf(s,"iter_out_%d",i);
      ofp = fopen(s,"w");
      start_val = i*25;
      end_val = start_val + 26;
      for (a=start_val;a<end_val;a++){
        printf("proc%d: i = %d\n",i,a);
        for (b=100-a; b>-1; b--){
          for (c=100-a-b; c>-1; c--){
            for (d=100-a-b-c; d>-1; d--){
              for (e=100-a-b-c-d; e>-1; e--){
                for (f=100-a-b-c-d-e; f>-1; f--){
                    
					  //If G - 1 >= the sum of the proportioned fluid contributions >= G+1, exclude that array from the final list. G = measured surface gas value of each species.
                    val = getVals(a,b,c,d,e,f,co2_vals);
                    if ((val < surface_gas_co2[1]) || (val >surface_gas_co2[2])){
                      continue;
                    }
                    
                    val = getVals(a,b,c,d,e,f,h2o_vals);
                    if ((val < surface_gas_h2o[1]) || (val >surface_gas_h2o[2])){
                      continue;
                    }
                    
                    val = getVals(a,b,c,d,e,f,so2_vals);
                    if ((val < surface_gas_so2[1]) || (val >surface_gas_so2[2])){
                      continue;
                    }
                    
                    //Prints all arrays where G - 1 <= sum of the proportioned fluid contributions <= G + 1 is true for all gas species. Output is saved to four files: one file for each forked process.
                    fprintf(ofp,"%d %d %d %d %d %d\n",a,b,c,d,e,f);
                        
                  }
                }
              }
            }
          }
        }
      
      fclose(ofp);
    }
  }
  wait(NULL);
  wait(NULL);
  wait(NULL);
  wait(NULL);
  return 0;
 }
      
    
  
