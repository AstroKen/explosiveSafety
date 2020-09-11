#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <bash/bash.h>
#include <array/array.h>

//checkInput
void checkInput(int);
//load file
void loadFile(char*, cmatrix, dvector);
//calc sum of W_{ei}
double sumWei(cmatrix*, dvector*);
//calc sum of W_{eo}
double sumWeo(cmatrix*, dvector*);
//calc sum of weight of propellants
double sumWp(dvector);
//calc Security Distance for Blast
double calcBSD(cmatrix, dvector);
//calc Security Distance for Scattered Object
double calcSSD(cmatrix, dvector);
//calc Security Distance for Fireball;
double calcEisenbergFSD(cmatrix, dvector);

double calcNASAFSD(cmatrix, dvector);

int main(int argc, char *argv[])
{
  char *fileName;
  int num;
  double BSD; // blast security distance
  double SSD;
  double eisenbergFSD;
  double nasaFSD;
  double FSD;

  cmatrix chemName;
  dvector chemWeight;

  // check number of files if 1, stop
  checkInput(argc);

  fileName = argv[1];

  // get number of line
  num = lc(fileName);

  // allocate memory fof chemName and chemWeight
  cmatlloc(&chemName, num - 2, 16);
  dveclloc(&chemWeight, num - 2);

  // load file argv[1]
  loadFile(fileName, chemName, chemWeight);

  //start iterative calculation
  BSD = calcBSD(chemName, chemWeight);
  SSD = calcSSD(chemName, chemWeight);
  eisenbergFSD = calcEisenbergFSD(chemName, chemWeight);
  nasaFSD = calcNASAFSD(chemName, chemWeight);
  if(nasaFSD > eisenbergFSD) FSD = nasaFSD;
  else FSD = eisenbergFSD;

  printf("BSD = %lf\n", BSD);
  printf("SSD = %lf\n", SSD);
  printf("FSD = %lf\n", FSD);
  return 0;
}

void checkInput(int argc)
{
  if(argc != 2)
  {
    printf("No file input.\n");
    printf("Exit");
    exit(1);
  }
}

void loadFile(char *fileName, cmatrix chemName, dvector chemWeight)
{
  FILE *fp;
  char data[256];
  printf("Loading the file %s.\n", fileName);
  double num = lc(fileName);
  readFile(&fp, fileName, 1);
  for(int i = 0; i < num; i++)
  {
    fgets(data, 256, fp);
    if(strncmp(data, "{", 1) == 0) continue;
    else if(strncmp(data, "}", 1) == 0) break;
    else sscanf(data, "%s %lf", chemName.v2[i - 1], &chemWeight.v[i - 1]);
  }
  fclose(fp);
  printf("File successfully loaded.\n");
}

double sumWei(cmatrix *chemName, dvector *chemWeight)
{
  double num = (*chemWeight).dim;
  double wei = 0;
  for(int i = 0; i < num; i++)
  {
    double wp = (*chemWeight).v[i];
    char *name = (*chemName).v2[i];
    if(strcmp(name, "LOX") == 0)
    {
      wei += 7.8 / pow(wp, 1e0 / 3e0) * wp;
    }
    else if(strcmp(name, "LCH4") == 0)
    {
      wei += 7.8 / pow(wp, 1e0 / 3e0) * wp;
    }
  }
  return wei;
}

double sumWeo(cmatrix *chemName, dvector *chemWeight)
{
  double num = (*chemWeight).dim;
  double weo = 0;
  for(int i = 0; i < num; i++)
  {
    double wp = (*chemWeight).v[i];
    char *name = (*chemName).v2[i];
    if(strcmp(name, "LOX") == 0)
    {
      weo += 7.8 / pow(wp, 1e0 / 3e0) * wp;
    }
    else if(strcmp(name, "LCH4") == 0)
    {
      weo += 7.8 / pow(wp, 1e0 / 3e0) * wp;
    }
  }
  return weo;
}

double calcBSD(cmatrix chemName, dvector chemWeight)
{
  double z;
  double i;
  double dp;
  double r;
  double rOld;
  double rInput;

  char data[256];

  // calc w_{ei} and w_{eo}
  double sumOfWei;
  double sumOfWeo;
  sumOfWei = sumWei(&chemName, &chemWeight);
  sumOfWeo = sumWeo(&chemName, &chemWeight);

  // load initial R
  printf("Input initial R\n>>>");
  fgets(data, 256, stdin);
  sscanf(data, "%lf", &rInput);

  //start iterative calculation
  r = rInput;
  for(int j = 1; ;j++)
  {
    rOld = r;
    z = r / pow(sumOfWei, 1e0 / 3e0);
    i = pow(sumOfWei, 1e0 / 3e0) * 3.67 * pow(z, -1.08 + 0.0072 * log(z));
    if(i <= 140) dp = 1.379;
    else if(i < 400 && i > 140) dp = 1.379 * pow(140 / i, 2.4e-1);
    else if(i >= 400) dp = 1.073;
    r = 74 / pow(dp, 1 / 1.41) * pow(sumOfWeo, 1e0 / 3e0);
    if(j % 10 == 0) printf("%d %lf\n", j, r);
    if(fabs(r - rOld) < 1e-6) break;
  }
  return r;
}

double sumWp(dvector chemWeight)
{
  double wp;
  int num = chemWeight.dim;
  for(int i = 0; i < num; i++) wp += chemWeight.v[i];
  return wp;
}

double calcSSD(cmatrix chemName, dvector chemWeight)
{
  double ssd;
  double wp = sumWp(chemWeight);
  ssd = 59 * pow(wp, 0.21);
  return ssd;
}

double calcEisenbergFSD(cmatrix chemName, dvector chemWeight)
{
  double fsd;
  double wp = sumWp(chemWeight);
  double a = 1.0;
  fsd = 12.134 * sqrt(a) * pow(wp, 0.4058);
  return fsd;
}

double calcNASAFSD(cmatrix chemName, dvector chemWeight)
{
  double fsd;
  double wp = sumWp(chemWeight);
  double i = 12560;
  double a = 1;
  fsd = sqrt(8.58 * pow(10, 6e0) * a * pow(wp, 2e0 / 3e0) / i);
  return fsd;
}
