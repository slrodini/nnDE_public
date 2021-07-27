#include <adam.h>
#include <default.h>
#include <network.h>
#include <ran2.h>
#include <time.h>

long idum = -3;

double foo_IN(int n, bool mode)
{
  int nL = 7;

  int arch[7] = {n, 10, 10, 10, 10, 10, 1};
  multilayerD net1 = multilD_init_net(nL, arch);
  multilD_setMode(&net1, mode);

  double x[n];
  for (int i = 0; i < n; i++)
  {
    x[i] = i * 0.1 + 0.1;
  }
  clock_t begin = clock();
  multilD_FullEvaluate(&net1, x);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  multilD_free_net(&net1);
  return time_spent;
}

double foo_D1_IN(int n)
{
  int nL = 7;

  int arch[7] = {n, 10, 10, 10, 10, 10, 1};
  multilayerD1 net1 = multilD1_init_net(nL, arch);

  double x[n];
  for (int i = 0; i < n; i++)
  {
    x[i] = i * 0.1 + 0.1;
  }
  clock_t begin = clock();
  multilD1_FullEvaluate(&net1, x);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  multilD1_free_net(&net1);
  return time_spent;
}

double foo_OUT(int n, bool mode)
{
  int nL = 7;

  int arch[7] = {4, 10, 10, 10, 10, 10, n};
  multilayerD net1 = multilD_init_net(nL, arch);
  multilD_setMode(&net1, mode);

  double x[4];
  for (int i = 0; i < 4; i++)
  {
    x[i] = i * 0.1 + 0.1;
  }
  clock_t begin = clock();
  multilD_FullEvaluate(&net1, x);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  multilD_free_net(&net1);
  return time_spent;
}

double foo_D1_OUT(int n)
{
  int nL = 7;

  int arch[7] = {4, 10, 10, 10, 10, 10, n};
  multilayerD1 net1 = multilD1_init_net(nL, arch);

  double x[4];
  for (int i = 0; i < 4; i++)
  {
    x[i] = i * 0.1 + 0.1;
  }
  clock_t begin = clock();
  multilD1_FullEvaluate(&net1, x);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  multilD1_free_net(&net1);
  return time_spent;
}

int main()
{
  double iter = 200.0;
  int maxIO;
  printf("Please type the number of maximum input/output to be tested:\n");
  scanf("%d", &maxIO);
  if (maxIO <= 1)
    maxIO = 2;
  FILE *fp = fopen("../data/time_nI.dat", "w");
  fprintf(fp, "#nOutput\tFull_Hessian\t Diagonal_Hessian\tFirst_Derivatives\n");
  for (int nI = 1; nI <= maxIO; nI++)
  {
    double time_spent1 = 0.0;
    double time_spent2 = 0.0;
    double time_spent3 = 0.0;
    for (int j = 0; j < iter; j++)
    {
      time_spent1 += foo_IN(nI, false);
      time_spent2 += foo_IN(nI, true);
      time_spent3 += foo_D1_IN(nI);
    }
    printf("Inputs: %d \n", nI);
    fprintf(fp, "%d\t%e\t%e\t%e\n", nI, time_spent1 / iter, time_spent2 / iter, time_spent3 / iter);
  }
  fclose(fp);

  fp = fopen("../data/time_nO.dat", "w");
  fprintf(fp, "#nOutput\tFull_Hessian\t Diagonal_Hessian\tFirst_Derivatives\n");
  for (int nI = 1; nI <= maxIO; nI++)
  {
    double time_spent1 = 0.0;
    double time_spent2 = 0.0;
    double time_spent3 = 0.0;
    for (int j = 0; j < iter; j++)
    {
      time_spent1 += foo_OUT(nI, false);
      time_spent2 += foo_OUT(nI, true);
      time_spent3 += foo_D1_OUT(nI);
    }
    printf("Outputs: %d \n", nI);
    fprintf(fp, "%d\t%e\t%e\t%e\n", nI, time_spent1 / iter, time_spent2 / iter, time_spent3 / iter);
  }
  fclose(fp);
}
