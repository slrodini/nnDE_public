#include <adam.h>
#include <default.h>
#include <network.h>
#include <ran2.h>
#include <time.h>

long idum = -3;

double foo(int n, bool mode)
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

double foo_D1(int n)
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

int main()
{
  double iter = 200.0;
  FILE *fp = fopen("../data/time_nI.dat", "w");
  fprintf(fp, "#nOutput\tFull_Hessian\t Diagonal_Hessian\tFirst_Derivatives\n");
  for (int nI = 1; nI <= 60; nI++)
  {
    double time_spent1 = 0.0;
    double time_spent2 = 0.0;
    double time_spent3 = 0.0;
    for (int j = 0; j < iter; j++)
    {
      time_spent1 += foo(nI, false);
      time_spent2 += foo(nI, true);
      time_spent3 += foo_D1(nI);
    }
    // time_spent2 = foo_thr(nI);
    printf("Inputs: %d  Time: %e \n", nI, time_spent1 / iter);
    fprintf(fp, "%d\t%e\t%e\t%e\n", nI, time_spent1 / iter, time_spent2 / iter, time_spent3 / iter);
  }
  fclose(fp);
}
