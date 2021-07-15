#include <default.h>
#include <network.h>
#include <ran2.h>
#include <time.h>
long idum = -2;

int main() {
  int nI = 2;
  int nO = 1;
  int arch[10] = {nI, 10, 10, 10, 10, 10, 10, 10, 10, nO};
  int nL = 10;
  int nPar = multilD_getNpar(nL, arch);
  multilayerD net = multilD_init_net(nL, arch);
  multilD_load_net(&net, "networkPar.dat");

  FILE *fp = fopen("data.dat", "w");
  double input[2];
  for (int i = 0; i < 1000; i++) {
    input[1] = 0.05;
    double x = -1.0 + 2.0 * (double)i / 999.0;
    input[0] = x;

    multilD_Evaluate(&net, input);
    double re = multilD_get(0, &net);
    fprintf(fp, "%lf\t%lf\t%lf\n", input[0], re, x * x * cos(M_PI * x));
  }

  // multilD_save_net(&net, "networkPar.dat");

  multilD_free_net(&net);

  return 0;
}
