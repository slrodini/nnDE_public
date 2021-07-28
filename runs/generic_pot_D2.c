#include <default.h>
#include <minimize.h>
#include <network.h>
#include <ran2.h>
#include <simAnnealing.h>
#include <time.h>
long idum = -2;

//struct to be fed to the adam minimizer
typedef struct
{
  int nX;
  double *x;
  multilayerD *net;
} minim_par;

// specified the potential to be used
#include <harmonic_potential.h>

//lower and upper limit of the xgrid
double xmin = -4.0;
double xmax = +4.0;

// wrappers for the multilD functions
double get(multilayerD *net) { return multilD_get(0, net); }
double get_grad(int j, multilayerD *net) { return multilD_get_grad(0, j, net); }
double get_d(multilayerD *net) { return multilD_get_d(0, 0, net); }
double get_grad_d(int j, multilayerD *net)
{
  return multilD_get_grad_d(0, 0, j, net);
}
double get_d2(multilayerD *net) { return multilD_get_d2(0, 0, 0, net); }
double get_grad_d2(int j, multilayerD *net)
{
  return multilD_get_grad_d2(0, 0, 0, j, net);
}

// function to compute both the loss functions and its gradient w.r.t. the network parameters
void Full_chi2(void *addPar, double *c2, double *grad)
{
  minim_par *mp = (minim_par *)addPar;
  multilayerD *net = mp->net;
  double res = 0.0;

  // auxiliary vectors for the gradient computation
  double psi2_grad[net->nPar];
  double hpsi_grad[net->nPar];

  for (int j = 0; j < net->nPar; j++)
  {
    grad[j] = 0.0;
    psi2_grad[j] = 0.0;
    hpsi_grad[j] = 0.0;
  }
  //initialize the energy and the norm
  double norm = 0.0;
  double en = 0.0;
  double dx = mp->x[1] - mp->x[0];

  double div = (double)mp->nX + 1.0;
  for (int i = 0; i < mp->nX; i++)
  {
    double x = mp->x[i];
    // evaluation of the output of the network and the gradient ath the point x
    multilD_FullEvaluate(net, &x);

    // get the wavefunction as the output of the network
    double psi = get(net);

    // compute the action of the Hamiltonian on the w.f.
    double h_psi = -0.5 * get_d2(net) + potential(x) * psi;
    // check
    if (isnan(h_psi) != 0)
    {
      printf("err %lf \t %lf \t %lf\n", get_d2(net), potential(x), x);
      exit(-1);
    }
    // simple rectangular integral approxiamtion
    en += psi * h_psi * dx;
    norm += psi * psi * dx;

    // integral of the gradients w.r.t. the parameters
    for (int j = 0; j < net->nPar; j++)
    {
      psi2_grad[j] += 2.0 * psi * get_grad(j, net) * dx;
      hpsi_grad[j] +=
          psi * dx *
          (-0.5 * get_grad_d2(j, net) + potential(x) * get_grad(j, net));
      hpsi_grad[j] += h_psi * get_grad(j, net) * dx;
    }
  }
  // computing the functional < psi | H | psi > / < psi | psi >
  en /= norm;
  // compute the difference of the norm from 1
  double diff = norm - 1.0;

  //save the loss function value
  *c2 = en + pow(diff, 2);

  // computing the gradient w.r.t. the parameters of the loss function (squared sum of the energy functional and the norm difference)
  for (int j = 0; j < net->nPar; j++)
  {
    grad[j] = (hpsi_grad[j] / norm - en * psi2_grad[j] / norm);
    grad[j] += 2.0 * diff * psi2_grad[j];
  }
}

double test_solution(void *addPar, double en)
{
  // Test the proposed network with the final energy using the eigenvalue problem (H-E)| psi > = 0
  minim_par *mp = (minim_par *)addPar;
  multilayerD *net = mp->net;
  double res = 0.0;
  //Finer grid of integration
  double n = 1e+5;
  double dx = (xmax - xmin) / (n - 1);
  for (int i = 0; i < n; i++)
  {
    double x = xmin + i * dx;
    multilD_FullEvaluate(net, &x);

    double psi = multilD_get(0, net);
    double h_psi = -0.5 * multilD_get_d2(0, 0, 0, net) + potential(x) * psi;

    res += pow(h_psi - en * psi, 2);
  }
  return res / n;
}

double get_energy(multilayerD *net)
{
  // Compute the energy via the energy functional < psi | H | psi > / < psi | psi >
  double res = 0.0;
  double norm = 0.0;
  double en = 0.0;
  double n = 1e+5;
  double dx = (xmax - xmin) / (n - 1);

  for (int i = 0; i < n; i++)
  {
    double x = xmin + i * dx;
    multilD_FullEvaluate(net, &x);

    double psi = multilD_get(0, net);
    double h_psi = -0.5 * multilD_get_d2(0, 0, 0, net) + potential(x) * psi;

    en += psi * h_psi * dx;
    norm += psi * psi * dx;
  }
  en /= norm;
  printf("Norm: %lf\n", norm);
  return en;
}

double get_norm(multilayerD *net)
{
  // Compute the norm (squared) of the w.f.
  double res = 0.0;
  double norm = 0.0;
  double n = 1e+5;
  double dx = (xmax - xmin) / (n - 1);

  for (int i = 0; i < n; i++)
  {
    double x = xmin + i * dx;
    multilD_FullEvaluate(net, &x);

    double psi = multilD_get(0, net);

    norm += psi * psi * dx;
  }
  return norm;
}

int main()
{
  clock_t begin = clock();
  int nI = 1;
  int nO = 1;
  // Specify the architecture
  int arch[8] = {nI, 10, 10, 10, 10, 10, 10, nO};
  int nL = sizeof(arch) / sizeof(arch[0]);
  // Get the number of parameters for the minimizer
  int nPar = multilD_getNpar(nL, arch);
  //Init the network
  multilayerD net = multilD_init_net(nL, arch);
  // Set the computation mode to diagonal Hessian
  multilD_setMode(&net, true);
  minim_par mp;
  mp.net = &net;

  // Init the integration grid
  mp.nX = 10000;
  mp.x = (double *)malloc(sizeof(double) * mp.nX);
  for (int i = 0; i < mp.nX; i++)
  {
    mp.x[i] = xmin + (xmax - xmin) * (double)(i) / ((double)mp.nX - 1.0);
  }

  // Perform the minimization with multiple pass of ADAM
  double chiMin = minimize(net.par, nPar, (void *)(&mp), Full_chi2);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time elapsed in s: %lf\n", time_spent);

  double energy = get_energy(&net);

  printf("Final enrgy: %lf\t Expected: %lf\n", energy, ground_state_energy());

  multilD_free_net(&net);

  return 0;
}
