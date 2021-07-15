#include <default.h>

typedef double (*multiL_act)(double);

typedef struct {
  int nI, nH, nO, nL, nPar;
  // architecture, arch[0]=nInput arch[back] = nOutput
  int *arch;
  int *offsetW, *offsetB;
  // hidden layers (just the linear comb + the activated neuron)
  double **li;
  double **s_li, **sd_li, **sd2_li, **sd3_li;
  double **ds_li_dx, **dsd_li_dx;

  double **d2s_li_dx2, **d2sd_li_dx2, **d2sd2_li_dx2;

  // double *res;
  double *res_d;
  double **res_d2;
  double *par;
  double *parGrad;
  double *parGrad_d;
  double *parGrad_d2;
  multiL_act *act_fun;
  multiL_act *d_act_fun;
  multiL_act *d2_act_fun;
  multiL_act *d3_act_fun;
  int maxNH;

  double *input;

} multilayerD;

typedef struct {
  int nI, nH, nO, nL, nPar;
  // architecture, arch[0]=nInput arch[back] = nOutput
  int *arch;
  int *offsetW, *offsetB;
  // hidden layers (just the linear comb + the activated neuron)
  double **li;
  double **s_li, **sd_li;

  // double *res;
  double *par;
  double *parGrad;
  multiL_act *act_fun;
  multiL_act *d_act_fun;
  int maxNH;

  double *input;

} multilayer;

// multilayer functions

int multil_getNpar(int nL, int *arch);
multilayer multil_init_net(int nL, int *arch);
void multil_set_act(multilayer *net, multiL_act *activ, multiL_act *d_activ);
void multil_save_net(multilayer *net, const char fileName[]);
void multil_load_net(multilayer *net, const char fileName[]);
void multil_free_net(multilayer *net);

void multil_Evaluate(multilayer *net, double *x);
void multil_EvaluateParGradient(multilayer *net, double *x);
void multil_FullEvaluate(multilayer *net, double *x);

double multil_get_grad(int out, int par, multilayer *net);
double multil_get(int out, multilayer *net);

// multilayerD functions
int multilD_getNpar(int nH, int *arch);
multilayerD multilD_init_net(int nL, int *arch);
void multilD_set_act(multilayerD *net, multiL_act *activ, multiL_act *d_activ,
                     multiL_act *d2_activ, multiL_act *d3_activ);
void multilD_set_act_one(multilayerD *net, int i, multiL_act activ,
                         multiL_act d_activ, multiL_act d2_activ,
                         multiL_act d3_activ);
void multilD_save_net(multilayerD *net, const char fileName[]);
void multilD_load_net(multilayerD *net, const char fileName[]);
void multilD_free_net(multilayerD *net);

void multilD_Evaluate(multilayerD *net);
void multilD_EvaluateParGradient(multilayerD *net);
void multilD_FullEvaluate(multilayerD *net, double *x);

double multilD_get_grad(int out, int par, multilayerD *net);
double multilD_get_grad_d(int out, int i, int par, multilayerD *net);
double multilD_get_grad_d2(int out, int i1, int i2, int par, multilayerD *net);
double multilD_get(int out, multilayerD *net);
double multilD_get_d(int out, int i, multilayerD *net);
double multilD_get_d2(int out, int i1, int i2, multilayerD *net);
