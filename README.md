# nnDE_public

This software provide the possibility of computing the first and second derivatives with respect to the input, as well as the gradient with respect of all the parameters of an arbitrary deep feed-forward neural network.

In the header file `network.h' can be found the supported network structures as well as all the relevant functions, which are briefly summarized here.

**Structures** \n
*multilayerD* network structure for second order derivatives computations\n
*multilayerD1* network structure for first order derivatives computations\n
*multilayer* network structure for just parameter gradients computations\n
*multiL_act* prototype for activation functions\n

**Functions**\n
We outline the functions for the *multilayerD* struct, being the functions for the other structures similar (*mutatis mutandis*)\n
The architecture of the network is provided as an array of integers, where each entries set the number of nodes in the specific layer.\n
*multilD_getNpar* assigned an architecture return the total number of network parameters, inputs: (number of  layers, architecture array)\n
*multilD_init_net* initialize the network, returns a copy of the network, inputs: (number of layers, architecture array)\n
*multilD_setMode* set the mode to full Hessian calculation or diagonal Hessian calculation, inputs: (pointer to the network, mode -> true = diagonal Hessian)\n
*multilD_set_act* set all the activatin functions of the hidden layers to the supplied one, inputs: (pointer to the network, ptr to activation function, ptr to derivative of a.f., ptr to double der. of a.f., ptr to triple der of a.f.)\n
*multilD_set_act_one* set one  a.f. to the supplied one, inputs: (ptr to the network, int indicating which layer has to be modified, ptr to activation function, ptr to derivative of a.f., ptr to double der. of a.f., ptr to triple der of a.f.)\n
*multilD_(save/load)_net* save or load, respectively, the network, inputs: (ptr to the network, file name)\n
*multilD_free_net* free the memory associated to the network, input: (ptr to the network)\n
*multilD_Evaluate* evaluate the output of the network, the input of the netwrok must be supplied separately, inputs: (ptr to the network)\n
*multilD_EvaluateParGradient* evaluate gradient w.r.t. the parameters of the network, the input of the netwrok must be supplied separately, inputs: (ptr to the network)\n
*multilD_FullEvaluate* evaluate the output of the network and the gradient w.r.t. the parameters of the network, inputs: (ptr to the network, input of the network); this function should be preferred over the previous two\n