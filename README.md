# nnDE_public

<h3> nnDE: a library for first and second order neural network derivatives</h3>

This software provide the possibility of computing the first and second derivatives with respect to the input, as well as the gradient with respect of all the parameters of an arbitrary deep feed-forward neural network.

In the header file `network.h' can be found the supported network structures as well as all the relevant functions, which are briefly summarized here.

**Structures** <br>
<code>multilayerD</code> network structure for second order derivatives computations<br>
<code>multilayerD1</code> network structure for first order derivatives computations<br>
<code>multilayer</code> network structure for just parameter gradients computations<br>
<code>multiL_act</code> prototype for activation functions<br>

**Functions**<br>
We outline the functions for the *multilayerD* struct, being the functions for the other structures similar (*mutatis mutandis*)<br>
The architecture of the network is provided as an array of integers, where each entries set the number of nodes in the specific layer.<br><br>
<code>multilD_getNpar</code> assigned an architecture return the total number of network parameters, inputs: (number of  layers, architecture array)<br>
<code>multilD_init_net</code> initialize the network, returns a copy of the network, inputs: (number of layers, architecture array)<br>
<code>multilD_setMode</code> set the mode to full Hessian calculation or diagonal Hessian calculation, inputs: (pointer to the network, mode -> true = diagonal Hessian)<br>
<code>multilD_set_act</code> set all the activatin functions of the hidden layers to the supplied one, inputs: (pointer to the network, ptr to activation function, ptr to derivative of a.f., ptr to double der. of a.f., ptr to triple der of a.f.)<br>
<code>multilD_set_act_one</code> set one  a.f. to the supplied one, inputs: (ptr to the network, int indicating which layer has to be modified, ptr to activation function, ptr to derivative of a.f., ptr to double der. of a.f., ptr to triple der of a.f.)<br>
<code>multilD_(save/load)_net</code> save or load, respectively, the network, inputs: (ptr to the network, file name)<br>
<code>multilD_free_net</code> free the memory associated to the network, input: (ptr to the network)<br>
<code>multilD_Evaluate</code> evaluate the output of the network, the input of the netwrok must be supplied separately, inputs: (ptr to the network)<br>
<code>multilD_EvaluateParGradient</code> evaluate gradient w.r.t. the parameters of the network, the input of the netwrok must be supplied separately, inputs: (ptr to the network)<br>
<code>multilD_FullEvaluate</code> evaluate the output of the network and the gradient w.r.t. the parameters of the network, inputs: (ptr to the network, input of the network); this function should be preferred over the previous two<br>
<code>multiD_get_grad</code> get the partial derivative of the output w.r.t. a parameter, inputs: (index of output, index of the chosen parameter, ptr to the netwrok)<br>
<code>multiD_get_grad_d</code> get partial derivative w.r.t. the input of the partial derivative of the output w.r.t. a parameter, inputs: (index of output, index of the input, index of the chosen parameter, ptr to the netwrok)<br>
<code>multiD_get_grad_d2</code> get double partial derivative w.r.t. the input of the partial derivative of the output w.r.t. a parameter, inputs: (index of output, index of the first input, indx of the second input, index of the chosen parameter, ptr to the netwrok)<br>
<code>multiD_get</code> get the output of the network, inputs: (index of output, ptr to the netwrok)<br>
<code>multiD_get_d</code> get the partial derivative of the output of the network w.r.t. the input, inputs: (index of output, index of the input, ptr to the netwrok)<br>
<code>multiD_get_d2</code> get the double partial derivative of the output of the network w.r.t. the input, inputs: (index of output, first index of the input, second index of the input, ptr to the netwrok)<br>

**Usage**
In the parent directory a make file is provided.<br> 
Run <code>make</code> to compile both the main files generic_pot_D1.c and generic_pot_D2.c that are in the runs directory.<br>
Run <code>make d*der</code> where *=1,2 to compile either generic_pot_D1.c or generic_pot_D2.c respectivley.<br>