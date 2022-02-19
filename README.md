<h1 align="center" style="display: block; font-size: 2.5em; font-weight: bold; margin-block-start: 1em; margin-block-end: 1em;">
<br><br><strong> Adiabatic master equation and quantum trajectories</strong>
</h1>

![Latest release](https://img.shields.io/github/v/release/aregtech/areg-sdk?label=%20%F0%9F%93%A3%20Latest%20release&style=flat&logoColor=b0c0c0&labelColor=363D44)

<!-- markdownlint-disable -->
<td><img src="https://img.shields.io/badge/MATLAB-R2020a-BLUE.svg" alt="MATLAB solution"/></td>

---
Code of adiabatic master equation and quantum trajectories used in paper
These are the codes of adiabatic master equation (ame.m) and quantum trajectories (aqt.m) used in paper:
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.97.022116

The Hamiltonian is a 8-qubit ferromagnetic Ising spin chain in a transverse field, with annealing schedule attached as 'DW1_parameters.txt'. This relatively complex problem illustrates a lot of subtle points in the implementation of ame and aqt, and can be easily reduced to smaller and simpler problems. 


![alt text](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/8-qubit_chain.png)






The master equation is:

![alt text](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/ame1.png)
![alt text](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/ame2.png)









Both ame.m and aqt.m use sparse matrix for computation,

![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/sparsem%20(1).png)

and rotation for every ode step.

![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/rotation%20(1).png)

aqt.m uses parallel computation, with sorting of jump operators,

![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/sortingjump%20(1).png)

and backtracking for error control.

![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/backtracking%20(1).png)



## Tutorial: UCL 4 qubit gadgets (and a single qubit)
The folder "tutorial_4_qubit_gadget" contains a tutorial on how to run the codes. 
![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/4-qubit.png)

It has a "onequbit_demo" for how to run codes on a single qubit. The pdf file "kawacode.pdf" contains a presentation on excecutions.
