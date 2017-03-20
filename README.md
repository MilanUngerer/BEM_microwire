# BEM_microwire
Boundary Element Method based code to solve Helmholz equation with BEM++ library. Scattering problem for a composite with microwire inclusions.

In the last time composite materials have been subject of multiple studies due to the great utility of the special properties that they present. In this work the study was oriented to a material composed of a dielectric matrix with inclusions of conductive materials called ferromagnetic micro-wires. These wires present the peculiarity of changing their electromagnetic properties in front of an external stimulus. When applying a magnetic field or a mechanical tension the permittivity and permeability change, thus changing the way in which the incident electromagnetic radiation is scattered.

This work is oriented to find the optimal configuration of the micro-wires to maximize the change of signal in order to design a stress sensor for dielectrics. The radiation will be of the microwave type and a mechanic tension stimulus will be used. In order to reach the proposed objectives a program was developed that allows to calculate the electromagnetic field scattered by the composite at any point in which it is desired to position the antenna.

For the calculation of the field the Helmholtz equation was used through its boundary integral formulation which allows a simpler manipulation for the calculate of the field and to realize the interaction between the different obstacles. Because this equation does not have an analytical solution it is necessary to approximate it by means of some numerical method, in this work we used the Boundary Element Method with, the library for Python, BEM ++.

The optimal configuration is studied under 2 experiments, the first one corresponds to varying the distance of the wires with screen configuration. Another alternative is to study how the received signal changes according the change of the volumetric composition of th wires in the composite.

In the case of a screen configuration it can be seen that there is a maximum signal change when the wires are closer together in the center, so it would be advisable to position the wires as close as possible. In addition, it can be seen that for the volumetric percentage have a maximum with 4 wires reaching a change of around -90%.

The main limitations of the program is that due to the particular properties of the inserted conductors it's difficult to achieve the convergence of the program. 
