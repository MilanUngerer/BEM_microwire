# BEM_microwire

## Abstract

In the last time composite materials have been subject of multiple studies due to the great utility of the special properties that they present. In this work the study was oriented to a material composed of a dielectric matrix with inclusions of conductive materials called ferromagnetic micro-wires. These wires present the peculiarity of changing their electromagnetic properties in front of an external stimulus. When applying a magnetic field or a mechanical tension the permittivity and permeability change, thus changing the way in which the incident electromagnetic radiation is scattered.

This work is oriented to find the optimal configuration of the micro-wires to maximize the change of signal in order to design a stress sensor for dielectrics. The radiation will be of the microwave type and a mechanic tension stimulus will be used. In order to reach the proposed objectives a program was developed that allows to calculate the electromagnetic field scattered by the composite at any point in which it is desired to position the antenna.

For the calculation of the field the Helmholtz equation was used through its boundary integral formulation which allows a simpler manipulation for the calculate of the field and to realize the interaction between the different obstacles. Because this equation does not have an analytical solution it is necessary to approximate it by means of some numerical method, in this work we used the Boundary Element Method with, the library for Python, BEM ++.

The optimal configuration is studied under 2 experiments, the first one corresponds to varying the distance of the wires with screen configuration. Another alternative is to study how the received signal changes according the change of the volumetric composition of th wires in the composite.

In the case of a screen configuration it can be seen that there is a maximum signal change when the wires are closer together in the center, so it would be advisable to position the wires as close as possible. In addition, it can be seen that for the volumetric percentage have a maximum with 4 wires reaching a change of around -90%.

The main limitations of the program is that due to the particular properties of the inserted conductors it's difficult to achieve the convergence of the program. 

## How this works?

This project consist in 3 main codes for the calculation of the electromagnetic field. The first one consist in the file that contain the info for the problem and the meshes that we will use, his name is 'info.txt'. The second one called 'conector.py' is a python file that when it's executed create a file named 'ejecutor.py' that will solve the problem, so the last step is execute the file 'ejecutor.py' and wait for a response. The work use mainly numpy and BEM++ libraries.

Besides exists two notebooks called 'Permittivity' and 'Permeability' that calculate the properties of the microwires.

For more information about the formulation of the codes visit https://github.com/MilanUngerer/BEM_microwire/blob/master/documentation/Tesis.pdf and https://github.com/MilanUngerer/BEM_microwire/blob/master/Code_instructions.ipynb
