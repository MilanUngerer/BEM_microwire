#Preambulo
import numpy as np
import bempp.api
omega = 2.*np.pi*10.e9
e0 = 8.854*1e-12*1e-18
mu0 = 4.*np.pi*1e-7*1e6
mue = (1.)*mu0
ee = (16.)*e0
mui = (-2.9214+0.5895j)*mu0
ei = (82629.2677-200138.2211j)*e0
k = omega*np.sqrt(e0*mu0)
lam = 2*np.pi/k
nm = np.sqrt((ee*mue)/(e0*mu0))
nc = np.sqrt((ei*mui)/(e0*mu0))
alfa_m = mue/mu0
alfa_c = mui/mue
antena = np.array([[1e4],[0.],[0.]])
print "Numero de onda exterior:", k
print "Indice de refraccion matriz:", nm
print "Indice de refraccion conductor:", nc
print "Numero de onda interior matriz:", nm*k
print "Numero de onda interior conductor:", nm*nc*k
print "Indice de transmision matriz:", alfa_m
print "Indice de transmision conductor:", alfa_c
print "Longitud de onda:", lam, "micras"

#Importando mallas
matriz = bempp.api.import_grid('/home/milan/matriz_12x12x300_E16772.msh')
grid_0 = bempp.api.import_grid('/home/milan/BH12_a5_l10_E5550.msh')
grid_1 = bempp.api.import_grid('/home/milan/BH22_a5_l10_E5550.msh')
grid_2 = bempp.api.import_grid('/home/milan/BH32_a5_l10_E5550.msh')
grid_3 = bempp.api.import_grid('/home/milan/BH42_a5_l10_E5550.msh')

#Funciones de dirichlet y neumann
def dirichlet_fun(x, n, domain_index, result):
	result[0] = 1.*np.exp(1j*k*x[0])
def neumann_fun(x, n, domain_index, result):
	result[0] = 1.*1j*k*n[0]*np.exp(1j*k*x[0])

#Operadores multitrazo
Ai_m = bempp.api.operators.boundary.helmholtz.multitrace_operator(matriz, nm*k)
Ae_m = bempp.api.operators.boundary.helmholtz.multitrace_operator(matriz, k)
Ai_0 = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_0,nm*nc*k)
Ae_0 = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_0,nm*k)
Ai_1 = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_1,nm*nc*k)
Ae_1 = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_1,nm*k)
Ai_2 = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_2,nm*nc*k)
Ae_2 = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_2,nm*k)
Ai_3 = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_3,nm*nc*k)
Ae_3 = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_3,nm*k)

#Transmision en Multitrazo
Ae_m[0,1] = Ae_m[0,1]*(1./alfa_m)
Ae_m[1,1] = Ae_m[1,1]*(1./alfa_m)
Ai_0[0,1] = Ai_0[0,1]*alfa_c
Ai_0[1,1] = Ai_0[1,1]*alfa_c
Ai_1[0,1] = Ai_1[0,1]*alfa_c
Ai_1[1,1] = Ai_1[1,1]*alfa_c
Ai_2[0,1] = Ai_2[0,1]*alfa_c
Ai_2[1,1] = Ai_2[1,1]*alfa_c
Ai_3[0,1] = Ai_3[0,1]*alfa_c
Ai_3[1,1] = Ai_3[1,1]*alfa_c

#Acople interior y exterior
op_m = (Ai_m + Ae_m)
op_0 = (Ai_0 + Ae_0)
op_1 = (Ai_1 + Ae_1)
op_2 = (Ai_2 + Ae_2)
op_3 = (Ai_3 + Ae_3)

#Espacios
dirichlet_space_m = Ai_m[0,0].domain
neumann_space_m = Ai_m[0,1].domain
dirichlet_space_0 = Ai_0[0,0].domain
neumann_space_0 = Ai_0[0,1].domain
dirichlet_space_1 = Ai_1[0,0].domain
neumann_space_1 = Ai_1[0,1].domain
dirichlet_space_2 = Ai_2[0,0].domain
neumann_space_2 = Ai_2[0,1].domain
dirichlet_space_3 = Ai_3[0,0].domain
neumann_space_3 = Ai_3[0,1].domain

#Operadores identidad
ident_m = bempp.api.operators.boundary.sparse.identity(neumann_space_m, neumann_space_m, neumann_space_m)
ident_0 = bempp.api.operators.boundary.sparse.identity(neumann_space_0, neumann_space_0, neumann_space_0)
ident_1 = bempp.api.operators.boundary.sparse.identity(neumann_space_1, neumann_space_1, neumann_space_1)
ident_2 = bempp.api.operators.boundary.sparse.identity(neumann_space_2, neumann_space_2, neumann_space_2)
ident_3 = bempp.api.operators.boundary.sparse.identity(neumann_space_3, neumann_space_3, neumann_space_3)

#Operadores diagonales
op_m[1,1] = op_m[1,1] + 0.5 * ident_m * ((alfa_m -1)/alfa_m)
op_0[1,1] = op_0[1,1] + 0.5 * ident_0* (alfa_c - 1)
op_1[1,1] = op_1[1,1] + 0.5 * ident_1* (alfa_c - 1)
op_2[1,1] = op_2[1,1] + 0.5 * ident_2* (alfa_c - 1)
op_3[1,1] = op_3[1,1] + 0.5 * ident_3* (alfa_c - 1)

#Operadores entre mallas
SLP_m_0 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_m, dirichlet_space_0, dirichlet_space_0, nm*k)
SLP_0_m = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_0, dirichlet_space_m, dirichlet_space_m, nm*k)
DLP_m_0 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_m, dirichlet_space_0, dirichlet_space_0, nm*k)
DLP_0_m = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_0, dirichlet_space_m, dirichlet_space_m, nm*k)
ADLP_m_0 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_m, neumann_space_0, neumann_space_0, nm*k)
ADLP_0_m = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_0, neumann_space_m, neumann_space_m, nm*k)
HYP_m_0 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_m, neumann_space_0, neumann_space_0, nm*k)
HYP_0_m = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_0, neumann_space_m, neumann_space_m, nm*k)
SLP_0_1 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_0, dirichlet_space_1, dirichlet_space_1, nm*k)
DLP_0_1 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_0, dirichlet_space_1, dirichlet_space_1, nm*k)
ADLP_0_1 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_0, neumann_space_1, neumann_space_1, nm*k)
HYP_0_1 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_0, neumann_space_1, neumann_space_1, nm*k)
SLP_0_2 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_0, dirichlet_space_2, dirichlet_space_2, nm*k)
DLP_0_2 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_0, dirichlet_space_2, dirichlet_space_2, nm*k)
ADLP_0_2 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_0, neumann_space_2, neumann_space_2, nm*k)
HYP_0_2 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_0, neumann_space_2, neumann_space_2, nm*k)
SLP_0_3 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_0, dirichlet_space_3, dirichlet_space_3, nm*k)
DLP_0_3 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_0, dirichlet_space_3, dirichlet_space_3, nm*k)
ADLP_0_3 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_0, neumann_space_3, neumann_space_3, nm*k)
HYP_0_3 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_0, neumann_space_3, neumann_space_3, nm*k)
SLP_m_1 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_m, dirichlet_space_1, dirichlet_space_1, nm*k)
SLP_1_m = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_1, dirichlet_space_m, dirichlet_space_m, nm*k)
DLP_m_1 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_m, dirichlet_space_1, dirichlet_space_1, nm*k)
DLP_1_m = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_1, dirichlet_space_m, dirichlet_space_m, nm*k)
ADLP_m_1 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_m, neumann_space_1, neumann_space_1, nm*k)
ADLP_1_m = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_1, neumann_space_m, neumann_space_m, nm*k)
HYP_m_1 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_m, neumann_space_1, neumann_space_1, nm*k)
HYP_1_m = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_1, neumann_space_m, neumann_space_m, nm*k)
SLP_1_0 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_1, dirichlet_space_0, dirichlet_space_0, nm*k)
DLP_1_0 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_1, dirichlet_space_0, dirichlet_space_0, nm*k)
ADLP_1_0 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_1, neumann_space_0, neumann_space_0, nm*k)
HYP_1_0 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_1, neumann_space_0, neumann_space_0, nm*k)
SLP_1_2 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_1, dirichlet_space_2, dirichlet_space_2, nm*k)
DLP_1_2 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_1, dirichlet_space_2, dirichlet_space_2, nm*k)
ADLP_1_2 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_1, neumann_space_2, neumann_space_2, nm*k)
HYP_1_2 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_1, neumann_space_2, neumann_space_2, nm*k)
SLP_1_3 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_1, dirichlet_space_3, dirichlet_space_3, nm*k)
DLP_1_3 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_1, dirichlet_space_3, dirichlet_space_3, nm*k)
ADLP_1_3 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_1, neumann_space_3, neumann_space_3, nm*k)
HYP_1_3 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_1, neumann_space_3, neumann_space_3, nm*k)
SLP_m_2 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_m, dirichlet_space_2, dirichlet_space_2, nm*k)
SLP_2_m = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_2, dirichlet_space_m, dirichlet_space_m, nm*k)
DLP_m_2 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_m, dirichlet_space_2, dirichlet_space_2, nm*k)
DLP_2_m = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_2, dirichlet_space_m, dirichlet_space_m, nm*k)
ADLP_m_2 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_m, neumann_space_2, neumann_space_2, nm*k)
ADLP_2_m = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_2, neumann_space_m, neumann_space_m, nm*k)
HYP_m_2 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_m, neumann_space_2, neumann_space_2, nm*k)
HYP_2_m = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_2, neumann_space_m, neumann_space_m, nm*k)
SLP_2_0 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_2, dirichlet_space_0, dirichlet_space_0, nm*k)
DLP_2_0 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_2, dirichlet_space_0, dirichlet_space_0, nm*k)
ADLP_2_0 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_2, neumann_space_0, neumann_space_0, nm*k)
HYP_2_0 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_2, neumann_space_0, neumann_space_0, nm*k)
SLP_2_1 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_2, dirichlet_space_1, dirichlet_space_1, nm*k)
DLP_2_1 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_2, dirichlet_space_1, dirichlet_space_1, nm*k)
ADLP_2_1 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_2, neumann_space_1, neumann_space_1, nm*k)
HYP_2_1 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_2, neumann_space_1, neumann_space_1, nm*k)
SLP_2_3 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_2, dirichlet_space_3, dirichlet_space_3, nm*k)
DLP_2_3 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_2, dirichlet_space_3, dirichlet_space_3, nm*k)
ADLP_2_3 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_2, neumann_space_3, neumann_space_3, nm*k)
HYP_2_3 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_2, neumann_space_3, neumann_space_3, nm*k)
SLP_m_3 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_m, dirichlet_space_3, dirichlet_space_3, nm*k)
SLP_3_m = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_3, dirichlet_space_m, dirichlet_space_m, nm*k)
DLP_m_3 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_m, dirichlet_space_3, dirichlet_space_3, nm*k)
DLP_3_m = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_3, dirichlet_space_m, dirichlet_space_m, nm*k)
ADLP_m_3 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_m, neumann_space_3, neumann_space_3, nm*k)
ADLP_3_m = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_3, neumann_space_m, neumann_space_m, nm*k)
HYP_m_3 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_m, neumann_space_3, neumann_space_3, nm*k)
HYP_3_m = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_3, neumann_space_m, neumann_space_m, nm*k)
SLP_3_0 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_3, dirichlet_space_0, dirichlet_space_0, nm*k)
DLP_3_0 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_3, dirichlet_space_0, dirichlet_space_0, nm*k)
ADLP_3_0 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_3, neumann_space_0, neumann_space_0, nm*k)
HYP_3_0 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_3, neumann_space_0, neumann_space_0, nm*k)
SLP_3_1 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_3, dirichlet_space_1, dirichlet_space_1, nm*k)
DLP_3_1 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_3, dirichlet_space_1, dirichlet_space_1, nm*k)
ADLP_3_1 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_3, neumann_space_1, neumann_space_1, nm*k)
HYP_3_1 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_3, neumann_space_1, neumann_space_1, nm*k)
SLP_3_2 = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_3, dirichlet_space_2, dirichlet_space_2, nm*k)
DLP_3_2 = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_3, dirichlet_space_2, dirichlet_space_2, nm*k)
ADLP_3_2 = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_3, neumann_space_2, neumann_space_2, nm*k)
HYP_3_2 = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_3, neumann_space_2, neumann_space_2, nm*k)

#Matriz de operadores
blocked = bempp.api.BlockedOperator(10,10)

#Diagonal
blocked[0,0] = op_m[0,0]
blocked[0,1] = op_m[0,1]
blocked[1,0] = op_m[1,0]
blocked[1,1] = op_m[1,1]
blocked[2,2] = op_0[0,0]
blocked[2,3] = op_0[0,1]
blocked[3,2] = op_0[1,0]
blocked[3,3] = op_0[1,1]
blocked[4,4] = op_1[0,0]
blocked[4,5] = op_1[0,1]
blocked[5,4] = op_1[1,0]
blocked[5,5] = op_1[1,1]
blocked[6,6] = op_2[0,0]
blocked[6,7] = op_2[0,1]
blocked[7,6] = op_2[1,0]
blocked[7,7] = op_2[1,1]
blocked[8,8] = op_3[0,0]
blocked[8,9] = op_3[0,1]
blocked[9,8] = op_3[1,0]
blocked[9,9] = op_3[1,1]

#Contribucion hilos-matriz
blocked[0,2] = DLP_0_m
blocked[0,3] = -SLP_0_m
blocked[1,2] = -HYP_0_m
blocked[1,3] = -ADLP_0_m
blocked[0,4] = DLP_1_m
blocked[0,5] = -SLP_1_m
blocked[1,4] = -HYP_1_m
blocked[1,5] = -ADLP_1_m
blocked[0,6] = DLP_2_m
blocked[0,7] = -SLP_2_m
blocked[1,6] = -HYP_2_m
blocked[1,7] = -ADLP_2_m
blocked[0,8] = DLP_3_m
blocked[0,9] = -SLP_3_m
blocked[1,8] = -HYP_3_m
blocked[1,9] = -ADLP_3_m

#Contribucion hilos-hilos
blocked[2,4] = DLP_1_0
blocked[2,5] = -SLP_1_0
blocked[3,4] = -HYP_1_0
blocked[3,5] = -ADLP_1_0

#Contribucion hilos-hilos
blocked[2,6] = DLP_2_0
blocked[2,7] = -SLP_2_0
blocked[3,6] = -HYP_2_0
blocked[3,7] = -ADLP_2_0

#Contribucion hilos-hilos
blocked[2,8] = DLP_3_0
blocked[2,9] = -SLP_3_0
blocked[3,8] = -HYP_3_0
blocked[3,9] = -ADLP_3_0

#Contribucion matriz-hilos
blocked[2,0] = -DLP_m_0
blocked[2,1] = SLP_m_0
blocked[3,0] = HYP_m_0
blocked[3,1] = ADLP_m_0

#Contribucion hilos-hilos
blocked[4,2] = DLP_0_1
blocked[4,3] = -SLP_0_1
blocked[5,2] = -HYP_0_1
blocked[5,3] = -ADLP_0_1

#Contribucion hilos-hilos
blocked[4,6] = DLP_2_1
blocked[4,7] = -SLP_2_1
blocked[5,6] = -HYP_2_1
blocked[5,7] = -ADLP_2_1

#Contribucion hilos-hilos
blocked[4,8] = DLP_3_1
blocked[4,9] = -SLP_3_1
blocked[5,8] = -HYP_3_1
blocked[5,9] = -ADLP_3_1

#Contribucion matriz-hilos
blocked[4,0] = -DLP_m_1
blocked[4,1] = SLP_m_1
blocked[5,0] = HYP_m_1
blocked[5,1] = ADLP_m_1

#Contribucion hilos-hilos
blocked[6,2] = DLP_0_2
blocked[6,3] = -SLP_0_2
blocked[7,2] = -HYP_0_2
blocked[7,3] = -ADLP_0_2

#Contribucion hilos-hilos
blocked[6,4] = DLP_1_2
blocked[6,5] = -SLP_1_2
blocked[7,4] = -HYP_1_2
blocked[7,5] = -ADLP_1_2

#Contribucion hilos-hilos
blocked[6,8] = DLP_3_2
blocked[6,9] = -SLP_3_2
blocked[7,8] = -HYP_3_2
blocked[7,9] = -ADLP_3_2

#Contribucion matriz-hilos
blocked[6,0] = -DLP_m_2
blocked[6,1] = SLP_m_2
blocked[7,0] = HYP_m_2
blocked[7,1] = ADLP_m_2

#Contribucion hilos-hilos
blocked[8,2] = DLP_0_3
blocked[8,3] = -SLP_0_3
blocked[9,2] = -HYP_0_3
blocked[9,3] = -ADLP_0_3

#Contribucion hilos-hilos
blocked[8,4] = DLP_1_3
blocked[8,5] = -SLP_1_3
blocked[9,4] = -HYP_1_3
blocked[9,5] = -ADLP_1_3

#Contribucion hilos-hilos
blocked[8,6] = DLP_2_3
blocked[8,7] = -SLP_2_3
blocked[9,6] = -HYP_2_3
blocked[9,7] = -ADLP_2_3

#Contribucion matriz-hilos
blocked[8,0] = -DLP_m_3
blocked[8,1] = SLP_m_3
blocked[9,0] = HYP_m_3
blocked[9,1] = ADLP_m_3

#Condiciones de borde
dirichlet_grid_fun_m = bempp.api.GridFunction(dirichlet_space_m, fun=dirichlet_fun)
neumann_grid_fun_m = bempp.api.GridFunction(neumann_space_m, fun=neumann_fun)

#Discretizacion lado izquierdo
blocked_discretizado = blocked.strong_form()

#Discretizacion lado derecho
rhs = np.concatenate([dirichlet_grid_fun_m.coefficients, neumann_grid_fun_m.coefficients,np.zeros(dirichlet_space_0.global_dof_count), np.zeros(neumann_space_0.global_dof_count), np.zeros(dirichlet_space_1.global_dof_count), np.zeros(neumann_space_1.global_dof_count), np.zeros(dirichlet_space_2.global_dof_count), np.zeros(neumann_space_2.global_dof_count), np.zeros(dirichlet_space_3.global_dof_count), np.zeros(neumann_space_3.global_dof_count)])

#Sistema de ecuaciones
import inspect
from scipy.sparse.linalg import gmres
array_it = np.array([])
array_frame = np.array([])
it_count = 0
def iteration_counter(x):
	global array_it
	global array_frame
	global it_count
	it_count += 1
	frame = inspect.currentframe().f_back
	array_it = np.append(array_it, it_count)
	array_frame = np.append(array_frame, frame.f_locals["resid"])
	print it_count, frame.f_locals["resid"]
print("Shape of matrix: {0}".format(blocked_discretizado.shape))
x,info = gmres(blocked_discretizado, rhs, tol=1e-5, callback = iteration_counter, maxiter = 50000)
print("El sistema fue resuelto en {0} iteraciones".format(it_count))
np.savetxt("Solucion.out", x, delimiter=",")

#Campo interior
interior_field_dirichlet_m = bempp.api.GridFunction(dirichlet_space_m, coefficients=x[:dirichlet_space_m.global_dof_count])
interior_field_neumann_m = bempp.api.GridFunction(neumann_space_m,coefficients=x[dirichlet_space_m.global_dof_count:dirichlet_space_m.global_dof_count + neumann_space_m.global_dof_count])

#Campo exterior
exterior_field_dirichlet_m = interior_field_dirichlet_m
exterior_field_neumann_m = interior_field_neumann_m*(1./alfa_m)

#Calculo campo en antena
slp_pot_ext_m = bempp.api.operators.potential.helmholtz.single_layer(dirichlet_space_m, antena, k)
dlp_pot_ext_m = bempp.api.operators.potential.helmholtz.double_layer(dirichlet_space_m, antena, k)
Campo_en_antena = (dlp_pot_ext_m * exterior_field_dirichlet_m - slp_pot_ext_m * exterior_field_neumann_m).ravel() + np.exp(1j*k*antena[0])
print "Valor del campo en receptor:", Campo_en_antena

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from matplotlib import rcParams
rcParams["font.family"] = "serif"
rcParams["font.size"] = 20
pyplot.figure(figsize = (15,10))
pyplot.title("Convergence")
pyplot.plot(array_it, array_frame, lw=2)
pyplot.xlabel("iteration")
pyplot.ylabel("residual")
pyplot.grid()
pyplot.savefig("Convergence.pdf")
