# %%
from sympy import *
from sympy.physics.quantum import TensorProduct
# %%
n_global = 4
def pauli_matrices(mu: int):
    if   mu == 0:
        sigma = eye(2)# sigma0
    elif mu == 1:
        sigma = Matrix(([0,1],[1,0]))# sigma1
    elif mu == 2:
        sigma = Matrix(([0,-1j],[1j,0]))# sigma2
    elif mu == 3:
        sigma = Matrix(([1,0],[0,-1]))# sigma3
    
    return sigma

# %%
def gamma_matrices_dirac(mu):
    # Return gamma matrices in the Dirac representation
    n = n_global
    if   mu == 0:
        gamma = TensorProduct(Matrix([[1,0],[0,-1]]), pauli_matrices(0))# gamma^0
    elif (mu == 1) or (mu == 2) or (mu == 3):
        gamma = TensorProduct(Matrix([[0,1],[-1,0]]) , pauli_matrices(mu))# gamma^mu
    elif mu == 5:
        gamma = 1j * prod([gamma_matrices_dirac(mu) for mu in range(n)])# gamma^5
    
    return gamma
    
def gamma_matrices_weyl(mu, sign_0 = 1, sign_1 = 1):
    # Return gamma matrices in som representation
    n = n_global
    if   mu == 0:
        gamma = sign_0 * TensorProduct(Matrix([[0,1],[1,0]]), pauli_matrices(0))# gamma^0
    elif (mu == 1) or (mu == 2) or (mu == 3):
        gamma = sign_1 * TensorProduct(Matrix([[0,1],[-1,0]]) , pauli_matrices(mu))# gamma^mu
    elif mu == 5:
        gamma = 1j * prod([gamma_matrices_weyl(mu) for mu in range(n)])# gamma^5
    
    return gamma

def gamma_matrices_majorana(mu):
    n = n_global
    if mu == 0:
        gamma = TensorProduct(Matrix([[0,1],[1,0]]), pauli_matrices(2))
    elif mu == 1:
        gamma = TensorProduct(Matrix([[1,0],[0,1]]), pauli_matrices(3) * 1j)
    elif mu == 2:
        gamma = TensorProduct(Matrix([[0,-1],[1,0]]), pauli_matrices(2))
    elif mu == 3:
        gamma = TensorProduct(Matrix([[1,0],[0,1]]), pauli_matrices(1) * (-1j))
    elif mu == 5:
        gamma = 1j * prod([gamma_matrices_majorana(mu) for mu in range(n)])# gamma^5
    
    return gamma

def gamma_matrices_pauli3(j):
    n = 3
    if j == 0:
        return pauli_matrices(3)
    if j == 1:
        return 1j * pauli_matrices(1)
    if j == 2:
        return 1j * pauli_matrices(2)
    if j == 5:
        return 1j * prod([gamma_matrices_pauli3(mu) for mu in range(n)])# gamma^5

def gamma_matrices(mu, rep = "Dirac"):
    if rep == "Dirac":
        return gamma_matrices_dirac(mu)
    
    elif rep == "Weyl":
        return gamma_matrices_weyl(mu)
    elif rep == "Majorana":
        return gamma_matrices_majorana(mu)
    elif rep == "Pauli3":
        return gamma_matrices_pauli3(mu)

# %%
def display_dirac_equation(rep="Dirac"):
    n = n_global
    x = symbols('x0:'+str(n))
    psi = Matrix([Function("\psi_"+str(j))(*x) for j in range(n)])

    _ans = zeros(n,1)

    for mu in range(n):    
        _ans += 1j * gamma_matrices(mu, rep="Dirac") @ psi.diff(x[mu])

    display(_ans)

# %%
def spin_matrices(j, rep = "Dirac"):
    if j == 1:
        spin_factor = gamma_matrices(2, rep=rep) @ gamma_matrices(3, rep=rep)
    if j == 2:
        spin_factor = gamma_matrices(3, rep=rep) @ gamma_matrices(1, rep=rep)
    if j == 3:
        spin_factor = gamma_matrices(1, rep=rep) @ gamma_matrices(2, rep=rep)

    return 1j * spin_factor / 2
    
# %%
def pauli_and_vector_sum(k):
    s = [pauli_matrices(j) * k[j] for j in range(1,n_global)]
    _ans = zeros(2,2)
    for _ in s:
        _ans += _
    return _ans
    
def eq4_43():
    a = symbols('a0:'+str(3+1), commutative=False)
    b = symbols('b0:'+str(3+1), commutative=False)

    sa = pauli_and_vector_sum(a)
    sb = pauli_and_vector_sum(b)
    
    return expand(sa * sb)
    
# %%
def commutator(A,B):
    return A @ B - B @ A
def anti_commutator(A,B):
    return A @ B + B @ A

def sigma_matrices(mu, nu, rep = "Dirac"):
    gamma_mu = gamma_matrices(mu, rep = rep)
    gamma_nu = gamma_matrices(nu, rep = rep)
    return 1j * commutator(gamma_mu, gamma_nu) / 2


# %%
