# %%
from sympy import *
from sympy.abc import a,b,c,u
v = Symbol('varphi') # v instead of phi

from IPython.core.display import display

from einsteinpy.symbolic import MetricTensor, ChristoffelSymbols, RiemannCurvatureTensor,RicciTensor, RicciScalar,EinsteinTensor

# %%
def one_sheet_hyperboloid(u,v):
    x = a * cosh(u) * cos(v)
    y = b * cosh(u) * sin(v)
    z = c * sinh(u)
    return Matrix([x,y,z]).subs(b,a).subs(c,a)
def two_sheet_hyperboloid(u,v):
    x = a * sinh(u) * cos(v)
    y = b * sinh(u) * sin(v)
    z = c * cosh(u)
    return Matrix([x,y,z]).subs(b,a).subs(c,a)
# %%
def J(f,args=(u,v)):# jacobian matrix
    x1, x2 = args[0], args[1]
    X = Matrix([x1,x2])
    return f(x1,x2).jacobian(X)

def jacobian_matrix(fc,x,args=None):# jacobian matrix
    if args is None:
        args = x
    return fc(x).jacobian(Matrix(args))

# surface definition end
# %%
def eta(n):# Minkowski metric
    x = eye(n)
    x[-1,-1] = -1
    return x

def flat_metric(args,time_zero = True):
    r, s = args[0], args[1]
    # the number of plus one = r
    # the number of minus one = s
    n = r + s
    x = eye(n)
    
    for j in range(s):
        if time_zero:
            k = j
        else:
            k = -(j+1)
        
        x[k,k] = - 1
    return x
# flat_metric(1,2) = diag(-1,-1,1)
# %%

def g1(f):# metric
    j = J(f)
    return simplify((j.T @ j))

def metric_from_jacobian(fc,x,_sign,__time_zero = True, args_jac=None):
    jac = jacobian_matrix(fc,x,args_jac)
    # jac.T is the transpose of j
    return simplify((jac.T @ flat_metric(_sign,time_zero = __time_zero) @ jac))

_sign_of_ads2 = (1,2)
_sign_of_ds2 = (2,1)
# %%
def scalar(metric_matrix, x, output=True):
    m_obj = MetricTensor(metric_matrix, x)
    # print('metric'); display(m_obj);
    Ric = RicciTensor.from_metric(m_obj)
    R = RicciScalar.from_riccitensor(Ric)
    
    if output == True:
        print("Ricci tensor")
        display(Ric.tensor().simplify())
        
        print("Ricci scalar")
        display(R.simplify())
    
    return Ric, R

def einstein_tensor(metric, x, output=True):
    ################
    Ric, R = (scalar(metric, x, output=True))
    Ric, R = Ric.simplify(), R.simplify()
    ################
    G = Ric - tensorproduct(R / 2,metric)
    G = G.simplify()
    if output == True:
        print("Einstein tensor")
        display(G)
    
    return (G)
# %%
import sympy as sp

def get_spherical_vars(n,angles_symbol='\\theta', subscript=0):
    # Define symbolic variables
    # angle
    if subscript == 0:
        spherical_vars = sp.symbols(angles_symbol + '_{0:%d}' % (n-1))
    if subscript == 1:
        spherical_vars = sp.symbols(angles_symbol + '_{1:%d}' % (n))
    return spherical_vars

def spherical_to_cartesian(n,radius = sp.Symbol("r"), _angles_symbol = '\\theta'):
    # Define symbolic variables
    # angle
    spherical_vars = get_spherical_vars(n,angles_symbol=_angles_symbol)
        
    # Create the coordinate mappings
    cartesian_coords = []
    polar_angle = spherical_vars[0]# theta_0
    cartesian_coords.append(radius * sp.cos(polar_angle))
    product_sin = sp.sin(polar_angle)
    
    for j in range(1, n - 1):
        cartesian_coords.append(radius * sp.cos(spherical_vars[j]))
        cartesian_coords[j] *= product_sin
        product_sin *= sp.sin(spherical_vars[j])
    
    cartesian_coords.append(radius * product_sin)

    print("Cartesian coordinates:")
    print(cartesian_coords)
    
    return sp.Matrix(cartesian_coords)

# %%
def spherical_to_de_Sitter(n,radius = sp.Symbol("a")):
    # Define symbolic variables
    cartesian_coords_sphere = spherical_to_cartesian(n,radius=radius)
    u = sp.Symbol("u")
    cartesian_coords_sphere *= sp.cosh(u)
    cartesian_coords_sphere = list(cartesian_coords_sphere)
    # add the last coordinate
    cartesian_coords_sphere.append(radius * sp.sinh(u))
    return sp.Matrix(cartesian_coords_sphere)

# %%
def get_cartesian_vars(n):
    return symbols('x1:'+str(n+1))
def cartesian_to_metric_sphere(n):
    x = get_cartesian_vars(n)
    # x sub n + 1
    x_np1_squared = sum([x[j]**2 for j in range(n)])
    return ([[KroneckerDelta(j,k) + x[j] * x[k] / (1 - x_np1_squared) for j in range(n)] for k in range(n)])
def cartesian_to_metric_hyperbola(n):
    x = get_cartesian_vars(n)
    # x sub n + 1
    x_np1_squared = sum([x[j]**2 for j in range(n)])
    return ([[KroneckerDelta(j,k) - x[j] * x[k] / (x_np1_squared - 1) for j in range(n)] for k in range(n)]) 
# %%
from itertools import product

def display_geodesic(matrix_metric, x):
    m_obj = MetricTensor(matrix_metric, x)
    Gamma = ChristoffelSymbols.from_metric(m_obj).simplify()
    
    # variables
    # X = Function("X")(*x)
    s = Symbol("s")
    n = len(x)
    X = [Function("X^"+str(j))(s) for j in range(n)]
    
    L = range(n)
    for mu in L:
        print("The geodesic equation with respect to", end=" ")
        display(x[mu])
        term1 = X[mu].diff(s).diff(s)
        for j,k in product(L,L):
            term1 += Gamma[mu,j,k] * X[j].diff(s) * X[k].diff(s)
        display(term1)
    