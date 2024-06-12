# %%

# %%
# from Metrics import *
#reload
from importlib import reload # >=python 3.4
import Metrics
reload( Metrics )
# from Metrics import * 

from Metrics import SchwarzschildMetric,geodesics

from einsteinpy.symbolic import  ChristoffelSymbols
from sympy import symbols,Function,pi,latex
from sympy.core.symbol import Symbol

from IPython.core.display import display
from sympy.abc import lamda
# %%

m_obj = SchwarzschildMetric(eta00=-1,rs=True)
g = m_obj.tensor()
# display(g)
# %%
ch = ChristoffelSymbols.from_metric(m_obj)
Γ = ch.tensor().simplify()
display(Γ)
display(latex(Γ))
# %%
λ = lamda
# U = Function('u')(λ)
# V = Function('varphi')(λ)
T = Function('t')(λ)
R = Function('r')(λ)
Q = Function('theta')(λ)
P = Function('varphi')(λ)
x = symbols('t r theta phi')# x0 = t
r = x[1]; θ = x[2];
X = [T,R,Q,P]
M = [j for j in range(len(X))]
ch
# %%

Γ = ch.tensor().subs(r,R).subs(θ,Q)
for μ in M:
    print(μ)
    # _ans = geodesics(μ,X,Γ).subs(Q.diff(λ),0).subs(Q,pi/2).subs(Symbol('r_s'),1).simplify()
    # _ans = geodesics(μ,X,Γ).subs(Q.diff(λ),0).subs(Q,pi/2).simplify()
    _ans = geodesics(μ,X,Γ).simplify()
    display(_ans)
    # print(_ans)
    display(latex(_ans))
# %%

# %%
# import os
# print("Current working directory: {0}".format(os.getcwd()))
# # %%

# path = '/Users/komatsubararyosuke/LearningSummary/Theory_of_Relativity_jp_book/factory/Metrics/'
# try:
#     os.chdir(path)
#     print("Current working directory: {0}".format(os.getcwd()))
# except:
#     pass
