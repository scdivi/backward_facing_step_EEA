#! /usr/bin/env python3

from nutils import mesh, function, solver
from nutils.expression_v2 import Namespace
import numpy, treelog, util

# geometry parameters
L = 4       # Length of total channel
H = 1       # Height of wide channel
l = 1       # Length of narrow channel
h = 0.5     # Height of narrow channel

# FEM approximation parameters
nelems = 8   # Number of elements
degree = 2   # Polynomial degree

# construct geometry
xverts = numpy.linspace(0,L,nelems+1)
yverts = numpy.linspace(0,H,nelems+1)
domain, geom = mesh.rectilinear([xverts, yverts])

# trim the boundary
x, y = geom
exact = abs(x/l+y/h)+abs(x/l-y/h)-2
domain = domain.trim(exact, maxrefine=0)

# construct namespaces
ns = Namespace()
ns.x = geom
ns.define_for('x', gradient='∇', normal='n', jacobians=('dV', 'dS'))
ns.basis = domain.basis('h-std', degree=degree).vector(domain.ndims)
ns.u = function.dotarg('lhs', ns.basis)

# parameters
ns.μ     = 1e-3 
ns.γ     = 500
ns.ε_ij  = '(∇_j(u_i) + ∇_i(u_j)) / 2'
ns.v_nij = '(∇_j(basis_ni) + ∇_i(basis_nj)) / 2'
ns.δ     = function.eye(domain.ndims)
ns.Σ     = function.ones([domain.ndims])
ns.σ_ij  = 'μ (2 ε_ij + γ ∇_k(u_k) δ_ij)' 

# inlet velocity
ns.umax  = .06
ns.h     = .5
ns.H     = H-.5
ns.uinx  = 'umax (x_1 - h) ((H - (x_1 - h)) / ((H / 2)^2))'

# constraints 
sqr  = domain.boundary['trimmed,bottom,top'].integral('u_i u_i dS' @ ns, degree=degree*2)
sqr += domain.boundary['left'].integral('(u_0 - uinx) (u_0 - uinx) dS' @ ns, degree=degree*2)
cons = solver.optimize('lhs', sqr, droptol=1e-15)

# residuals
res  = domain.integral('2 μ ε_ij v_nij dV' @ ns, degree=degree*2)
res += domain.integral('γ μ ∇_i(u_i) ∇_j(basis_nj) dV' @ ns, degree=degree*2)
lhs  = solver.solve_linear('lhs', res, constrain=cons)

# number of degrees of freedom
ndofs = len(ns.basis)
treelog.user(f'Number of DOFs: {ndofs}')

# quantity of interest
ns.xc     = function.Array.cast([1.25, 0.25])
ns.r      = 0.2
ns.pi     = numpy.pi
ns.weight = '( 1 / (2 pi r^2) )  exp( - ((x_i - xc_i) (x_i - xc_i)) / ( 2 r^2 ) )'  
ns.Q      = 'u_0 weight'
qoi       = domain.integrate('Q dV' @ ns, degree=(degree+1)*2, arguments=dict(lhs=lhs))
treelog.user(f'Quantity of interest: {qoi}')

# post-process
util.plot_solution(domain, ns, lhs=lhs)
