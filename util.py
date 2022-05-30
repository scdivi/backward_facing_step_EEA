from matplotlib import collections, colors
import matplotlib.pyplot as plt
from nutils import function, solver, sparse

def plot_solution(domain, ns, **arguments):
    
    # streamlines
    ns.streambasis = domain.basis('h-std', degree=2)[1:]  # remove first dof to obtain non-singular system
    ns.stream = function.dotarg('streamdofs', ns.streambasis)  # stream function
    ns.ε = function.levicivita(2)
    sqr = domain.integral('Σ_i (u_i - ε_ij ∇_j(stream))^2 dV' @ ns, degree=4)
    arguments['streamdofs'] = solver.optimize('streamdofs', sqr, arguments=arguments)  # compute streamlines

    # post-processing
    bezier  = domain.sample('bezier', 9)
    x, umag, stream = bezier.eval(['x_i', 'sqrt(u_i u_i)', 'stream'] @ ns, **arguments)
    
    # export plots
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect = 'equal')
    ax.autoscale(enable=True, axis='both', tight=True)
    im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, umag, shading='gouraud', cmap='jet')
    fig.colorbar(im, orientation='horizontal')
    ax.tricontour(x[:,0], x[:,1], bezier.tri, stream, 16, colors='k', linestyles='dotted', linewidths=.5, zorder=9)
    ax.add_collection(collections.LineCollection(x[bezier.hull], colors='w', linewidths=.5, alpha=.2))
    ax.axis('off')
    plt.show()
    
def plot_mesh(domain, geom):
    
    # post-processing
    bezier   = domain.sample('bezier', 2)
    x, nulls = bezier.eval([geom, 0.])
    
    # export plots
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect = 'equal')
    ax.autoscale(enable=True, axis='both', tight=True)
    cmap = colors.ListedColormap("darkgray")
    ax.tripcolor(x[:,0], x[:,1], bezier.tri, nulls, shading='gouraud', cmap=cmap)
    ax.add_collection(collections.LineCollection(x[bezier.hull], colors='k', linewidths=.5, alpha=.2))
    ax.axis('off')
    plt.show()

def plot_dualsolution(domain, ns, **arguments):

    # compute values 
    bezier  = domain.sample('bezier', 9)
    x, dual = bezier.eval(['x_i', 'sqrt(z_i z_i)'] @ ns,**arguments) 
    
    # export plots
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect = 'equal')
    ax.autoscale(enable=True, axis='both', tight=True)
    im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, dual, shading='gouraud', cmap='jet')
    fig.colorbar(im, orientation='horizontal')
    ax.add_collection(collections.LineCollection(x[bezier.hull], colors='w', linewidths=.5, alpha=.2))
    plt.show()

def plot_dualerror(domain, ns, **arguments):

    # compute values 
    bezier  = domain.sample('bezier', 9)
    #x, dual = bezier.eval(['x_i', 'sqrt((z_i - psi_i) (z_i - psi_i))'] @ ns,**arguments) 
    x, dual = bezier.eval(['x_i', 'indicatorfunc'] @ ns,**arguments) 
    
    # export plots
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect = 'equal')
    ax.autoscale(enable=True, axis='both', tight=True)
    im = ax.tripcolor(x[:,0], x[:,1], bezier.tri, dual, shading='gouraud', cmap='jet')
    fig.colorbar(im, orientation='horizontal')
    ax.add_collection(collections.LineCollection(x[bezier.hull], colors='w', linewidths=.5, alpha=.2))
    plt.show()

def plot_convergence(ndofs, qoi):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # axis
    ax.set_xlabel(r'#dofs')
    ax.set_ylabel(r'Quantity of interest')
    # log plots
    ax.loglog(ndofs, qoi, 'bo-')
    # grids
    ax.grid(b=True, which='minor', axis="x",linestyle=':')
    ax.grid(b=True, which='minor', axis="y",linestyle=':')
    # tight layout
    fig.set_tight_layout(True)
    plt.show()

def get_support_vector_basis(domain, basis, mask):
    '''determine local support in the form of a dof -> transform chains mapping.'''
    f = function.kronecker(basis.sum(range(1, basis.ndim)), pos=domain.f_index, length=len(domain), axis=0)
    supp = [[] for i in range(len(basis))]
    for ielem, idof in zip(*sparse.indices(domain.sample('gauss', 1).integrate_sparse(f))):
        supp[idof].append(domain.transforms[ielem])

    refine = set()
    for idof in mask.nonzero()[0]:
        transforms = supp[idof]
        minlen = min(len(trans) for trans in transforms)
        refine.update(trans for trans in transforms if len(trans) == minlen)

    return refine
