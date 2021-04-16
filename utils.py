import meshio
import numpy as np
from collections import deque
from dolfin import (Mesh,MeshEditor,edges,vertex_to_dof_map,
                   Function,FunctionSpace,FiniteElement,
                   assemble,dx,Constant,solve,Matrix,Vector,
                   TrialFunction,TestFunction,inner,grad,action,
                   PointSource)

def import_from_gmsh(fname):
    "Convert from gmsh to dolfin"
    
    # read with meshio
    msh = meshio.read(fname)
    
    # create a DOLFIN mesh (assuming 2d)
    gdim,tdim = 2,2
    mm = Mesh()
    editor = MeshEditor()
    editor.open(mm,"triangle",gdim,tdim)

    npt = msh.points.shape[0]
    nc  = msh.get_cells_type("triangle").shape[0]

    editor.init_vertices_global(npt,npt)
    editor.init_cells_global(nc,nc)

    for i,p in enumerate(msh.points):
        editor.add_vertex(i,p[:2])

    for i,c in enumerate(msh.get_cells_type("triangle")):
        editor.add_cell(i,c)

    editor.close()
    
    # domains
    md = mm.domains()
    md.init(tdim)
    markers = {}
    
    if 'gmsh:physical' not in msh.cell_data_dict:
        # no markers at all
        return mm, markers

    phy = msh.cell_data_dict['gmsh:physical']
    if 'triangle' in phy:
        for eid,val in enumerate(phy['triangle']):
            md.set_marker((eid,val),2)

    if 'line' in phy:
        mm.init(0,1)
        p2e = mm.topology()(0,1)

        for l,k in zip(msh.get_cells_type("line"),phy['line']):
            e = set(p2e(l[0])).intersection(p2e(l[1])).pop()
            md.set_marker((e,k),1)

    if 'vertex' in phy:
        for eid,val in zip(msh.get_cells_type("vertex"),phy['vertex']):
            md.set_marker((eid[0],val),0)

    # names
    markers = tuple({n:v.item() for n,(v,d) in msh.field_data.items() if d == dim} for dim in range(tdim+1))
    
    return mm,markers


def curvilinear_coordinate_1d(mb,p0=0,function=None):
    "Returns parametrization of a curve"
    
    # edge-to-vertex connectivity
    EE = np.zeros((mb.num_cells(),mb.num_vertices()),dtype=bool)
    for e in edges(mb):
        EE[e.index(),e.entities(0)] = True

    # vertex-to-vertex connectivity (via edges)
    PP = EE.T @ EE
    np.fill_diagonal(PP,False)
    mmap = -np.ones(PP.shape[0],dtype=int)

    # order vertices
    mmap[0] = p0
    for k in range(PP.shape[0]-1):
        neig = np.where(PP[mmap[k],:])[0]
        mmap[k+1] = neig[1] if neig[0] in mmap else neig[0]

    # cumulative length of edges
    l = np.linalg.norm(np.diff(mb.coordinates()[mmap,:],axis=0),axis=1)
    s = np.r_[0,np.cumsum(l)]

    if function is None:
        P1e = FiniteElement("CG",mb.ufl_cell(),1)
        Ve  = FunctionSpace(mb,P1e)
        function = Function(Ve)
        function.vector()[vertex_to_dof_map(Ve)[mmap]] = s
        return function
    else:
        Ve = function.function_space()
        function.vector()[vertex_to_dof_map(Ve)[mmap]] = s

        
def eikonal_1d(mb,p0=0,function=None):
    "Compute distance from p0 on set of edges"
    
    # edge-to-vertex connectivity
    EE = np.zeros((mb.num_cells(),mb.num_vertices()),dtype=bool)
    for e in edges(mb):
        EE[e.index(),e.entities(0)] = True

    # vertex-to-vertex connectivity (via edges)
    PP = EE.T @ EE
    np.fill_diagonal(PP,False)

    # initial solution is inf everywhere
    sol = np.empty(PP.shape[0])
    sol.fill(np.inf)

    # initial conditions
    active  = deque([p0])
    sol[p0] = 0.0
    
    # fast marching on edges
    x = mb.coordinates()
    while active:
        curr = active.pop()
        neig = np.where(PP[curr,:])[0]
        ll = sol[curr] + np.linalg.norm(x[neig,:] - x[curr,:],axis=1)
        up = neig[ll < sol[neig]]
        active.extend(up)
        sol[neig] = np.minimum(sol[neig],ll)

    # return solution
    if function is None:
        P1e = FiniteElement("CG",mb.ufl_cell(),1)
        Ve  = FunctionSpace(mb,P1e)
        function = Function(Ve)
        function.vector()[vertex_to_dof_map(Ve)] = sol
        return function
    else:
        Ve = function.function_space()
        function.vector()[vertex_to_dof_map(Ve)] = sol

        
def uniform_distribution(mb,p0,delta,function=None):
    "Uniform distribution centered at p0 and width 2*delta"

    # compute the distance
    if function is None:
        function = eikonal_1d(mb,p0=p0)
    else:
        eikonal_1d(mb,p0=p0,function=function)

    # distribution
    dval = function.vector().get_local()
    dist = np.zeros_like(dval)
    dist[dval <= delta] = 1.0
    function.vector().set_local(dist)
    
    # normalize
    area = assemble(function*dx)
    function.vector()[:] /= area

    
def gaussian_distribution(mb,mu,sigma,function=None,lumping=True,nsteps=100):
    "Gaussian distribution via heat equation"
    
    tend = 0.5*sigma**2
    dt = Constant(tend/nsteps,name="smooth")

    # prepare the problem
    P1e = FiniteElement("CG",mb.ufl_cell(),1)
    Ve  = FunctionSpace(mb,P1e)

    u,v = TrialFunction(Ve),TestFunction(Ve)
    uold = Function(Ve)
    
    if lumping:
        # diffusion
        K = assemble(dt*inner(grad(u),grad(v))*dx)
        # we use mass lumping to avoid negative values
        Md = assemble(action(u*v*dx,Constant(1.0)))
        # full matrix (divide my mass)
        M = Matrix(K)
        M.zero()
        M.set_diagonal(Md)
        A = M + K
    else:
        a = u*v*dx + dt*inner(grad(u),grad(v))*dx
        L = uold*v*dx
        A = assemble(a)

    # initial conditions
    dist = function or Function(Ve)

    dist.vector().zero()
    PointSource(Ve,mu,1.0).apply(dist.vector())

    # iterations
    for t in range(nsteps):
        uold.assign(dist)
        if lumping:
            solve(A,dist.vector(),M*uold.vector())
        else:
            b = assemble(L)
            solve(A,dist.vector(),b)

    # normalize
    area = assemble(dist*dx)
    dist.vector()[:] /= area
    
    if function is None:
        return dist
