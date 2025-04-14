import numpy as np
from pyscf import gto, scf, dft, tddft
#from pyscf.geomopt.berny_solver import optimize
from pyscf.geomopt.geometric_solver import optimize

import pyscf
#from pyscf.tools.mo_mapping import mo_comps

#from benchmarking_utils import setup_logger, get_cpu_timings

#log = setup_logger()

# Define molecule geometry
mol = gto.Mole()
mol.atom = '''
O 0 0 0
H 0.958 0 0
H -0.958 0 0
'''
mol.basis = '6-31g'  # Adjust basis set as needed
mol.build()

#cpu0 = get_cpu_timings()

# Ground state calculation with DFT
# Choose the desired functional (e.g., B3LYP)
mf = dft.RKS(mol).density_fit(auxbasis='def2-universal-jkfit')
mf.xc = 'B3LYP'
mf.grids.becke_scheme=dft.gen_grid.original_becke
mf.run()
#cpu0 = log.timer('H2O %s RKS'%basis, *cpu0)

#mf = dft.RKS(mol).x2c().set(xc='pbe0').run()
## Switch to xcfun because 3rd order GGA functional derivative is not
## available in libxc
#mf._numint.libxc = dft.xcfun
postmf = tddft.TDA(mf).run(nstates=4)
# PySCF-1.6.1 and newer supports the .Gradients method to create a grad
# object after grad module was imported. It is equivalent to call the
# .nuc_grad_method method.
#from pyscf import grad
#g = postmf.Gradients()
#g.kernel(state=1)

#postmf.nuc_grad_method().kernel()
postmf.nuc_grad_method().as_scanner(state=1).optimizer(solver='geomeTRIC').kernel()

## Optimize S0 geometry
#geom = optimize(mf)
#geom.run()

#test = geom.RKS().set(xc = 'B3LYP', conv_tol = 1e-6) \
	#.density_fit(auxbasis='def2-universal-jkfit') \
	#.apply(pyscf.scf.addons.remove_linear_dep_) \
	#.run()

## Print results
#print("S0 energy:", mf.e_tot)
#print ('test S0',test.e_tot)
#print('Delta', abs(mf.e_tot - test.e_tot))

## Excited state calculations using TDDFT
#mytd = tddft.TDA(mf)
#mytd.nstates = 2  # Calculate S1 and T1

### Optimize S1
##mytd.density = mf.make_rdm1()  # Use S0 density for TDDFT
##geom.restart()
##geom.init_guess = mytd.kernel()[1][0]  # Use S1 excitation for guess
##geom.max_step = 0.01  # Adjust step size for excited states
##geom.tol_grad = 5e-4  # Adjust tolerance for excited states
##geom.run()

### Optimize T1
##geom.restart()
##geom.init_guess = mytd.kernel()[2][0]  # Use T1 excitation for guess
##geom.run()

## Print results
#print("S0 energy:", mf.e_tot)

