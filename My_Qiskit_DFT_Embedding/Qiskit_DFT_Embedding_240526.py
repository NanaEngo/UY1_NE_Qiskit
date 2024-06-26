from qiskit_nature.second_q.drivers import MethodType, PySCFDriver
from qiskit_nature.second_q.operators import ElectronicIntegrals, PolynomialTensor
from qiskit_nature.second_q.operators.symmetric_two_body import S8Integrals
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit_nature.second_q.problems import ElectronicBasis, ElectronicStructureProblem
from qiskit_nature.second_q.mappers import JordanWignerMapper, ParityMapper, TaperedQubitMapper
from qiskit_algorithms import NumPyEigensolver
from qiskit_nature.second_q.algorithms import ExcitedStatesEigensolver

import numpy as np

# setup driver
omega = 1.0
driver = PySCFDriver(
    atom="O 0.0 0.0 0.115; H 0.0 0.754 -0.459; H 0.0 -0.754 -0.459",
    #atom=open('PSPCz.xyz').read(),
    basis="6-31g*",
    method=MethodType.RKS,
    xc_functional=f"ldaerf + lr_hf({omega})",
    xcf_library="xcfun",
)

'''
1. Run the reference calculation
'''
driver.run_pyscf()
mf = driver._calc

'''
 2. Build the total Hamiltonian outside of a Problem instance
'''

from qiskit_nature.second_q.formats.qcschema import QCSchema
from qiskit_nature.second_q.formats import qcschema_translator

qcschema = driver.to_qcschema(include_dipole=True)
hamiltonian = qcschema_translator._get_mo_hamiltonian_direct(qcschema)

e_nuc = hamiltonian.nuclear_repulsion_energy
hamiltonian.electronic_integrals.alpha += PolynomialTensor({"": e_nuc})
hamiltonian.nuclear_repulsion_energy = None

'''
3. Prepare the active space by initializing it with the total problem size information
'''

# Somes variables needed
total_num_particles = driver._mol.nelec
total_num_spatial_orbitals = driver._mol.nao # hamiltonian.register_length
total_num_electrons = driver._mol.nelectron

# Alpha-spin
num_alpha = total_num_particles[0]
orbital_occupations = np.asarray([1.0] * num_alpha
                                 + [0.0] * (total_num_spatial_orbitals - num_alpha))
# Beta-spin
num_beta = total_num_particles[1]
orbital_occupations_b = np.asarray([1.0] * num_beta
                                   + [0.0] * (total_num_spatial_orbitals - num_beta))
# Initialize the transformer
transformer =  ActiveSpaceTransformer(2,2)

# Prepare the active space
transformer.prepare_active_space(total_num_particles, total_num_spatial_orbitals,
                                 occupation_alpha=orbital_occupations,
                                 occupation_beta=orbital_occupations_b,)
# Determine the active space
as_orbitals = transformer._determine_active_space(total_num_electrons, total_num_spatial_orbitals)[0]

# Use the active space transformer to reduce the total Hamiltonian to the active space
reduced_hamiltonian = transformer.transform_hamiltonian(hamiltonian)

'''
4. Setup the Problem instance
'''

reduced_hamiltonian.constants["inactive energy"] = mf.e_tot - driver._mol.energy_nuc()
problem = ElectronicStructureProblem(reduced_hamiltonian)
problem.basis = ElectronicBasis.MO
problem.num_spatial_orbitals = transformer._num_spatial_orbitals
problem.num_particles = transformer._num_electrons
problem.orbital_energies = mf.mo_energy[as_orbitals]
problem.reference_energy = mf.e_tot

'''
5. Mapping the Hamiltonian to spin space
'''

mapper = ParityMapper(num_particles=problem.num_particles)
mapper = problem.get_tapered_mapper(mapper)

'''
6. Setup solver
'''

algo = NumPyEigensolver(k=100)
algo.filter_criterion = problem.get_default_filter_criterion()

solver = ExcitedStatesEigensolver(mapper, algo)
result = solver.solve(problem)

print(f"Total ground state energy = {result.total_energies[0]:.4f}")
print(f"Total first excited state energy = {result.total_energies[1]:.3f}")
print(f"Total second excited state energy = {result.total_energies[2]:.3f}")
