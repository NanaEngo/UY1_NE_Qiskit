import numpy as np
from qiskit_nature.second_q.drivers import MethodType, PySCFDriver
%
from pyscf import ao2mo
from qiskit_nature.second_q.operators.symmetric_two_body import fold
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.operators import ElectronicIntegrals, PolynomialTensor
%
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

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

mo_coeff, mo_coeff_b = driver._expand_mo_object(mf.mo_coeff, array_dimension=3)
h1_a = mf.get_hcore()
h1_a = np.dot(np.dot(mo_coeff.T, h1_a), mo_coeff)
if mo_coeff_b is not None:
    h1_b = np.dot(np.dot(mo_coeff_b.T, h1_a), mo_coeff_b)
else:
    h1_b = None

# h2_aa = driver._mol.intor("int2e", aosym=8)
h2_aa = fold(ao2mo.full(driver._mol, mo_coeff, aosym=4))
if mo_coeff_b is not None:
    h2_bb = fold(ao2mo.full(driver._mol, mo_coeff_b, aosym=4))
    h2_ba = fold(ao2mo.general(
            driver._mol,
            [mo_coeff_b, mo_coeff_b, mo_coeff, mo_coeff],
            aosym=4,
        )
    )
else:
    h2_bb = None
    h2_ba = None

hamiltonian = ElectronicEnergy.from_raw_integrals(
    h1_a,
    h2_aa,
    h1_b,
    h2_bb,
    h2_ba
)
from qiskit_nature.second_q.formats.qcschema import QCSchema
from qiskit_nature.second_q.formats import qcschema_translator

qcschema = driver.to_qcschema(include_dipole=True)
hamiltonian = qcschema_translator._get_mo_hamiltonian_direct(qcschema)
hamiltonian.nuclear_repulsion_energy = driver._mol.energy_nuc()

# To included the constant nuclear_repulsion_energy in the second_q_op
e_nuc = hamiltonian.nuclear_repulsion_energy
hamiltonian.electronic_integrals.alpha += PolynomialTensor({"": e_nuc})
hamiltonian.nuclear_repulsion_energy = None

'''
3. Prepare the active space by initializing it with the total problem size information
'''
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer, BasisTransformer
import numpy as np

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
transformer._determine_active_space(total_num_electrons, total_num_spatial_orbitals)[0]

# Use the active space transformer to reduce the total Hamiltonian to the active space
reduced_hamiltonian = transformer.transform_hamiltonian(hamiltonian)

'''
4. Setup the Problem instance
'''
from qiskit_nature.second_q.problems import ElectronicBasis, ElectronicStructureProblem

reduced_hamiltonian.constants["inactive energy"] = mf.e_tot - driver._mol.energy_nuc()

problem = ElectronicStructureProblem(reduced_hamiltonian)
problem.basis = ElectronicBasis.MO
problem.num_spatial_orbitals = transformer._num_spatial_orbitals
problem.num_particles = (transformer._num_electrons // 2, transformer._num_electrons // 2)
problem.orbital_energies = mf.mo_energy[as_orbitals]
problem.reference_energy = mf.e_tot

'''
Problem Electronic properties
'''
from qiskit_nature.second_q.properties import (
    AngularMomentum,
    ElectronicDensity,
    Magnetization,
    ParticleNumber,
    ElectronicDipoleMoment
)

norb = len(as_orbitals) # number of active orbitals
problem.properties.particle_number = ParticleNumber(norb)
problem.properties.magnetization = Magnetization(norb)
problem.properties.angular_momentum = AngularMomentum(norb)
problem.properties.electronic_density = ElectronicDensity.from_orbital_occupation(
    problem.orbital_occupations,
    problem.orbital_occupations_b
)

'''
5. mapping the Hamiltonian to spin space
'''
from qiskit_nature.second_q.mappers import JordanWignerMapper, ParityMapper, TaperedQubitMapper

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
