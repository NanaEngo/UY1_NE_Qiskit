import numpy as np
from qiskit_nature.second_q.drivers import MethodType, PySCFDriver
#
from pyscf import ao2mo
from qiskit_nature.second_q.operators.symmetric_two_body import fold
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.operators import ElectronicIntegrals, PolynomialTensor
#
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from pyscf.tools import cubegen
#
from qiskit_nature.second_q.problems import ElectronicBasis, ElectronicStructureProblem
from typing import cast
#
from qiskit_nature.second_q.properties import (
    AngularMomentum,
    ElectronicDensity,
    Magnetization,
    ParticleNumber,
    ElectronicDipoleMoment
)
from qiskit_nature.second_q.operators import Tensor
#
from qiskit_nature.second_q.mappers import ParityMapper, TaperedQubitMapper
#
from qiskit_algorithms import NumPyEigensolver
from qiskit_nature.second_q.algorithms import ExcitedStatesEigensolver

# setup driver
omega = 1.0
driver = PySCFDriver(
    #atom="O 0.0 0.0 0.115; H 0.0 0.754 -0.459; H 0.0 -0.754 -0.459",
    atom=open('PSPCz.xyz').read(),
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
mf.dump_scf_summary()

'''
 2. Build the total Hamiltonian outside of a Problem instance
 Use of QCSchema
'''
from qiskit_nature.second_q.formats.qcschema import QCSchema
from qiskit_nature.second_q.formats import qcschema_translator 
from qiskit_nature.second_q.operators import PolynomialTensor

qcschema = driver.to_qcschema(include_dipole=True)
hamiltonian = qcschema_translator._get_mo_hamiltonian_direct(qcschema)

hamiltonian.nuclear_repulsion_energy = driver._mol.energy_nuc()

# '''
#  2. Build the total Hamiltonian outside of a Problem instance

# one-body integrals should be in the form
# h[p,q]= \int \phi_p(x)* (T + V_{ext}) \phi_q(x) dx

# two-body integrals should be in the form
# h[p,q,r,s] = \int \phi_p(x) * \phi_q(y) * V_{elec-elec} \phi_r(y) \phi_s(x) dxdy

# Using molecular orbitals \phi_j(x) = \sum_{ij} A_i(x) mo_coeff_{i,j} where A_i(x) are the atomic orbitals.
# '''

# mo_coeff, mo_coeff_b = driver._expand_mo_object(mf.mo_coeff, array_dimension=3)
# hcore = mf.get_hcore() # equivalent to int1e_kin + int1e_nuc (in ao)
# h1_a = np.dot(np.dot(mo_coeff.T, hcore), mo_coeff)
# #h1_a = mo_coeff.T @ hcore @ mo_coeff in mo
# if mo_coeff_b is not None:
#     h1_b = np.dot(np.dot(mo_coeff_b.T, hcore), mo_coeff_b)
#     #h1_b =  mo_coeff_b.T @ hcore @ mo_coeff_b
# else:
#     h1_b = None

# h2_aa = driver._mol.intor("int2e", aosym=8) # in ao
# h2_aa = fold(ao2mo.full(driver._mol, mo_coeff, aosym=4)) # in mo
# if mo_coeff_b is not None:
#     h2_bb = fold(ao2mo.full(driver._mol, mo_coeff_b, aosym=4))
#     h2_ba = fold(ao2mo.general(
#             driver._mol,
#             [mo_coeff_b, mo_coeff_b, mo_coeff, mo_coeff],
#             aosym=4,
#         )
#     )
# else:
#     h2_bb = None
#     h2_ba = None

# #h2 = ao2mo.kernel(driver._mol.intor("int2e"), mo_coeff)
# #h2 = ao2mo.restore(1, h2, len(mo_coeff))

# #h2 = h2.T(0, 2, 3, 1)
# #if isinstance(h2, tuple):
#     #h2_aa, h2_ab, h2_bb = h2
# #else:
#     #h2_aa, h2_ab, h2_bb = h2, None, None
# #h2_ba = h2_ab.T

# hamiltonian = ElectronicEnergy.from_raw_integrals(
#     h1_a,
#     h2_aa,
#     h1_b,
#     h2_bb,
#     h2_ba
# )
# hamiltonian.nuclear_repulsion_energy = driver._mol.energy_nuc()

# # # To included the constant nuclear_repulsion_energy in the second_q_op
# # e_nuc = hamiltonian.nuclear_repulsion_energy
# # hamiltonian.electronic_integrals.alpha += PolynomialTensor({"": e_nuc})

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

# Output cube files for active orbitals that can read by Jmol. Avogadro, Gabedit
for i in as_orbitals:
    cubegen.orbital(
        driver._mol, 
        f'Test_DFT_{i+1}.cube', 
        mf.mo_coeff[:, i])
    
'''
4. Setup the Problem instance
'''
reduced_hamiltonian.constants["inactive energy"] = mf.e_tot - driver._mol.energy_nuc()
print(as_orbitals)
print(reduced_hamiltonian.second_q_op())

problem = ElectronicStructureProblem(reduced_hamiltonian)
problem.basis = ElectronicBasis.MO
problem.reference_energy = mf.e_tot
problem.num_spatial_orbitals = transformer._num_spatial_orbitals

problem.orbital_occupations = np.diag(
    cast(np.ndarray, transformer._active_density.alpha["+-"])
)[transformer._active_alpha_indices]
problem.orbital_occupations_b = np.diag(
    cast(np.ndarray, transformer._active_density.beta["+-"])
)[transformer._active_beta_indices]
problem.num_particles = (
    round(sum(problem.orbital_occupations)),
    round(sum(problem.orbital_occupations_b)),
)

mo_energy, mo_energy_b = driver._expand_mo_object(driver._calc.mo_energy)
if mo_energy is not None:
    problem.orbital_energies = mo_energy[transformer._active_alpha_indices]
if mo_energy_b is not None:
    problem.orbital_energies_b = mo_energy_b[transformer._active_beta_indices]

''' 
5. Setup Problem Electronic properties
'''
for prop in problem.properties:
    if isinstance(prop, (Magnetization, ParticleNumber)):
        problem.properties.add(prop.__class__(problem.num_spatial_orbitals))
    
    elif isinstance(prop, ElectronicDipoleMoment):
        problem.properties.electronic_dipole_moment = (
            _transform_electronic_dipole_moment(
                prop,
                transformer._density_total,
                transformer._active_density,
                transformer._active_basis,
                transformer.__class__.__name__,
            )
        )

    elif isinstance(prop, ElectronicDensity):
        transformed = transformer._active_basis.transform_electronic_integrals(prop)
        problem.properties.electronic_density = ElectronicDensity(
            transformed.alpha, transformed.beta, transformed.beta_alpha
        )

    elif isinstance(prop, AngularMomentum):
        if prop.overlap is None:
            # only the size needs to be changed
            problem.properties.add(prop.__class__(problem.num_spatial_orbitals))
            continue

        if isinstance(transformer._active_basis.coefficients, ElectronicIntegrals):
            coeff_alpha = transformer._active_basis.coefficients.alpha["+-"]
            coeff_beta: Tensor
            if transformer._active_basis.coefficients.beta.is_empty():
                coeff_beta = coeff_alpha
            else:
                coeff_beta = transformer._active_basis.coefficients.beta["+-"]

            problem.properties.angular_momentum = AngularMomentum(
                problem.num_spatial_orbitals,
                coeff_alpha.transpose() @ prop.overlap @ coeff_beta,
            )

norb = problem.num_spatial_orbitals
problem.properties.particle_number = ParticleNumber(norb)
problem.properties.magnetization = Magnetization(norb)
problem.properties.angular_momentum = AngularMomentum(norb)
problem.properties.electronic_density = ElectronicDensity.from_orbital_occupation(
    problem.orbital_occupations, 
    problem.orbital_occupations_b
)

'''
6. mapping the Hamiltonian to spin space
'''
mapper = ParityMapper(num_particles=problem.num_particles)
mapper = problem.get_tapered_mapper(mapper)
Red_hamil_z2qubit = mapper.map(reduced_hamiltonian.second_q_op())
print(f"Number of items in the PM Z2 Pauli list:", len(Red_hamil_z2qubit))
Red_hamil_z2qubit

'''
7. Setup solver
'''

algo = NumPyEigensolver(k=10)
algo.filter_criterion = problem.get_default_filter_criterion()

solver = ExcitedStatesEigensolver(mapper, algo)
result = solver.solve(problem)

print(f"Total ground state energy = {result.total_energies[0]:.4f}")
print(f"Total first excited state energy = {result.total_energies[1]:.3f}")
print(f"Total second excited state energy = {result.total_energies[2]:.3f}")
