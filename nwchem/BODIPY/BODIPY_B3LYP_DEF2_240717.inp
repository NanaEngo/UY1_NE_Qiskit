import mlatom as ml

# Get the initial guess for the molecules to optimize
initmol = ml.data.molecule.from_xyz_string('''21

B           -0.00019437583225       -1.16316739355963       -0.00176466041749
N            1.24631725164681       -0.24138805849423       -0.00080783267345
C            2.51739477533339       -0.61997176132243       -0.00075140337776
C            3.35689094927452        0.50356665438456       -0.00082246625734
C            2.54079373059925        1.61117059666300       -0.00120562482099
C            1.20575157674696        1.13778760226353       -0.00108087155394
C            0.00058911568866        1.82086825837058       -0.00124861634612
C           -1.20493021235489        1.13842232496486       -0.00080649652061
N           -1.24621945325567       -0.24073087404503       -0.00010064286439
C           -2.51749640395899       -0.61864315229661        0.00071302235781
C           -3.35639906343006        0.50533841214423        0.00057932193959
C           -2.53972124748802        1.61251062336190       -0.00048409441898
F           -0.00007497888277       -1.95161629895753        1.12633462988970
F           -0.00072588381522       -1.94780131847676       -1.13265938071381
H            2.78162898425872       -1.66019437319538       -0.00056421851775
H            4.42858745196628        0.48215630454488       -0.00072457564702
H            2.83148118049264        2.64333386128906       -0.00135193565684
H            0.00087029015745        2.90047753858275       -0.00170385373193
H           -2.78227557332534       -1.65872624312978        0.00143771204265
H           -4.42810657767131        0.48449521810039        0.00115450061887
H           -2.82986760835016        2.64482665500762       -0.00091044163019
''')

# define DFT method
b3lyp = ml.models.methods(method='b3lyp/def2-tzvp')

# Optimize the geometry with the choosen optimizer:
geomopt = ml.optimize_geometry(model=b3lyp, initial_molecule=initmol, program='geometric')

# Get the final geometry
final_mol = geomopt.optimized_molecule

# Do vibration analysis and thermochemistry calculation
freq = ml.thermochemistry(model=b3lyp, molecule=final_mol, program='pyscf')
# or vibration analysis only
# freq = ml.freq(model=mymodel, molecule=final_mol, program='')

# Save the molecule with vibration analysis and thermochemistry results
final_mol.dump(filename='final_mol.json',format='json')

# Check vibration analysis
print("Mode     Frequencies     Reduced masses     Force Constants")
print("           (cm^-1)            (AMU)           (mDyne/A)")
for ii in range(len(final_mol.frequencies)):
    print("%d   %13.4f   %13.4f   %13.4f"%(ii,final_mol.frequencies[ii],final_mol.reduced_masses[ii],final_mol.force_constants[ii]))

# Check thermochemistry results
print(f"Zero-point vibrational energy: {final_mol.ZPE} Hartree")
print(f"Enthalpy at 298 K: {final_mol.H} Hartree")
print(f"Gibbs Free energy at 298 K: {final_mol.G} Hartree")
print(f"Heat of formation at 298 K: {final_mol.DeltaHf298} Hartree")
