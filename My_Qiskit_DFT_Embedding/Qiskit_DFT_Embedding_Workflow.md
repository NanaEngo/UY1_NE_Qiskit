# DFT Embedding with Qiskit-nature workflow

1. Effectuer les calculs mean-field (les calculs de référence)

2. Obtenir le Hamiltonien complet du système hors de l'instance `ElectronicStructureProblem`, à partir des intégrales à 1- et 2-corps

3. Définir un espace actif en utilisant `ActiveSpaceTransformer` et en déduire l'Hamiltonien réduit et l'instance `ElectronicStructureProblem` associé

4. Effectuer le mapping dans l'espace qubit

5. Effectuer le calculs des états excités à partir du solver choisi.

De facon plus détaillée on a

1. **Reference Calculation**
    `driver.run_pyscf()` effectue le calcul Hartree-Fock ou DFT, afin d'obtenir les données nécessaires à la construction de l'Hamiltonien.

2. **Hamiltonian Construction**
    `hamiltonian = ElectronicEnergy(...)` crée un objet `ElectronicEnergy`, représentant l'Hamiltonien électronique.

    `ElectronicIntegrals(...)` est utilisé pour construire les intégrales, cad les termes à 1-électron (+-) et à 2-électrons (++--).

    `S8Integrals(...)` peut-être utilisé pour calculer les intégrales à 2-électrons sous une forme symétrique (format S8).

3. **Active Space Preparation**
    `active_orbitals` qui spécifie l'ensemble des orbitales qui seront prises en compte dans l'espace actif.

    `ActiveSpaceTransformer` initialise avec le nombre d'orbitales spatiales et d'électrons dans l'espace actif.

    `transformer.prepare_active_space(...)` configure le transformateur en fonction des orbitales actives sélectionnées.

4. **Hamiltonian Reduction**
    `reduced_hamiltonian = transformer.transform_hamiltonian(hamiltonian)` transforme l'Hamiltonien complet du système en un Hamiltonien de l'espace actif, réduisant ainsi la taille du problème.

5. **ElectronicStructureProblem**

    `ElectronicStructureProblem` crée une instance Problem avec `reduced_hamiltonian`.

    Renseigner, pour le Problem, d'autres propriétés pertinentes (comme les énergies orbitales et l'énergie de référence) .

6. **Qubit Mapping and Reduction**
    `ParityMapper` est utilisé pour mapper les opérateurs fermioniques aux opérateurs qubit et on effectue la réduction $\mathbb{Z}_2$.

7. **Excited States Calculation**
    `NumPyEigensolver` est utilisé pour résoudre l'Hamiltonien qubit afin d'obtenir les états excités.

    `ExcitedStatesEigensolver` est configuré pour utiliser `NumPyEigensolver`.

## ElectronicEnergy calculation

La classe `ElectronicEnergy` implémente le Hamiltonien

$$ H_{el} = \sum_{p, q} h_{pq} a^\dagger_p a_q
         + \sum_{p, q, r, s} g_{pqrs} a^\dagger_p a^\dagger_q a_r a_s ,$$

où $h_{pq}$ et $g_{pqrs}$ sont les intégrales électroniques à un et deux corps,
stocké dans le conteneur `qiskit_nature.second_q.operators.ElectronicIntegrals`.
Lorsqu'il s'agit de coefficients séparés pour les électrons de spin $\alpha$- et $\beta$,
l'Hamiltonien à spin sans restriction (UHF) peut être obtenu à partir de celui ci-dessus de manière simple.

On peut construire une instance de cet Hamiltonien de plusieurs manières.

1. Avec une instance existante de `qiskit_nature.second_q.operators.ElectronicIntegrals` :

``` python
   intégrales : ElectronicIntegrals = ...

   hamiltonien = ElectronicEnergy (intégrales, constantes = {"nuclear_repulsion_energy": 1.0})
```
2. À partir d’un ensemble brut de matrices de coefficients intégraux
    * h1_a: the alpha-spin one-body coefficients.
    * h2_aa: the alpha-alpha-spin two-body coefficients.
    * h1_b: the beta-spin one-body coefficients.
    * h2_bb: the beta-beta-spin two-body coefficients.
    * h2_ba: the beta-alpha-spin two-body coefficients.:

``` python
    # en supposant que les intégrales à un et deux corps ont été préalablement calculées
    h1_a, h2_aa, h1_b, h2_bb, h2_ba = ...

    hamiltonien = ElectronicEnergy.from_raw_integrals(h1_a, h2_aa, h1_b, h2_bb, h2_ba)
    hamiltonian.nuclear_repulsion_energy = 1.0
```

Il est à noter que nous avons spécifié l'énergie de répulsion nucléaire comme un décalage d'énergie constant dans ce qui précède. Ce terme ne sera pas inclus dans l'opérateur qubit mappé puisqu'il s'agit d'un terme de décalage constant et qu'il ne nécessite pas d'erreurs lors de la mesure sur un appareil quantique. Il est cependant possible d'inclure des termes d'énergie constante à l'intérieur du conteneur `qiskit_nature.second_q.operators.ElectronicIntegrals`, si l'on souhait qu'il soit inclus dans l'opérateur qubit, lors du mapping de l'opérateur de 2e quantification à l'espace qubit.
 ``` python
    from qiskit_nature.second_q.operators import PolynomialTensor

        e_nuc = hamiltonian.nuclear_repulsion_energy
        hamiltonian.electronic_integrals.alpha += PolynomialTensor({"": e_nuc})
        hamiltonian.nuclear_repulsion_energy = None
```

## Calcul des matrices de coefficients intégraux (voir def to_qcschema() de pysfdriver)

```python
import numpy as np
from pyscf import ao2mo
from qiskit_nature.second_q.operators.symmetric_two_body import fold

mo_coeff, mo_coeff_b = driver._expand_mo_object(mf.mo_coeff, array_dimension=3)
h1_a = mf.get_hcore()
h1_a = np.dot(np.dot(mo_coeff.T, h1_a), mo_coeff)
if mo_coeff_b is not None:
    h1_b = np.dot(np.dot(mo_coeff_b.T, h1_a), mo_coeff_b)
else:
    h1_b = None


h2_aa = driver._mol.intor("int2e", aosym=8)
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
```

## Calcul simplifié de l'hamiltonian à partir d'une instance QCSchema

```python
from qiskit_nature.second_q.formats.qcschema import QCSchema
from qiskit_nature.second_q.formats import qcschema_translator

qcschema = driver.to_qcschema(include_dipole=True)
hamiltonian = qcschema_translator._get_mo_hamiltonian_direct(qcschema)
```

