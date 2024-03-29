---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: qiskit-env
    language: python
    name: qiskit-env
---

# QISKIT Lab 4 - Fermionic Hamiltonian

1. **S. G. Nana Engo**, serge.nana-engo@facsciences-uy1.cm
    * Department of Physics, Faculty of Science, University of Yaounde I
1. **J-P. Tchapet Njafa**, jean-pierre.tchapet-njafa@univ-maroua.cm
    * Department of Physics, Faculty of Science, University of Maroua
1. **P. Djorwe**, djorwepp@gmail.com
    * Department of Physics, Faculty of Science, University of Ngaoundere
       
March 2023


$
\newcommand{\HH}{\mathtt{H}}  
\newcommand{\ad}{a^\dagger}  
$


## Bases chimiques

Un ensemble de base est un ensemble de fonctions, appelées **[fonctions de base](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Quantum_Mechanics/17%3A_Quantum_Calculations/ab_initio_Basis_Sets)**, telles que toute orbitale moléculaire électronique peut être approchée
comme une **combinaison linéaire de ses orbitales atomiques (LCAO, Linear) combination of atomic orbitals**.

Les deux classes d'orbitales de base approximatives couramment utilisées sont les **Slater-tyoe orbitals (STOs)** basées sur le déterminant de Slater, et les **orbitales cartésiennes de type Gaussien (GTO)**. Ces deux types de fonctions de base peut être combiné comme **STO-nG (Slater-type orbital-n Gaussians)**, où n est le nombre de gaussiennes utilisées pour faire les approximations.

Détaillons la structure de ces deux classes.

### Orbitales de type Slater
 Ce sont les fonctions d'état radiales de la forme

 \begin{align*}
& R_\ell(r) = A(\ell,\alpha) r^\ell e^{-\alpha r},
& A(\ell,\alpha) = (2\alpha)^{\ell+1} \sqrt{\frac{2\alpha}{(2\ell+2)!}},
\end{align*}

où,
 * $\ell\leq n$ est le nombre quantique de moment quantique orbital et $n$ le nombre quantique principal;
 * $r$ est la distance entre l'électron et le noyau atomique;
 * $\alpha$ est l'exposant orbital qui contrôle la vitesse à laquelle la densité de l'orbite s'annule en fonction de la distance nucléaire $r$;
 * $A(\ell,\alpha)$ est le facteur de normalisation.

Pour une orbitale $1s$, $\ell=0$ et
\begin{align*}
&A(0,\alpha) = 2\alpha^{3/2}, &R_0(r) =2\alpha^{3/2}e^{-\alpha r}.
\end{align*}

Un code python permettant de visualiser cette fonction est le suivant:

```python vscode={"languageId": "python"}
import numpy as np
import matplotlib.pyplot as plt
```

```python vscode={"languageId": "python"}
x = np.linspace(-5,5,num=1000)
r = abs(x)
alpha = 1.0

R = 2*alpha**(1.5)*np.exp(-alpha*r)

plt.figure(figsize=(4,3))
plt.plot(x,R,label="STO 1s H")
plt.legend()

```

### Do it yourself

Ecriver un code python pour visualiser STO 2s H.


Un code python permettant visualiser la fonction d'état spatiale antisymétrique pour la molécule d'hydrogène comme une combinaison linéaire de la partie radiale précédente de l'orbitale de Slater pour un atome d'hydrogène est la suivante:

```python vscode={"languageId": "python"}
x = np.linspace(-7,7,num=1000)
r1 = abs(x+2.5)
r2 = abs(x-2.5)
alpha = 1.0

R = 2*alpha**(1.5)*np.exp(-alpha*r1)-2*alpha**(1.5)*np.exp(-
alpha*r2)

plt.figure(figsize=(4,3))
plt.plot(x,R,label="Antisymmetric STO H2")
plt.legend()

```

### Orbitales de type GTO

Ce sont les fonctions d'état radiales de la forme

 \begin{equation*}
 R_\ell(r) = B(\ell,\alpha) r^\ell e^{-\alpha r^2},
\end{equation*}

où,
 * $\ell\leq n$ est le nombre quantique de moment quantique orbital et $n$ le nombre quantique principal;
 * $r$ est la distance entre l'électron et le noyau atomique;
 * $\alpha$ est l'exposant orbital qui contrôle la vitesse à laquelle la densité de l'orbite s'annule en fonction de la distance nucléaire $r$;
 * $B(\ell,\alpha)$ est le facteur de normalisation.

En pratique, nous approximons la partie radiale d'un STO avec une combinaison linéaire de fonctions Gaussiennes primitives, appelée **fonction Gaussienne contractée**. Les ensembles de base STO-nG incluent une fonction Gaussienne contractée par orbitale atomique.

Un code python permettant de visualiser la fonction STO-3G pour l'orbite 1 s de l'atome d'hydrogène fonction est le suivant:

```python vscode={"languageId": "python"}
x = np.linspace(-7,7,num=1000)
r = abs(x)
c = [0.444635,0.535328,0.154329]
alpha = [0.109818,0.405771,2.227660]

psi = 0
for k in range(3):
    psi += c[k]*(2*alpha[k]/np.pi)**0.75 * np.exp(-alpha[k]*r**2)

plt.figure(figsize=(5,3))
plt.plot(x,psi,label="STO-3G 1s H")
plt.legend()

```

Un code python permettant de visualiser la fonction d'état spatiale antisymétrique pour la molécule d'hydrogène comme une combinaison linéaire de la partie radiale précédente de la fonction STO-3G pour l'orbite $1s$ d'un atome d'hydrogène est le suivant :

```python vscode={"languageId": "python"}
x = np.linspace(-7,7,num=1000)
r1 = abs(x+2.5)
r2 = abs(x-2.5)
c = [0.444635,0.535328,0.154329]
alpha = [0.109818,0.405771,2.227660]

psi = 0
for k in range(3):
    psi += c[k]*(2*alpha[k]/np.pi)**0.75 * np.exp(-
alpha[k]*r1**2) \
- c[k]*(2*alpha[k]/np.pi)**0.75 * np.exp(-alpha[k]*r2**2)

plt.figure(figsize=(5,3))
plt.plot(x,psi,label="Antisymmetric STO-3G H2")
plt.legend()
```

## Construction d'un opérateur Hamiltonien fermionique 

L'Hamiltonien  est exprimé dans la base des solutions de la méthode HF, également appelées Orbitales Moléculaires (OM) :

$$
\hat{H}_{elec}=\sum_{pq} h_{pq} a^{\dagger}_p a_q + 
\frac{1}{2} \sum_{pqrs} h_{pqrs}  a^{\dagger}_p a^{\dagger}_q a_r  a_s
$$
avec 
* les **intégrales à 1 électron**
$$
h_{pq} = \int \phi^*_p(r) \left( -\frac{1}{2} \nabla^2 - \sum_{I} \frac{Z_I}{R_I- r} \right)   \phi_q(r)dr
$$
qui décrivent l’énergie cinétique des électrons individuels et leurs interactions avec les champs électriques du noyau ;
* et **intégrales à 2 électrons**
$$
h_{pqrs} = \int \frac{\phi^*_p(r_1)  \phi^*_q(r_2) \phi_r(r_2)  \phi_s(r_1)}{|r_1-r_2|}dr_1dr_2,
$$
décrivent les interactions entre les électrons.

Les MO ($\phi_u$) peuvent être occupés ou virtuels (inoccupés). Un MO peut contenir 2 électrons. Cependant, dans ce qui suit, nous travaillons en fait avec des orbitales de spin qui sont associées à un spin up ($\alpha$) d'électron spin down ($\beta$). Ces deux spins sont également communément désignés par $\alpha$ et $\beta$, respectivement. Ainsi, les orbitales de spin peuvent contenir un électron ou être inoccupées.

### Représentation interne

On peut avoir une idée de ce que les termes intégraux signifie en observant les opérateurs de création et d’annihilation qui les composent. Par exemple, $h_{pq} a^{\dagger}_p a_q$ décrit les sauts d’électron de l’orbital de rotation $q$ à l’orbital de rotation $p$. De même, le terme $ h_{pqrs}  a^{\dagger}_p a^{\dagger}_q a_r a_s$ (pour un p,q,r,s distinct) décrit deux électrons dans des orbitaux de rotation $r$ et $s$ se dispersant et se terminant par des orbitaux de rotation $p$ et $q$. Si $r=q$ et $p=s$ alors $ h_{prrp}  a^{\dagger}_p a^{\dagger}_r a_r  a_p = h_{prrp}n_pn_r$ donne la pénalité d’énergie associée aux deux électrons étant proches l’un de l’autre, mais ne décrit pas un processus dynamique.

Dans l'expression de l'Hamiltonien, il existe au maximum $N^2+N^4$ coefficients. Cependant, un grand nombre de ces coefficients peuvent être collectés, car ils correspondent au même opérateur. Par exemple, dans le cas où $p,q,r,s$ sont des indices distincts, on peut utiliser les règles d'anti-commutation pour indiquer que :
$$
\ad_p\ad_qa_ra_s = -\ad_q\ad_pa_ra_s = -\ad_p\ad_qa_sa_r = \ad_q\ad_pa_sa_r.
$$}
En outre, comme $\HH$ est Hermitien, tout opérateur fermionique non Hermitien, par exemple $h_{pqrs}\ad_p\ad_qa_ra_s$, a un conjugué Hermitien qui se trouve également dans $\HH$. Afin d'indexer de manière unique les groupes de termes caractérisés par ces symétries, nous définissons un ordre canonique sur les indices $(i_1,\cdots,i_n,j_1,\cdots,j_m)$ de toute suite de $n+m$ opérateurs fermioniques $ \ad_{i_1}\cdots \ad_{i_n}a_{j_1}\cdots a_{j_m}$ comme suit :

* Tous les opérateurs de création $\ad_{i_\cdot}$ sont placés avant tous les opérateurs d'annihilation $a_{j_\cdot}$.

* Tous les index des opérateurs de création sont triés par ordre croissant, c'est-à-dire $i_1< i_2< \cdots < i_n$.

* Tous les indices des opérateurs d'annihilation sont triés par ordre décroissant, c'est-à-dire $j_1> j_2 \cdots > j_m$.

* L'index le plus à gauche est inférieur ou égal à l'index le plus à droite, c'est-à-dire $i_1\le j_m$.

Identifions cet ensemble d'indices canoniquement ordonnés comme
$$
(i_1,\cdots,i_n,j_1,\cdots,j_m) \in S_{n,m}.
$$
Avec cet ordre canonique, l'Hamiltonien fermionique peut être exprimé comme
$$
\HH=\frac12\sum_{(p,q)\in S_{1,1}}h'_{pq}(\ad_pa_q+\ad_qa_p)
+\frac12\sum_{(p,q,r,s)\in S_{2,2}}h'_{pqrs}(\ad_p\ad_qa_ra_s+\ad_s\ad_ra_qa_p),
$$
avec des intégrales à un et deux électrons convenablement adaptées $h'_{pq}$ et $h'_{pqrs}$, respectivement.



## Qiskit Nature

Ce tutoriel utilise essentiellement le package Qiskit Nature dont la conception abstraite est donnée par la figure ci-dessous.

<center><img src="Graphics/Qiskit_Nature_overview.png" width="500" /></center>

Le package se divise en deux concepts avec chacun trois piliers chacun :

* `problems`,  qui sont des représentations de problèmes scientifiques auxquels on cherche une solution;
    * _Electronic Structure Problem_ représentant le problème de l'équation de Schrödinger électronique des systèmes moléculaires;

    * _Vibrational Structure Problem_ représentant le problème posé par l'Hamiltonien de Watson des systèmes moléculaires;

    * _Lattice Model Problem_ représentant les problèmes définis sur des treillis;

* `algorithms`, qui fournissent les moyens de trouver des solutions auxdits problèmes;

    * _Ground State Solver_ pour trouver l'état fondamental d'un problème;

    * _Excited States Solver_ pour trouver les états excités d'un problème;

    * _Hamiltonian Simulation_ pour simuler la dynamique d'un problème (pas encore implémenté).


### Articulation des modules de Qiskit Nature

La bibliothèque Qiskit Nature comprend différents modules s'articulant autour de :

- chargement de données à partir de pilotes (drivers) de chimie (PySCF, Psi4, Gaussian, etc.) ou de formats de fichiers;
- construction et manipulation d'opérateurs de seconde quantification;
- traduction de la seconde quantification à l'espace qubit;
- une bibliothèque de circuits quantiques d'analyses ciblées en sciences naturelles;
- algorithmes et utilitaires spécifiques aux sciences naturelles pour utiliser  les algorithmes de `Qiskit Terra` plus faciles;
- et beaucoup plus.

Par exemple, les pilotes (drivers) de chimie, lorsqu'ils sont fournis avec une configuration moléculaire, renverront des intégrales à 1 ($h_{pq}$) et 2 ($h_{pqrs}$) corps, ainsi que d'autres données qui sont efficacement calculées de manière classique. Ces données de sortie d'un pilote peuvent ensuite être utilisées comme entrée dans Qiskit Nature qui contient une logique capable de les traduire sous une forme adaptée aux algorithmes quantiques. La conversion crée d'abord un `FermionicOperator` qui doit ensuite être mappé, par ex. par une cartographie de **Jordan Wigner**, à un opérateur qubit prêt pour le calcul quantique.

> pip install qiskit-nature[pyscf] -U

Nous allons utiliser le plugin `qiskit_nature_pyscf` qui couple PySCF et Qiskit Nature.  C'est un solveur [FCI](https://en.wikipedia.org/wiki/Full_configuration_interaction) (Full Configuration Interaction) basé sur Qiskit qui permet à un utilisateur de PySCF (Python-based Simulations of Chemistry Framework) de tirer parti des algorithmes quantique implémentés dans Qiskit pour être utilisés à la place de leurs homologues classiques (dans un esprit similaire à l'intégration NWChemEx). Ce plugin est assez recent (début 2023).
> pip install qiskit-nature-pyscf -U



### Cas de la molécule d'hygrogène

<center><img src="Graphics/h2.png" width="150"/></center>

Dans ce qui suit, nous utilisons le driver PySCF, pouronstruction d'un opérateur Hamiltonien fermioniqu la molécule d'hydrogène à la longueur de la liaison d'équilibre (0,735 angström) à l'état singulet et sans charge.


#### Information sur la structure moléculaire

1. [Qiskit Nature V0.5](https://qiskit.org/documentation/nature/migration/00b_Electronic_structure_with_v0.5.html) a introduit la classe `qiskit_nature.second_q.formats.molecule_info.MoleculeInfo` qui stocke les informations moléculaires sous le format
    ```
    MoleculeInfo(
        symbols: 'Sequence[str]',the ordered sequence of atoms which make up this molecule

        coords: 'Sequence[tuple[float, float, float]]', The XYZ coordinates of the atoms

        multiplicity: 'int' = 1, the multiplicity of the molecule (`= 2 * spin + 1`)

        charge: 'int' = 0, he total charge of the molecule

        units: 'DistanceUnit' = <DistanceUnit.ANGSTROM: 'Angstrom'>, the distance unit in which the XYZ coordinates are stored

        masses: 'Sequence[float] | None' = None, the sequence of masses for all atoms part of the molecule
    ) 
    ```
    Par la suite, `qiskit_nature.second_q.drivers.PySCFDriver.from_molecule` la met sous une forme requise par  `PySCFDriver`:
    ```
    PySCFDriver.from_molecule(
    molecule: 'MoleculeInfo', the molecular information
    
    basis: 'str', the basis set
    method: 'MethodType' = <MethodType.RHF: 'rhf'>, the SCF method type
    )
    ```

 2. Mais, si on utilise le plugin `qiskt_nature_pyscf`, l initialisation de ls structure moléculaire va se faire avec `pyscf.gto.mole` sous le format 
 > gto.M(atom='H 0 0 0; F 0 0 1', basis='6-31g')

Comme ce format étant aussi celui de `qiskit_nature.second_q.drivers.PySCFDriver`, c'est lui nous allons préférentiellement utiliser.

```python vscode={"languageId": "python"}
# Pour les données moléculaires
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
```

```python vscode={"languageId": "python"}
PySCFDriver?
```

```python vscode={"languageId": "python"}
H2_driver = PySCFDriver(
    atom="H 0 0 0; H 0 0 0.735",
    basis="sto3g",
    charge=0,
    spin=0,
    unit=DistanceUnit.ANGSTROM,
)
```

Qiskit fournit une classe utile nommée `ElectronicStructureProblem`, qui appelle le driver (qui contient déjà les informations moléculaires stockées à l'étape précédente) pour construire les orbitales moléculaires. Nous effectuons un calcul HF pour la base 1 et 2 STO-3G.

 Il s'agit en réalité représenter le problème de l'equation de Schrödinguer électronique, $\mathtt{H}_{\rm el}|\Psi\rangle = E_{\rm el}|\Psi\rangle$. $\mathtt{H}_{\rm el}$ est l'Hamiltonien de la classe `qiskit_nature.second_q.hamiltonians.ElectronicEnergy`

```python vscode={"languageId": "python"}
H2_problem = H2_driver.run()
print(H2_problem)
```

### `ElectronicStructureProblem` et ses composants

Examinons cette instance de problème et ses composants.

#### Hamiltonien `ElectronicEnergy`

L'aspect le plus important est l'Hamiltonien interne ; dans ce cas, un Hamiltonien `ElectronicEnergy`. Cette classe est capable de générer un opérateur de seconde quantification à partir des intégrales à 1 ($h_{pq}$) et 2 ($h_{pqrs}$) corps qu'un code classique a calculées pour nous.

> **NB :** La classe de conteneur pour les coefficients intégraux (`PolynomialTensor`) nécessite que les termes à 2 corps soient fournis dans **l'ordre du physicien** !

Ces tenseurs sont soumis à certaines attentes, à savoir :

* pour les intégrales électroniques $\alpha$ (spins-up) et $\beta$ (spins-down), seules les clés suivantes sont autorisées : "", "+-", "++--";
* pour les intégrales électroniques $\alpha\beta$ (spins-down-up) la seule clé autorisée est "++--".


```python vscode={"languageId": "python"}
H2_hamiltonian = H2_problem.hamiltonian

H2_coefficients = H2_hamiltonian.electronic_integrals
print(H2_coefficients.alpha)
```

Nous utilisons la méthode `second_q_ops()`, pour obtenir les integrales $h_{pq}$ et $h_{pqrs}$ en seconde quantication.

```python vscode={"languageId": "python"}
H2_fermionic_op = H2_hamiltonian.second_q_op()
print(H2_fermionic_op)
```

Notez qu'il s'agit purement de l'Hamiltonien **électronique** du système. Cela signifie que l'_énergie de répulsion nucléaire_ n'est pas incluse. Au lieu de cela, Qiskit Nature ajoutera ce décalage d'énergie constant dans une étape de post-traitement, afin de calculer l'énergie totale de votre système. Pour savoir comment inclure l'énergie de répulsion nucléaire dans cet opérateur, veuillez vous référer à la documentation de la classe `ElectronicEnergy` [ici](https://qiskit.org/documentation/nature/stubs/qiskit<_>nature.second<_>q.hamiltonians.ElectronicEnergy.html).

Par exemple que les termes $\alpha$ sont
```
 -1.2563390730032498 * ( +_0 -_0 )
 -0.47189600728114245 * ( +_1 -_1 )
```
et les termes $\beta$ sont
```
 -1.2563390730032498 * ( +_2 -_2 )
 -0.47189600728114245 * ( +_3 -_3 )
```


```python vscode={"languageId": "python"}
H2_hamiltonian.nuclear_repulsion_energy 
```




#### Plus d'attributs du `ElectronicStructureProblem`

Examinons quelques attributs supplémentaires de notre instance "problème".
 <!-- qui renvoie une liste de seconds opérateurs quantifiés : 
* opérateur hamiltonien, 
* opérateur du nombre total de particules, 
* opérateur moment angulaire total,
* opérateur d'aimantation totale 
* et,opérateur du dipôle x, y, z,  si disponible. -->

```python vscode={"languageId": "python"}
H2_problem.molecule
```

```python vscode={"languageId": "python"}
H2_problem.reference_energy
```

```python vscode={"languageId": "python"}
H2_problem.num_particles
```

```python vscode={"languageId": "python"}
H2_problem.num_spatial_orbitals
```

```python vscode={"languageId": "python"}
H2_problem.basis
```

```python vscode={"languageId": "python"}
H2_problem.orbital_energies
```

```python vscode={"languageId": "python"}
H2_problem.num_alpha
```

```python vscode={"languageId": "python"}
H2_problem.num_beta
```

```python vscode={"languageId": "python"}
H2_problem.num_spin_orbitals
```

Définissons la fonction `get_particle_number()` qui permet d'obtenir les propriétés d'une instance `problem` d'une structure électronique.

```python vscode={"languageId": "python"}
def get_particle_number(problem):
    print("Hydrogen molecule, basis: sto3g, Hartree-Fock calculation")
    print(f"Number of alpha electrons: {problem.num_alpha}")
    print(f"Number of beta electrons: {problem.num_beta}")
    print(f"Number of spin orbitals: {problem.num_spin_orbitals}")
 
```

Nous appelons la fonction `get_particle_number()` pour obtenir et imprimer les propriétés du nombre de particules comme suit :

```python vscode={"languageId": "python"}
get_particle_number(H2_problem)
```

Pour comprendre l'expression du Hamiltonien en seconde quantification affichée, il faut comprendre la classe `FermionicOp`.

#### `FermionicOp` Class

Un `FermionicOp` représente une somme pondérée de termes d'opérateurs de création/annihilation fermioniques. Ces termes sont codés sous forme d'étiquettes éparses, qui sont des chaînes constituées d'une liste d'expressions séparées par des espaces. Chaque expression doit ressembler à `[+-]_<index>`, où `<index>` est un entier non négatif représentant l'index du mode fermionique où l'opération `+` (création) ou `-` (annihilation) doit être effectuée. La valeur de l'indice est liée par le nombre d'orbitales de spin (`num_spin_orbitals`) de l'opérateur. Puisque les indices Python commencent par 0, la valeur maximale qu'un indice peut prendre est donnée par `num_spin_orbitals`-1.

Par exemple, `(+_0 -_0)`$\equiv a_0^\dagger a_0$, `( +_2 +_3 -_2 -_3 )`$\equiv a_2^\dagger a_3^\dagger a_2 a_3$.


Imprimons seulement les 10 premiers termes de l'opérateur fermionique de la molécule:

```python vscode={"languageId": "python"}
print("\n".join(str(H2_fermionic_op).splitlines()[:12] + ["..."]))
```

Utilisons la méthode `FermionicOp.to_matrix` pour obtenir une représentation matricielle de la molécule d'hydrogène dans la base de Fock, où les états de base sont classés dans l'ordre croissant des chaines de bits, comme 0000, 0001, ..., 1111.

```python vscode={"languageId": "python"}
print(H2_fermionic_op.to_matrix())
```

Ce résultat indique que l'opérateur Hamiltonien fermionique de la molécule d'hydrogène contient 

1. 4 opérateurs nombre de particules

$h_{pq}$ | $a_p^\dagger a_q$ | +0 | -0 | +1 | -1|  +2 | -2 | +3 |  -3 
---------|-------------------|----|----|----|---|-----|----|----|-----
$-0.4719$| $a_3^\dagger a_3$ |    |    |    |   |     |    |  + | +  |
$-1.2563$| $a_2^\dagger a_2$ |    |    |    |   | +   | +  |    |    |
$-0.4719$| $a_1^\dagger a_1$ |    |    |  + | + |     |    |    |    |
$-1.2563$| $a_0^\dagger a_0$ |  + | +  |    |   |     |    |    |    |

2. 10 opérateurs d'échange entre 2 électrons

$h_{pqrs}$|$a_p^\dagger a_q^\dagger a_r a_s$ | +0 | -0 | +1 | -1|  +2 | -2 | +3 |  -3 
---------|------------------------------|----|----|----|----|---|-----|----|--------- 
$+0.1809$| $a_0^\dagger a_2^\dagger a_3 a_1$ |  + |    |    | + |  +  |    |    | +  |
$-0.1809$| $a_0^\dagger a_3^\dagger a_2 a_1$ |  + |    |    |   |     | +  | +  |    |
$-0.1809$| $a_1^\dagger a_2^\dagger a_3 a_0$ |    | +  |  + |   |  +  |    |    | +  |
$+0.1809$| $a_1^\dagger a_3^\dagger a_2 a_0$ |    | +  |  + |   |     |  + | +  |    |
$+0.4836$| $a_2^\dagger a_3^\dagger a_3 a_2$ |    |    |    |   | +   |  + |  + | +  |
$+0.6986$| $a_1^\dagger a_3^\dagger a_3 a_1$ |    |    | +  | + |     |    |  + | +  |
$+0.6646$| $a_1^\dagger a_2^\dagger a_2 a_1$ |    |    |  + | + |  +  | +  |    |    |
$+0.6646$| $a_0^\dagger a_3^\dagger a_3 a_0$ |  + | +  |    |   |     |    |  + | +  |
$+0.6757$| $a_0^\dagger a_2^\dagger a_2 a_0$ | +  | +  |    |   |  +  | +  |    |    |
$+0.4836$| $a_0^\dagger a_1^\dagger a_1 a_0$ | +  | *  | +  | + |     |    |    |    |



### Utilisation du plugin `qiskiy_nature_pyscf`

```python vscode={"languageId": "python"}
from pyscf import gto, scf, mcscf
from pyscf.mcscf import avas #AVAS method to construct mcscf active space

from qiskit_nature_pyscf import QiskitSolver
```

- Initialisons la structure moléculaire

```python vscode={"languageId": "python"}
H2_mol = gto.M(atom="H 0 0 0; H 0 0 .735", basis="sto-3g")
```

- Effectuons les calculs HF avec `pyscf.scf.RHF` (Restricted).

```python vscode={"languageId": "python"}
H2_h_f = scf.RHF(H2_mol).run()
```

- Effectuons les calculs post-HF avec `pyscf.mcscf.CASCI` (Complete active space configuration interaction) pour améliorer les solutions de l'équation de Schrödinger précédent.

```python vscode={"languageId": "python"}
# To obtain norb, the number of (active) orbitals and nelec, the number of (active) electrons. 
norb, nel, mo =avas.avas(H2_h_f,['Li 2s','H 1s'])

H2_cas = mcscf.CASCI(H2_h_f, norb, nelec)

H2_cas
```

Il ne reste plus qu'à intégrer un algorithme quantique à notre simulation! Nous allons donc, lors du prochain tutorial ou lab session, élaborer cet algorithme.


## Cas de la molécule d'hydride de lithium

<center><img src="Graphics/Lithium_hydride.png" width="150"/></center>

Do it yourself

```
LiH_molecule = Molecule(geometry=[['Li', [0., 0., 0.]], ['H', [0., 0., 1.5474]]], charge=0, multiplicity=1)
```


```python vscode={"languageId": "python"}
import qiskit.tools.jupyter

%qiskit_version_table
```
