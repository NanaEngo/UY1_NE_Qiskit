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
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

## 

<!-- #region -->
# UE 4268 - QISKIT LAB 4 - Evalution formative


**`vos noms et prenoms`**, `matricule` et `email` 


Department of Physics - Faculty of Science - University of Yaoundé I

`Nom du Laboratoire`

`Date`

 **Duree: 120 min**
<!-- #endregion -->

$$
\newcommand{\br}{\mathbf{r}}  
\newcommand{\ad}{a^\dagger}  
$$

<!-- #region -->
## 1 - Atome d'Helium



Considérons l’exemple d’un système de deux électrons situés dans le champ coulombien d’un
noyau $Ze_{0}$. C’est le cas de l’atome d’hélium.


 **1. Rappeler l’expression générale du Hamiltonien de cet atome.**

<!-- #endregion -->



<!-- #region -->

**2. Espace symétrique et antisymétrique des orbitales de type STO.**

 
   _Soit $\phi_{1S}(r)= 2\alpha^{3/2}e^{-\alpha r}$ la fonction d'etat d'un electron a la position r._

 (a) Ecrire les fonctions d'etat de type Slater symétriques et antisymétriques normalisées de He.
<!-- #endregion -->




 (b) Ecrire un code python permettant de visualiser ces fonctions.


```python
import numpy as np
import matplotlib.pyplot as plt
```

```python
#Reponse
```

### **3. Espace symétrique et antisymétrique des orbitales de type GTO**

  _Soit $\phi_{1S}(r)= (\frac{2\alpha}{\pi})^{\frac{3}{4}}e^{-\alpha r^2}$ la fonction d'etat d'un electron a la position r._

 (a) Ecrire les fonctions d'etat de type Gaussienne symétriques et antisymétriques normalisées de He.
 

```python
#Reponse
```

(b) Ecrire un code python permettant de visualiser ces fonctions.


```python
#Reponse
```

### 2 - Relations d'anti-commutation

```python
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.operators.commutators import anti_commutator
```

   Verifier les egalités suivantes :
   
   * $\{a_1, \ad_1\}= 1$
   
   * $\{\ad_2, \ad_3\}= 0$
   
   * $\{\ad_2, \ad_3\}= 0$
   * $\ad_0a_1\ad_1a_2 = \ad_0a_2 - \ad_0\ad_1a_1\ad_2$
   
   * $\ad_0a_4+\ad_0a_2=0$

```python
#Reponse
```

<!-- #region -->
## 3 - Problème de structure électronique : Cas de la molécule d’eau (H2O)



  On donne les coordonnées xyz de la molécule
  
  **H : 0.968877 0.012358 0.000000**
  
  **O : -0.019830 -0.025588 0.000000**
  
  **H : -0.229801 0.941311 0.000000**


<!-- #endregion -->

```python
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
```


  * Reprendre les étapes du problème électronique du cas de la molécule d’hydrogène pour
traiter le cas de la molécule d’eau.


```python
#reponse
```

```python

```

* Imprimer les 6 premiers termes de l'opérateur fermionique de la molécule $H_2O$:

```python
#reponse
```

* Utiliser la méthode `FermionicOp.to_matrix` pour obtenir une représentation matricielle de la molécule d'eau dans la base de Fock.

```python
#reponse
```

```python

```


Lien utile : https://github.com/qiskit-community/qiskit-community-tutorials/blob/master/chemistry/QubitMappings.ipynb?short_path=70b1d78

```python
import qiskit.tools.jupyter

%qiskit_version_table
```
