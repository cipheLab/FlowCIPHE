# FlowCIPHE

Librayrie R du CIPHE contenant des fonction pour le traitement de fichier FCS. La majorité des fonctions tourne autour de flowCore et utilise ces fonctions

## Liste des fonctions 

- read.FCS.CIPHE : charge un fichier FCS dans R en objet flowFrame ou return NULL en cas d'érreur

- catch.create.FCS.from.CSV : créer un flowFrame a partir d'une table (csv) ou return NULL en cas d'erreur.

- found.spill.CIPHE : Lit dans un flowFrame la partie description et cherche une matrice carré contenant en nom de colonne des parametres existant. On présuppose qu'il s'agit alors de la matrice de comp et la fonction retourne le nom du keyword et le nombre de dimension.

- compensate.CIPHE : reprend la fonction compensate de flowCore mais retourne le flowFrame et affiche un warning en cas d'absence de matrice de compensation (au lieu de générer un erreur).

- delete.column.FCS.CIPHE : Suprime un marqueur, une dimension, un anticorps d'un flowFrame (utilisé dans la deconcatenation).

- enrich.FCS.CIPHE : Ajoute un marqueur, une dimension dans un flowFrame (utilisé dans la concatenation).

- concatenante.CIPHE : Créer un flowFrame à partir d'une liste de flowFrame en ajoutant un nouveau parametre indiquand le numero du FCS utile pour deconcatener le fichier.

- deconcatenate.CIPHE : en cours décriture....

- clean.tails.FCS.CIPHE : permet de supprimer les valeur négatives (une intensité n'est pas censé être negative) ou celle qui sont "stacker" a la valeur canal max. A faire avant transformation !!




## How to install
```
install.packages("flowCore")
install.packages("devtools")
library("devtools")
install_github("cipheLab/FlowCIPHE")

```
