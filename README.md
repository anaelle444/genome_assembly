# Genome_assembly

brought to you by an amazing team : Najat, Tara and Anaelle! 


Objectif de TP : Création d'un outil d'assemblage de génome sur un jeu de données de tailles réduites 

Méthode Possible : 
* Méthode Gloutonne 
* Méthode OLC
* Méthode Eulérienne 

=> Méthode Eulérienne avec Graphe de Bruijn et graphe de k-mers 

![texte alternatif](image/schema.png)



DeBruijnAssembler
=================

Petit assembleur basé sur un graphe de De Bruijn et la recherche de chemins eulériens
implémenté en C++.

Usage
-----

Compilation (dans le répertoire `genome_assembly`):

    make

Exécution:

    ./assembler reads.fastq.fq k out.fa

où `k` est la taille du k-mer (entier >= 2). Le programme lit un fichier FASTQ simple
et produit un fichier FASTA avec les contigs reconstruits.

Principes
---------

- On extrait tous les k-mers des reads.
- Les nœuds du graphe sont des (k-1)-mers; chaque k-mer ajoute une arête dirigée
  du préfixe (k-1) vers le suffixe (k-1).
- On applique l'algorithme d'Hierholzer pour trouver des chemins eulériens par composant
  connexe, puis on reconstruit des contigs en réassemblant ces (k-1)-mers.

Limites et améliorations possibles
---------------------------------

- Pas de nettoyage/filtrage des k-mers de faible couverture.
- Pas de compactage (unitigs) avancé ni de suppression d'embranchements.
- Pas d'utilisation de structures mémoire-compactes — adapté aux petits jeux de données.
# Genome_assembly

brought to you by an amazing team : Najat, Tara and Anaelle! 


