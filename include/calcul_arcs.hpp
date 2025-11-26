//
//  calcul_arcs.hpp
//  
//
//  Created by Anaelle Ji-Seun Joo on 20/11/2025.
//

#ifndef calcul_arcs_hpp
#define calcul_arcs_hpp

#include <string>
#include <vector>
#include <utility>  // Pour std::pair

// Recherche un k-mer dans une liste triée
// kmer : le k-mer à rechercher
// L : liste triée de k-mers
// Retourne : l'indice du k-mer dans L, ou -1 s'il n'est pas trouvé
int rechercher(const std::string& kmer, const std::vector<std::string>& L);

// Calcule les arcs du graphe de De Bruijn
// L : liste de k-mers (ordre alphabétique)
// k : taille des k-mers
// Retourne : liste des arcs (paires d'indices des k-mers)
std::vector<std::pair<int, int>> calculArcs(const std::vector<std::string>& L, int k);

#endif /* calcul_arcs_hpp */
