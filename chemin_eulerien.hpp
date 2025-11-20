//
//  chemin_eulerien.hpp
//  
//
//  Created by Anaelle Ji-Seun Joo on 20/11/2025.
//

#ifndef chemin_eulerien_hpp
#define chemin_eulerien_hpp

#include "graphe_bruijn.hpp"
#include <string>
#include <vector>

// Trouve un chemin eulérien et assemble la séquence
// T : graphe orienté de De Bruijn eulérien
// L : liste des k-mers
// k : taille des k-mers
// Retourne : séquence assemblée
std::string cheminEulerienEtAssemblage(GrapheBruijn& T, 
                                       const std::vector<std::string>& L, 
                                       int k);

#endif /* chemin_eulerien_hpp */
