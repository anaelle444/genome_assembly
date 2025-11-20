//
//  kmer_extract.hpp
//  
//
//  Created by Anaelle Ji-Seun Joo on 14/11/2025.
//

#ifndef kmer_extract_hpp
#define kmer_extract_hpp

#include <string>
#include <vector>

// Trie une liste de k-mers par ordre alphabétique
std::vector<std::string> trier(std::vector<std::string> unsorted_list);

// Extrait tous les k-mers d'un ensemble de séquences
// k : taille des k-mers
// F : ensemble de séquences (reads)
// Retourne : liste triée de k-mers
std::vector<std::string> kmerExtract(int k, const std::vector<std::string>& F);

#endif /* kmer_extract_hpp */
