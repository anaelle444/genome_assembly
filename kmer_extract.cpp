//
//  kmer_extract.cpp
//  
//
//  Created by Anaelle Ji-Seun Joo on 14/11/2025.
//

#include "kmer_extract.hpp"
#include <algorithm>
#include <vector>
#include <string>

using namespace std;

// Trie une liste de k-mers par ordre alphabétique
std::vector<string> trier(std::vector<string> unsorted_list) {
    std::sort(unsorted_list.begin(), unsorted_list.end());
    return unsorted_list;
}

// KmerExtract
// Entrée :
//   k : un entier (taille des k-mers)
//   F : un ensemble de mots avec un alphabet A de taille n
// Sortie : une liste triée de k-mers issus des mots de F
std::vector<string> kmerExtract(int k, const std::vector<string>& F) {
    std::vector<string> L;  // Liste de k-mers
    
    // Pour chaque séquence dans F
    for (size_t i = 0; i < F.size(); i++) {
        // Pour chaque position possible dans la séquence
        for (size_t j = 0; j <= F[i].length() - k; j++) {
            // Extraire le k-mer de position j à j+k-1
            string kmer = F[i].substr(j, k);
            L.push_back(kmer);
        }
    }
    
    // Trier la liste de k-mers
    L = trier(L);
    
    return L;
}
