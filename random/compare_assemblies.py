#!/usr/bin/env python3
"""
Script pour comparer les résultats QUAST entre l'assemblage 51mers et Minia
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Configuration des graphiques
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

# Données extraites des rapports QUAST
data = {
    'Métrique': [
        '# contigs',
        'Largest contig',
        'Total length',
        'Genome fraction (%)',
        '# misassemblies',
        '# local misassemblies',
        'Unaligned length',
        '# mismatches per 100 kbp',
        '# indels per 100 kbp',
        'N50',
        'NA50',
        'Duplication ratio'
    ],
    '51mers (votre assembleur)': [
        1, 16159, 16159, 88.413, 0, 5, 6757, 2329.29, 521.17, 16159, 9402, 1.001
    ],
    'Minia': [
        1, 9936, 9936, 93.524, 0, 0, 0, 0.00, 0.00, 9936, 9936, 1.000
    ]
}

df = pd.DataFrame(data)

# Référence
reference_length = 10624

# Création d'une figure avec plusieurs sous-graphiques
fig = plt.figure(figsize=(16, 12))

# 1. Comparaison des longueurs et couverture
ax1 = plt.subplot(3, 3, 1)
metrics = ['Total length', 'Largest contig', 'N50', 'NA50']
idx = df[df['Métrique'].isin(metrics)].index
x = np.arange(len(metrics))
width = 0.35

values_51 = df.loc[idx, '51mers (votre assembleur)'].values
values_minia = df.loc[idx, 'Minia'].values

ax1.bar(x - width/2, values_51, width, label='51mers', alpha=0.8, color='steelblue')
ax1.bar(x + width/2, values_minia, width, label='Minia', alpha=0.8, color='coral')
ax1.axhline(y=reference_length, color='green', linestyle='--', label='Référence', linewidth=2)
ax1.set_ylabel('Longueur (bp)')
ax1.set_title('Comparaison des longueurs', fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(metrics, rotation=45, ha='right')
ax1.legend()
ax1.grid(axis='y', alpha=0.3)

# 2. Couverture du génome
ax2 = plt.subplot(3, 3, 2)
genome_coverage = [df.loc[3, '51mers (votre assembleur)'], df.loc[3, 'Minia']]
colors = ['steelblue', 'coral']
bars = ax2.bar(['51mers', 'Minia'], genome_coverage, color=colors, alpha=0.8)
ax2.axhline(y=100, color='green', linestyle='--', label='Couverture complète', linewidth=2)
ax2.set_ylabel('Couverture (%)')
ax2.set_title('Fraction du génome assemblée', fontweight='bold')
ax2.set_ylim([0, 105])
for i, bar in enumerate(bars):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
             f'{height:.2f}%', ha='center', va='bottom', fontweight='bold')
ax2.legend()
ax2.grid(axis='y', alpha=0.3)

# 3. Nombre de contigs
ax3 = plt.subplot(3, 3, 3)
n_contigs = [df.loc[0, '51mers (votre assembleur)'], df.loc[0, 'Minia']]
ax3.bar(['51mers', 'Minia'], n_contigs, color=colors, alpha=0.8)
ax3.set_ylabel('Nombre de contigs')
ax3.set_title('Nombre de contigs', fontweight='bold')
ax3.set_ylim([0, max(n_contigs) + 0.5])
for i, v in enumerate(n_contigs):
    ax3.text(i, v, str(int(v)), ha='center', va='bottom', fontweight='bold')
ax3.grid(axis='y', alpha=0.3)

# 4. Misassemblies
ax4 = plt.subplot(3, 3, 4)
misassemblies = [df.loc[4, '51mers (votre assembleur)'], df.loc[4, 'Minia']]
local_mis = [df.loc[5, '51mers (votre assembleur)'], df.loc[5, 'Minia']]
x_pos = np.arange(2)
width = 0.35
ax4.bar(x_pos - width/2, misassemblies, width, label='Misassemblies', alpha=0.8, color='red')
ax4.bar(x_pos + width/2, local_mis, width, label='Local misassemblies', alpha=0.8, color='orange')
ax4.set_ylabel('Nombre')
ax4.set_title('Erreurs d\'assemblage', fontweight='bold')
ax4.set_xticks(x_pos)
ax4.set_xticklabels(['51mers', 'Minia'])
ax4.legend()
ax4.grid(axis='y', alpha=0.3)

# 5. Séquences non alignées
ax5 = plt.subplot(3, 3, 5)
unaligned = [df.loc[6, '51mers (votre assembleur)'], df.loc[6, 'Minia']]
bars = ax5.bar(['51mers', 'Minia'], unaligned, color=['red', 'green'], alpha=0.8)
ax5.set_ylabel('Longueur non alignée (bp)')
ax5.set_title('Séquences non alignées', fontweight='bold')
for i, bar in enumerate(bars):
    height = bar.get_height()
    ax5.text(bar.get_x() + bar.get_width()/2., height,
             f'{int(height)}', ha='center', va='bottom', fontweight='bold')
ax5.grid(axis='y', alpha=0.3)

# 6. Taux d'erreurs
ax6 = plt.subplot(3, 3, 6)
mismatches = [df.loc[7, '51mers (votre assembleur)'], df.loc[7, 'Minia']]
indels = [df.loc[8, '51mers (votre assembleur)'], df.loc[8, 'Minia']]
x_pos = np.arange(2)
width = 0.35
ax6.bar(x_pos - width/2, mismatches, width, label='Mismatches', alpha=0.8, color='purple')
ax6.bar(x_pos + width/2, indels, width, label='Indels', alpha=0.8, color='pink')
ax6.set_ylabel('Erreurs par 100 kbp')
ax6.set_title('Taux d\'erreurs', fontweight='bold')
ax6.set_xticks(x_pos)
ax6.set_xticklabels(['51mers', 'Minia'])
ax6.legend()
ax6.grid(axis='y', alpha=0.3)
ax6.set_yscale('log')

# 7. Ratio de duplication
ax7 = plt.subplot(3, 3, 7)
duplication = [df.loc[11, '51mers (votre assembleur)'], df.loc[11, 'Minia']]
bars = ax7.bar(['51mers', 'Minia'], duplication, color=colors, alpha=0.8)
ax7.axhline(y=1.0, color='green', linestyle='--', label='Pas de duplication', linewidth=2)
ax7.set_ylabel('Ratio de duplication')
ax7.set_title('Duplication du génome', fontweight='bold')
ax7.set_ylim([0.99, 1.01])
for i, bar in enumerate(bars):
    height = bar.get_height()
    ax7.text(bar.get_x() + bar.get_width()/2., height,
             f'{height:.3f}', ha='center', va='bottom', fontweight='bold')
ax7.legend()
ax7.grid(axis='y', alpha=0.3)

# 8. Tableau récapitulatif des métriques clés
ax8 = plt.subplot(3, 3, 8)
ax8.axis('tight')
ax8.axis('off')

summary_metrics = [
    ['Métrique', '51mers', 'Minia', 'Meilleur'],
    ['Genome fraction (%)', '88.41', '93.52', 'Minia'],
    ['# misassemblies', '0', '0', 'Égalité'],
    ['# local misassemblies', '5', '0', 'Minia'],
    ['Mismatches/100kbp', '2329.29', '0.00', 'Minia'],
    ['Indels/100kbp', '521.17', '0.00', 'Minia'],
    ['Unaligned length', '6757', '0', 'Minia'],
    ['N50', '16159', '9936', '51mers'],
]

table = ax8.table(cellText=summary_metrics, cellLoc='center', loc='center',
                  colWidths=[0.35, 0.2, 0.2, 0.25])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 2)

# Colorer l'en-tête
for i in range(4):
    table[(0, i)].set_facecolor('#4472C4')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Colorer les cellules "Meilleur"
for i in range(1, len(summary_metrics)):
    if summary_metrics[i][3] == 'Minia':
        table[(i, 3)].set_facecolor('#90EE90')
    elif summary_metrics[i][3] == '51mers':
        table[(i, 3)].set_facecolor('#87CEEB')
    else:
        table[(i, 3)].set_facecolor('#FFFFE0')

ax8.set_title('Tableau récapitulatif', fontweight='bold', pad=20)

# 9. Score global (à gauche)
ax9 = plt.subplot(3, 3, 9)
ax9.axis('off')

# Calcul d'un score qualitatif
score_51 = f"""
ÉVALUATION QUALITATIVE:

51mers (votre assembleur):
+ Bon: Pas de misassemblies majeures
+ Bon: N50 élevé (16159 bp)
+ Bon: Un seul contig
- Problème: 88.41% du génome couvert
- Problème: 6757 bp non alignés
- Problème: Taux d'erreurs modéré
  (2329.29 mismatches/100kbp)
  (521.17 indels/100kbp)
- Problème: 5 misassemblies locales

Minia:
+ Excellent: 93.52% du génome couvert
+ Excellent: Pas d'erreurs détectées
+ Excellent: Pas de séquences non alignées
+ Bon: Un seul contig
~ N50 légèrement plus faible (9936 bp)

CONCLUSION: Minia produit un assemblage
de meilleure qualité avec moins d'erreurs
et une meilleure couverture.
"""

ax9.text(0.1, 0.5, score_51, fontsize=9, verticalalignment='center',
         family='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('comparison_51mers_vs_minia.png', dpi=300, bbox_inches='tight')
print("✓ Graphique de comparaison sauvegardé: comparison_51mers_vs_minia.png")

# Créer un rapport détaillé en texte
report = f"""
==================================================================================
RAPPORT DE COMPARAISON - ANALYSE QUAST
Assemblage 51mers vs Minia
==================================================================================

GÉNOME DE RÉFÉRENCE:
- Longueur: {reference_length} bp
- GC%: 44.14%

==================================================================================
1. STATISTIQUES GÉNÉRALES
==================================================================================

Métrique                          | 51mers    | Minia     | Meilleur
----------------------------------|-----------|-----------|----------
Nombre de contigs                 | 1         | 1         | Égalité
Longueur totale (bp)              | 16,159    | 9,936     | Minia (plus proche de ref)
Plus grand contig (bp)            | 16,159    | 9,936     | 51mers
Fraction du génome assemblée (%)  | 88.41     | 93.52     | Minia ✓✓
GC content (%)                    | 44.58     | 44.05     | Minia (plus proche de 44.14%)

ANALYSE: Minia couvre 93.52% du génome contre 88.41% pour 51mers.
Votre assembleur produit un contig plus long mais cela inclut des erreurs.

==================================================================================
2. QUALITÉ DE L'ASSEMBLAGE (Métriques d'alignement)
==================================================================================

Métrique                          | 51mers    | Minia     | Meilleur
----------------------------------|-----------|-----------|----------
# Misassemblies                   | 0         | 0         | Égalité ✓
# Local misassemblies             | 5         | 0         | Minia ✓✓
Longueur non alignée (bp)         | 6,757     | 0         | Minia ✓✓
Plus grand alignement (bp)        | 9,402     | 9,936     | Minia ✓
Longueur totale alignée (bp)      | 9,402     | 9,936     | Minia ✓

ANALYSE CRITIQUE:
- Votre assembleur a 5 misassemblies locales (erreurs mineures d'assemblage)
- 6,757 bp ne s'alignent pas sur la référence (41.8% de votre assemblage!)
- Minia a un alignement parfait avec 0 erreur

==================================================================================
3. PRÉCISION (Taux d'erreurs)
==================================================================================

Métrique                          | 51mers    | Minia     | Meilleur
----------------------------------|-----------|-----------|----------
Mismatches par 100 kbp            | 2,329.29  | 0.00      | Minia ✓✓✓
Indels par 100 kbp                | 521.17    | 0.00      | Minia ✓✓✓
Ratio de duplication              | 1.001     | 1.000     | Minia ✓

ANALYSE CRITIQUE:
- Votre assembleur a un TAUX D'ERREUR ÉLEVÉ:
  * ~2.3% de mismatches (environ 1 erreur tous les 43 bp!)
  * ~0.52% d'indels
- Minia a une précision PARFAITE (0 erreur)
- Ces erreurs expliquent pourquoi votre contig est plus long

==================================================================================
4. MÉTRIQUES DE CONTIGUÏTÉ
==================================================================================

Métrique                          | 51mers    | Minia     | Commentaire
----------------------------------|-----------|-----------|-------------
N50                               | 16,159    | 9,936     | 51mers plus élevé
NG50                              | 16,159    | 9,936     | 51mers plus élevé
NA50 (aligned)                    | 9,402     | 9,936     | Minia ✓
NGA50 (aligned)                   | 9,402     | 9,936     | Minia ✓

ANALYSE:
- Le N50 de 51mers est artificiellement gonflé par les erreurs
- Quand on ne considère QUE les séquences correctement alignées (NA50),
  Minia est légèrement meilleur (9,936 vs 9,402)
- NA50 et NGA50 sont les métriques les plus importantes car elles
  excluent les erreurs d'assemblage

==================================================================================
5. SYNTHÈSE ET CONCLUSION
==================================================================================

POINTS FORTS de votre assembleur (51mers):
✓ Structure simple (1 seul contig)
✓ Pas de misassemblies majeures
✓ N50 élevé (mais voir ci-dessous)
✓ Meilleure couverture qu'avec 31mers (88.41% vs 80.89%)

POINTS FAIBLES de votre assembleur (51mers):
✗ Couverture insuffisante: seulement 88.41% du génome
✗ Taux d'erreurs ÉLEVÉ (2329.29 mismatches/100kbp)
✗ 6,757 bp non alignés (séquences erronées ou non présentes dans référence)
✗ 5 misassemblies locales
✗ Précision médiocre comparée à Minia

POINTS FORTS de Minia:
✓✓✓ PRÉCISION PARFAITE (0 erreur)
✓✓ Excellente couverture (93.52% du génome)
✓✓ Pas de séquences non alignées
✓✓ Pas de misassemblies
✓ Structure simple (1 seul contig)

CONCLUSION GÉNÉRALE:
--------------------
Votre assembleur avec k=51 produit un assemblage de qualité MOYENNE comparé
à Minia. Les principales problématiques sont:

1. COUVERTURE: Vous manquez ~11.6% du génome (vs 6.5% pour Minia)
   Amélioration notable par rapport à k=31 (88.41% vs 80.89%)

2. PRÉCISION: Votre taux d'erreurs reste élevé, probablement dû à:
   - Des k-mers erronés non filtrés (erreurs de séquençage)
   - Des répétitions mal gérées
   - Des bulles dans le graphe de De Bruijn non résolues
   - Un chemin eulérien qui traverse des arêtes erronées
   NOTE: k=51 réduit significativement les erreurs vs k=31
         (2329 vs 6407 mismatches/100kbp)

3. SÉQUENCES ERRONÉES: 6,757 bp ne s'alignent pas sur la référence,
   ce qui représente 41.8% de votre assemblage!

RECOMMANDATIONS POUR AMÉLIORER:
-------------------------------
1. Filtrer les k-mers de faible couverture (probablement des erreurs)
2. Implémenter une meilleure résolution des bulles dans le graphe
3. Tester différentes valeurs de k (21, 31, 51, 71, 91)
4. Implémenter une correction d'erreurs avant l'assemblage
5. Améliorer l'algorithme de choix du chemin eulérien

==================================================================================
ÉVALUATION FINALE:
==================================================================================

Assembleur 51mers: 5/10 - Qualité moyenne, amélioration vs 31mers
Minia:             10/10 - Excellente qualité, assemblage de référence

Minia est supérieur en termes de qualité et précision.
Votre assembleur fonctionne et k=51 améliore les résultats,
mais des améliorations restent nécessaires.

==================================================================================
"""

# Sauvegarder le rapport
with open('rapport_comparaison_quast.txt', 'w', encoding='utf-8') as f:
    f.write(report)

print("✓ Rapport détaillé sauvegardé: rapport_comparaison_quast.txt")

# Afficher le rapport
print("\n" + report)

# Créer un tableau pandas pour export
df_export = pd.DataFrame({
    'Métrique': data['Métrique'],
    '51mers': data['51mers (votre assembleur)'],
    'Minia': data['Minia']
})

df_export.to_csv('comparison_table.csv', index=False)
print("✓ Tableau de comparaison sauvegardé: comparison_table.csv")

plt.show()
