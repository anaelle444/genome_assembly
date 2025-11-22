#!/usr/bin/env python3
"""
Script pour comparer les résultats d'assemblage QUAST
Compare différentes tailles de k-mers et l'assemblage Minia
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# Configuration du style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)
plt.rcParams['font.size'] = 10

# Définir les dossiers à analyser
assemblies = {
    '7-mers': '7mers',
    '11-mers': '11mers',
    '21-mers': '21mers',
    '31-mers': '31mers',
    '91-mers': '91mers',
    'Minia': 'minia_results'
}

def parse_quast_report(report_path):
    """Parse un fichier report.tsv de QUAST"""
    data = {}
    try:
        with open(report_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip()
                    # Convertir en nombre si possible
                    try:
                        if '.' in value:
                            data[key] = float(value)
                        else:
                            data[key] = int(value)
                    except (ValueError, AttributeError):
                        # Gérer les valeurs spéciales comme '-' ou texte
                        if value == '-':
                            data[key] = None
                        else:
                            data[key] = value
    except FileNotFoundError:
        print(f"Fichier non trouvé: {report_path}")
        return None
    return data

def collect_all_data(base_path, assemblies):
    """Collecte les données de tous les assemblages"""
    all_data = []
    
    for name, folder in assemblies.items():
        report_path = os.path.join(base_path, folder, 'report.tsv')
        data = parse_quast_report(report_path)
        if data:
            data['Assembly_Name'] = name
            all_data.append(data)
    
    return pd.DataFrame(all_data)

def create_comparison_plots(df, output_dir):
    """Crée les graphiques de comparaison"""
    
    # Créer le dossier de sortie s'il n'existe pas
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Comparaison GC% vs Référence
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: GC content comparison
    ax1 = axes[0, 0]
    x_pos = range(len(df))
    width = 0.35
    
    gc_assembly = df['GC (%)'].values
    gc_ref = df['Reference GC (%)'].fillna(0).values
    
    ax1.bar([x - width/2 for x in x_pos], gc_assembly, width, label='Assemblage', color='skyblue')
    ax1.bar([x + width/2 for x in x_pos], gc_ref, width, label='Référence', color='lightcoral')
    ax1.set_xlabel('Assemblage', fontweight='bold')
    ax1.set_ylabel('GC (%)', fontweight='bold')
    ax1.set_title('Comparaison du contenu GC', fontweight='bold', fontsize=12)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(df['Assembly_Name'], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Genome Fraction
    ax2 = axes[0, 1]
    genome_frac = df['Genome fraction (%)'].fillna(0).values
    colors = ['#2ecc71' if x > 90 else '#f39c12' if x > 50 else '#e74c3c' for x in genome_frac]
    bars = ax2.bar(df['Assembly_Name'], genome_frac, color=colors, alpha=0.7)
    ax2.set_xlabel('Assemblage', fontweight='bold')
    ax2.set_ylabel('Fraction du génome (%)', fontweight='bold')
    ax2.set_title('Couverture du génome de référence', fontweight='bold', fontsize=12)
    ax2.set_xticklabels(df['Assembly_Name'], rotation=45, ha='right')
    ax2.axhline(y=90, color='green', linestyle='--', alpha=0.5, label='Excellent (>90%)')
    ax2.axhline(y=50, color='orange', linestyle='--', alpha=0.5, label='Moyen (>50%)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Ajouter les valeurs sur les barres
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    # Plot 3: N50 et NGA50
    ax3 = axes[1, 0]
    n50_vals = df['N50'].fillna(0).values
    nga50_vals = df['NGA50'].fillna(0).values
    
    x_pos = range(len(df))
    width = 0.35
    ax3.bar([x - width/2 for x in x_pos], n50_vals, width, label='N50', color='steelblue')
    ax3.bar([x + width/2 for x in x_pos], nga50_vals, width, label='NGA50', color='darkorange')
    ax3.set_xlabel('Assemblage', fontweight='bold')
    ax3.set_ylabel('Longueur (bp)', fontweight='bold')
    ax3.set_title('N50 et NGA50 (qualité de contiguïté)', fontweight='bold', fontsize=12)
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(df['Assembly_Name'], rotation=45, ha='right')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Nombre de contigs et erreurs
    ax4 = axes[1, 1]
    ax4_twin = ax4.twinx()
    
    num_contigs = df['# contigs'].values
    misassemblies = df['# misassemblies'].fillna(0).values
    
    x_pos = range(len(df))
    line1 = ax4.plot(x_pos, num_contigs, 'o-', color='blue', linewidth=2, 
                     markersize=8, label='Nombre de contigs')
    line2 = ax4_twin.plot(x_pos, misassemblies, 's-', color='red', linewidth=2,
                          markersize=8, label='Misassemblies')
    
    ax4.set_xlabel('Assemblage', fontweight='bold')
    ax4.set_ylabel('Nombre de contigs', fontweight='bold', color='blue')
    ax4_twin.set_ylabel('Nombre de misassemblies', fontweight='bold', color='red')
    ax4.set_title('Fragmentation et erreurs d\'assemblage', fontweight='bold', fontsize=12)
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(df['Assembly_Name'], rotation=45, ha='right')
    ax4.tick_params(axis='y', labelcolor='blue')
    ax4_twin.tick_params(axis='y', labelcolor='red')
    ax4.grid(True, alpha=0.3)
    
    # Combiner les légendes
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax4.legend(lines, labels, loc='upper left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'comparison_overview.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Graphique de vue d'ensemble sauvegardé: comparison_overview.png")
    plt.close()
    
    # 2. Graphique détaillé des statistiques d'alignement
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Longueur totale vs Référence
    ax1 = axes[0, 0]
    total_len = df['Total length'].values
    ref_len = df['Reference length'].fillna(0).values
    
    x_pos = range(len(df))
    width = 0.35
    ax1.bar([x - width/2 for x in x_pos], total_len, width, label='Longueur assemblage', color='mediumseagreen')
    ax1.bar([x + width/2 for x in x_pos], ref_len, width, label='Longueur référence', color='salmon', alpha=0.7)
    ax1.set_xlabel('Assemblage', fontweight='bold')
    ax1.set_ylabel('Longueur (bp)', fontweight='bold')
    ax1.set_title('Longueur totale vs Référence', fontweight='bold', fontsize=12)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(df['Assembly_Name'], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Plus grand contig
    ax2 = axes[0, 1]
    largest = df['Largest contig'].values
    bars = ax2.bar(df['Assembly_Name'], largest, color='purple', alpha=0.6)
    ax2.set_xlabel('Assemblage', fontweight='bold')
    ax2.set_ylabel('Longueur (bp)', fontweight='bold')
    ax2.set_title('Taille du plus grand contig', fontweight='bold', fontsize=12)
    ax2.set_xticklabels(df['Assembly_Name'], rotation=45, ha='right')
    ax2.grid(True, alpha=0.3)
    
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height):,}', ha='center', va='bottom', fontsize=8)
    
    # Plot 3: Taux d'erreurs
    ax3 = axes[1, 0]
    mismatches = df['# mismatches per 100 kbp'].fillna(0).values
    indels = df['# indels per 100 kbp'].fillna(0).values
    
    x_pos = range(len(df))
    width = 0.35
    ax3.bar([x - width/2 for x in x_pos], mismatches, width, label='Mismatches', color='indianred')
    ax3.bar([x + width/2 for x in x_pos], indels, width, label='Indels', color='lightcoral')
    ax3.set_xlabel('Assemblage', fontweight='bold')
    ax3.set_ylabel('Erreurs par 100 kbp', fontweight='bold')
    ax3.set_title('Taux d\'erreurs (mismatches et indels)', fontweight='bold', fontsize=12)
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(df['Assembly_Name'], rotation=45, ha='right')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Duplication ratio
    ax4 = axes[1, 1]
    dup_ratio = df['Duplication ratio'].fillna(1).values
    colors_dup = ['#27ae60' if x <= 1.05 else '#e67e22' if x <= 1.2 else '#c0392b' for x in dup_ratio]
    bars = ax4.bar(df['Assembly_Name'], dup_ratio, color=colors_dup, alpha=0.7)
    ax4.set_xlabel('Assemblage', fontweight='bold')
    ax4.set_ylabel('Ratio de duplication', fontweight='bold')
    ax4.set_title('Ratio de duplication (idéal ≈ 1.0)', fontweight='bold', fontsize=12)
    ax4.set_xticklabels(df['Assembly_Name'], rotation=45, ha='right')
    ax4.axhline(y=1.0, color='green', linestyle='--', linewidth=2, label='Idéal (1.0)')
    ax4.axhline(y=1.05, color='orange', linestyle='--', alpha=0.5, label='Acceptable (<1.05)')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    for bar in bars:
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'detailed_stats.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Graphique détaillé sauvegardé: detailed_stats.png")
    plt.close()
    
    # 3. Heatmap de toutes les métriques normalisées
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Sélectionner les métriques importantes
    metrics = [
        'Genome fraction (%)',
        'N50',
        'NGA50',
        'Largest contig',
        '# contigs',
        '# misassemblies',
        '# mismatches per 100 kbp',
        '# indels per 100 kbp',
        'Duplication ratio'
    ]
    
    # Créer un DataFrame pour la heatmap
    heatmap_data = df[['Assembly_Name'] + metrics].set_index('Assembly_Name')
    
    # Normaliser les données (0-1) pour chaque métrique
    # Pour les métriques où "plus c'est mieux", normaliser normalement
    # Pour les métriques où "moins c'est mieux", inverser
    normalized_data = heatmap_data.copy()
    
    # Métriques à maximiser (plus = mieux)
    maximize_metrics = ['Genome fraction (%)', 'N50', 'NGA50', 'Largest contig']
    for metric in maximize_metrics:
        if metric in normalized_data.columns:
            max_val = normalized_data[metric].max()
            min_val = normalized_data[metric].min()
            if max_val > min_val:
                normalized_data[metric] = (normalized_data[metric] - min_val) / (max_val - min_val)
            else:
                normalized_data[metric] = 1.0
    
    # Métriques à minimiser (moins = mieux)
    minimize_metrics = ['# contigs', '# misassemblies', '# mismatches per 100 kbp', 
                        '# indels per 100 kbp']
    for metric in minimize_metrics:
        if metric in normalized_data.columns:
            max_val = normalized_data[metric].max()
            min_val = normalized_data[metric].min()
            if max_val > min_val:
                normalized_data[metric] = 1 - (normalized_data[metric] - min_val) / (max_val - min_val)
            else:
                normalized_data[metric] = 1.0
    
    # Duplication ratio (proche de 1 est mieux)
    if 'Duplication ratio' in normalized_data.columns:
        normalized_data['Duplication ratio'] = 1 - abs(normalized_data['Duplication ratio'] - 1)
    
    # Remplacer les NaN par 0
    normalized_data = normalized_data.fillna(0)
    
    # Créer la heatmap
    sns.heatmap(normalized_data.T, annot=True, fmt='.2f', cmap='RdYlGn', 
                cbar_kws={'label': 'Score normalisé (0=pire, 1=meilleur)'},
                linewidths=0.5, ax=ax)
    ax.set_title('Comparaison normalisée des métriques d\'assemblage\n(Vert=Meilleur, Rouge=Moins bon)', 
                 fontweight='bold', fontsize=13, pad=20)
    ax.set_xlabel('Assemblage', fontweight='bold')
    ax.set_ylabel('Métrique', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'heatmap_comparison.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Heatmap de comparaison sauvegardée: heatmap_comparison.png")
    plt.close()
    
    # 4. Graphique radar/spider pour comparaison globale
    from math import pi
    
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
    
    # Utiliser les données normalisées
    categories = list(normalized_data.columns)
    N = len(categories)
    
    # Calculer les angles pour chaque axe
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]
    
    # Tracer chaque assemblage
    colors_radar = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
    
    for idx, (assembly, row) in enumerate(normalized_data.iterrows()):
        values = row.values.tolist()
        values += values[:1]
        ax.plot(angles, values, 'o-', linewidth=2, label=assembly, 
                color=colors_radar[idx % len(colors_radar)])
        ax.fill(angles, values, alpha=0.15, color=colors_radar[idx % len(colors_radar)])
    
    # Configurer les labels
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, size=10)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], size=8)
    ax.set_title('Comparaison globale des assemblages\n(Plus la zone est grande, meilleur est l\'assemblage)', 
                 fontweight='bold', fontsize=12, pad=20)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
    ax.grid(True)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'radar_comparison.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Graphique radar sauvegardé: radar_comparison.png")
    plt.close()

def create_summary_table(df, output_dir):
    """Crée un tableau récapitulatif en HTML et CSV"""
    
    # Sélectionner les colonnes importantes
    summary_cols = [
        'Assembly_Name',
        '# contigs',
        'Total length',
        'Largest contig',
        'N50',
        'NGA50',
        'GC (%)',
        'Genome fraction (%)',
        '# misassemblies',
        '# mismatches per 100 kbp',
        '# indels per 100 kbp',
        'Duplication ratio'
    ]
    
    summary_df = df[summary_cols].copy()
    
    # Arrondir les valeurs
    for col in summary_df.columns:
        if summary_df[col].dtype in ['float64', 'float32']:
            summary_df[col] = summary_df[col].round(2)
    
    # Sauvegarder en CSV
    csv_path = os.path.join(output_dir, 'assembly_comparison_summary.csv')
    summary_df.to_csv(csv_path, index=False, sep='\t')
    print(f"✓ Tableau récapitulatif sauvegardé: assembly_comparison_summary.csv")
    
    # Créer une version HTML stylisée
    html_path = os.path.join(output_dir, 'assembly_comparison_summary.html')
    
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <title>Comparaison des assemblages</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
            }
            h1 {
                color: #2c3e50;
                text-align: center;
            }
            table {
                border-collapse: collapse;
                width: 100%;
                margin: 20px 0;
                background-color: white;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            th {
                background-color: #3498db;
                color: white;
                padding: 12px;
                text-align: left;
                font-weight: bold;
            }
            td {
                padding: 10px;
                border-bottom: 1px solid #ddd;
            }
            tr:hover {
                background-color: #f5f5f5;
            }
            .best {
                background-color: #d4edda;
                font-weight: bold;
            }
            .note {
                margin-top: 20px;
                padding: 15px;
                background-color: #fff3cd;
                border-left: 4px solid #ffc107;
            }
        </style>
    </head>
    <body>
        <h1>📊 Comparaison des assemblages - Résultats QUAST</h1>
        <div class="note">
            <strong>Note:</strong> Les cellules en vert clair indiquent les meilleures valeurs pour chaque métrique.
        </div>
    """
    
    # Convertir le DataFrame en HTML
    html_table = summary_df.to_html(index=False, classes='data', border=0)
    
    # Identifier les meilleures valeurs pour chaque colonne
    for col in summary_cols[1:]:  # Skip Assembly_Name
        if col in ['# contigs', '# misassemblies', '# mismatches per 100 kbp', '# indels per 100 kbp']:
            # Pour ces métriques, la valeur la plus basse est la meilleure
            best_val = summary_df[col].min()
        elif col == 'Duplication ratio':
            # Pour le ratio de duplication, le plus proche de 1.0 est le meilleur
            best_val = summary_df[col].iloc[(summary_df[col] - 1).abs().argsort()[0]]
        else:
            # Pour les autres, la valeur la plus élevée est la meilleure
            best_val = summary_df[col].max()
        
        # Ajouter la classe "best" aux cellules avec les meilleures valeurs
        html_table = html_table.replace(f'<td>{best_val}</td>', f'<td class="best">{best_val}</td>')
    
    html_content += html_table
    html_content += """
        <div class="note" style="background-color: #d1ecf1; border-left-color: #17a2b8; margin-top: 30px;">
            <h3>Interprétation des métriques:</h3>
            <ul>
                <li><strong>Genome fraction (%):</strong> Pourcentage du génome de référence couvert (plus élevé = meilleur)</li>
                <li><strong>N50:</strong> Longueur du contig tel que 50% de l'assemblage est dans des contigs ≥ cette longueur</li>
                <li><strong>NGA50:</strong> N50 aligné sur le génome de référence (plus fiable)</li>
                <li><strong># misassemblies:</strong> Nombre d'erreurs d'assemblage détectées (plus bas = meilleur)</li>
                <li><strong>GC (%):</strong> Contenu en GC de l'assemblage (devrait être proche de la référence: 44.14%)</li>
                <li><strong>Duplication ratio:</strong> Ratio de duplication (idéal ≈ 1.0)</li>
            </ul>
        </div>
    </body>
    </html>
    """
    
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✓ Rapport HTML sauvegardé: assembly_comparison_summary.html")

def main():
    # Chemin de base (utiliser le symlink ASCII-safe)
    base_path = '/home/najat/quast_ascii/resultats'
    output_dir = os.path.join(base_path, 'comparison_plots')
    
    print("=" * 60)
    print("Analyse comparative des assemblages QUAST")
    print("=" * 60)
    print()
    
    # Collecter les données
    print("📊 Collecte des données QUAST...")
    df = collect_all_data(base_path, assemblies)
    
    if df.empty:
        print("❌ Erreur: Aucune donnée n'a pu être collectée!")
        return
    
    print(f"✓ {len(df)} assemblages analysés: {', '.join(df['Assembly_Name'].tolist())}")
    print()
    
    # Créer les visualisations
    print("📈 Génération des graphiques de comparaison...")
    create_comparison_plots(df, output_dir)
    print()
    
    # Créer le tableau récapitulatif
    print("📋 Création du tableau récapitulatif...")
    create_summary_table(df, output_dir)
    print()
    
    print("=" * 60)
    print("✅ Analyse terminée!")
    print(f"📁 Tous les résultats sont dans: {output_dir}")
    print()
    print("Fichiers générés:")
    print("  • comparison_overview.png - Vue d'ensemble des principales métriques")
    print("  • detailed_stats.png - Statistiques détaillées d'alignement")
    print("  • heatmap_comparison.png - Heatmap des métriques normalisées")
    print("  • radar_comparison.png - Graphique radar de comparaison globale")
    print("  • assembly_comparison_summary.csv - Tableau récapitulatif (CSV)")
    print("  • assembly_comparison_summary.html - Rapport HTML interactif")
    print("=" * 60)

if __name__ == "__main__":
    main()
