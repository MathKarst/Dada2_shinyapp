# DADA2 Shiny

*[English version](README.md)*

Interface Shiny pas-à-pas pour le pipeline DADA2 (Illumina, single- ou paired-end),
pensée pour plusieurs marqueurs (16S, 18S, rbcL, COI, ITS).

## Installation

```r
install.packages(c("shiny", "DT", "ggplot2", "future", "promises", "BiocManager"))
BiocManager::install(c("dada2", "phyloseq", "Biostrings", "decontam", "DECIPHER"))
```

cutadapt (recommandé pour retirer les amorces) : `pip install cutadapt` ou
`conda install -c bioconda cutadapt`. Sous Windows, passer par conda ou WSL et
indiquer le chemin de l'exécutable dans l'onglet Amorces.

## Lancement

```r
shiny::runApp("chemin/vers/dada2shiny")   # ou ouvrir app.R dans RStudio > Run App
```

Pour déboguer (tout dans le processus principal, sans calcul en arrière-plan) :
`options(dada2shiny.async = FALSE)` avant `runApp()`.

## Organisation

```
dada2shiny/
├── app.R                  interface Shiny
└── pipeline/
    ├── pipeline.R         fonctions DADA2 indépendantes de Shiny (utilisables en script)
    ├── viz.R              abondances relatives et Krona (utilisables en script)
    └── presets.R          préréglages par marqueur et bases de référence
```

## Étapes

1. **Données** : détection et appariement R1/R2 *par nom d'échantillon* ; multi-runs
   (un sous-dossier par run) ; single-end ; métadonnées CSV ; reprise d'une session.
2. **Qualité** : profils de qualité.
3. **Amorces** : préréglages (16S V4, V3-V4, 18S V9, rbcL diatomées, COI Leray, ITS2),
   contrôle de la présence des amorces dans chaque orientation, retrait par cutadapt
   (mélanges d'amorces et IUPAC acceptés) ou par longueur fixe.
4. **Filtrage** : `filterAndTrim`, avec calcul du chevauchement attendu et alertes
   (truncLen plus long que les reads, chevauchement insuffisant, reads qui dépassent l'amplicon).
5. **Modèle d'erreur** : appris *run par run* ; option « qualités binnées » pour NovaSeq/NextSeq.
6. **Débruitage et fusion** : run par run, puis `mergeSequenceTables` ; alerte si < 50 % de fusion.
7. **Longueur et chimères**.
8. **Contaminants et filtres** (optionnel) : decontam (prévalence/fréquence), seuil de tag-jump,
   minimum de reads par ASV et par échantillon.
9. **Taxonomie** : `assignTaxonomy` (FASTA DADA2) ou IDTAXA (trainingSet .RData) ; six bases
   téléchargées et mises en cache automatiquement (voir plus bas) ; tout autre fichier en local
   ou en upload.

Puis : **Suivi des reads**, **Export** (table, FASTA, phyloseq, script R reproductible,
session) et **Exploration** (voir ci-dessous).

## Exploration

L'objet phyloseq (table finale + taxonomie + métadonnées) est construit automatiquement à
l'ouverture de l'onglet ; les ASV y portent les mêmes identifiants que dans les exports.

- **Abondances relatives** : barres empilées ou heatmap, au niveau des ASV (sans taxonomie,
  avec option d'annoter chaque ASV par son genre) ou agrégées à n'importe quel rang. Les taxons
  les moins abondants sont regroupés dans « Autres », les non assignés dans « Non assigné ».
  Échantillons moyennés par une variable de métadonnées, facettes, ordre par nom, par
  similarité (Bray-Curtis, UPGMA) ou par taxon dominant. Export PDF, PNG et CSV.
- **Krona** : graphique interactif généré directement en R (KronaTools n'est pas nécessaire),
  un jeu de données par échantillon, par groupe ou cumulé, avec ou sans les ASV au dernier
  niveau. Téléchargeable en HTML. Le moteur Krona est chargé depuis marbl.github.io : il faut
  une connexion internet pour l'afficher.
- **Diversité alpha**, **NMDS Bray-Curtis**, et l'objet phyloseq téléchargeable en RDS.

## Tout lancer

Le bouton **Tout lancer** (bandeau du haut, visible depuis tous les onglets) enchaîne les
étapes avec les réglages actuels des onglets :

- il part de la première étape non faite ou dont les paramètres ont changé (pastille orange),
  ou des amorces si « Tout recalculer » est coché ;
- il s'arrête à l'étape choisie dans « jusqu'à » ;
- il s'arrête à la première erreur ; les étapes déjà terminées restent sauvegardées ;
- la taxonomie est sautée si aucune référence n'est indiquée ;
- « Arrêter après l'étape en cours » interrompt l'enchaînement proprement.

Sur un nouveau jeu de données, mieux vaut faire les premières étapes pas à pas (profils de
qualité, modèle d'erreur, taux de fusion), puis utiliser Tout lancer une fois les réglages validés.

## Amorces et préréglages

- **Fichier d'amorces** (onglet Amorces) :
  - FASTA : le sens est lu dans le nom (`>rbcL_708F_1`, `>R3_2`, `>515F-Y`, `>806RB`, `>COI_fwd`…) ;
  - CSV/TSV : une colonne `sequence` et une colonne `direction` (F/R, forward/reverse).
    C'est la solution à privilégier quand les noms ne portent pas le sens (ex. `fITS7`, `ITS4`).
- **Mes préréglages** : « Enregistrer mes réglages » sauvegarde tous les réglages des onglets
  Amorces à Taxonomie (amorces, troncature, maxEE, modèle d'erreur, chimères, filtres,
  référence…) sous un nom. Ils apparaissent ensuite dans la liste des préréglages, sous
  « Mes préréglages ». Ils sont stockés hors du dossier de l'appli (dossier de configuration
  de R, `tools::R_user_dir("dada2shiny", "config")`) et survivent donc aux mises à jour.

## Bases de référence

Téléchargées et mises en cache automatiquement : SILVA 138.2 (16S), RDP 19 (16S),
Greengenes2 2024.09 (16S), PR2 5.0.0 (18S), SILVA 138.2 Fungi 18S et 26S.

À télécharger soi-même puis utiliser via « Fichier local » (liens dans l'onglet Taxonomie) :
Diat.barcode (rbcL), UNITE (ITS), MIDORI2 (COI et autres gènes mitochondriaux).

## Fonctionnement

- Les étapes longues tournent dans un processus séparé : l'interface reste utilisable et le
  **Journal** se met à jour en direct.
- Relancer une étape efface les résultats des étapes suivantes. La barre de statut indique
  les étapes faites (vert), en cours (bleu) et celles dont les paramètres ont changé depuis
  l'exécution (orange).
- Après chaque étape, le dossier de sortie contient `dada2_session.rds` (pour reprendre le
  projet) et `dada2_log.txt` (le journal).
- Le script R exporté utilise les paramètres *réellement utilisés* à chaque étape et tourne
  sans l'appli.

## À vérifier

- Les préréglages sont des points de départ pour du 2x250/2x300 : à ajuster selon les profils
  de qualité. Les séquences d'amorces sont à vérifier contre les publications d'origine,
  en particulier le mélange rbcL (Vasselon et al. 2017).
- COI : il n'existe pas de base BOLD officielle au format DADA2 (voir MIDORI2 ou une base maison).
  rbcL : Diat.barcode au format DADA2, en adaptant les rangs à ce fichier.
- Pas adapté aux reads longs Nanopore/PacBio sans adaptation.

## Citer

Si tu utilises cet outil, cite aussi les méthodes qu'il emploie : DADA2 (Callahan et al. 2016,
*Nature Methods*), cutadapt (Martin 2011, *EMBnet.journal*), decontam (Davis et al. 2018,
*Microbiome*), phyloseq (McMurdie & Holmes 2013, *PLoS ONE*), Krona (Ondov et al. 2011,
*BMC Bioinformatics*), IDTAXA (Murali et al. 2018, *Microbiome*) le cas échéant, ainsi que la
base de référence utilisée.

## Licence

MIT, voir [LICENSE](LICENSE).
