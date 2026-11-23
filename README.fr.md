# DADA2 Shiny

*[English version](README.md)*

Interface Shiny pas-à-pas pour le pipeline DADA2 (Illumina, single- ou paired-end),
pensée pour plusieurs marqueurs (16S, 18S, rbcL, COI, ITS).

Le dépôt contient deux applis qui partagent les mêmes modules R (`pipeline/`) :
- **`app.R`** : amplicons Illumina avec DADA2 ;
- **`nanopore/app.R`** : amplicons Oxford Nanopore, du FASTQ à la table de séquences,
  entièrement en R, donc utilisable sous Windows (voir [Appli nanopore](#appli-nanopore)).

## Récupérer le code

Avec git :

```
git clone https://github.com/MathKarst/Dada2_shinyapp.git
```

Sans git : sur la [page du dépôt](https://github.com/MathKarst/Dada2_shinyapp), cliquer sur
**Code → Download ZIP**, puis décompresser. Le dossier s'appelle alors `Dada2_shinyapp-main`
(on peut le renommer). Mises à jour : `git pull` dans le dossier, ou retélécharger le ZIP.

## Installation

```r
install.packages(c("shiny", "DT", "ggplot2", "future", "promises", "curl", "BiocManager"))
BiocManager::install(c("dada2", "phyloseq", "Biostrings", "pwalign", "decontam", "DECIPHER", "ShortRead"))
```

`phyloseq`, `Biostrings`, `decontam` et `DECIPHER` sont facultatifs pour l'appli Illumina :
les fonctions correspondantes sont désactivées s'ils manquent. L'appli nanopore a besoin de
`DECIPHER` et `ShortRead`. `pwalign` contient les fonctions d'alignement par paires que les
versions récentes de Bioconductor (>= 3.19) ont sorties de `Biostrings` ; avec une version plus
ancienne, l'appli utilise `Biostrings`.

cutadapt (recommandé pour retirer les amorces) : `pip install cutadapt` ou
`conda install -c bioconda cutadapt`. Sous Windows, passer par conda ou WSL et
indiquer le chemin de l'exécutable dans l'onglet Amorces.

## Lancement

```r
shiny::runApp("chemin/vers/Dada2_shinyapp")   # ou ouvrir app.R dans RStudio > Run App
```

Pour déboguer (tout dans le processus principal, sans calcul en arrière-plan) :
`options(dada2shiny.async = FALSE)` avant `runApp()`.

## Organisation

```
Dada2_shinyapp/
├── app.R                  appli Illumina / DADA2
├── nanopore/app.R         appli Oxford Nanopore
└── pipeline/
    ├── pipeline.R         fonctions DADA2 indépendantes de Shiny (utilisables en script)
    ├── viz.R              abondances relatives et Krona (utilisables en script)
    ├── blast.R            BLAST à distance sur les serveurs du NCBI (utilisable en script)
    ├── curation.R         réplicats, filtre taxonomique, raréfaction, mock, préparation de base
    ├── nanopore.R         nanopore : QC, orientation des amorces, regroupement, consensus, polissage
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
   minimum de reads par ASV et par échantillon, réplicats PCR (détection des PCR ratées et fusion).
9. **Taxonomie** : `assignTaxonomy` (FASTA DADA2) ou IDTAXA (trainingSet .RData) ; six bases
   téléchargées et mises en cache automatiquement (voir plus bas) ; tout autre fichier en local
   ou en upload.
10. **Filtre taxonomique** (optionnel) : retirer ou ne garder que certains taxons, retirer les ASV
    non assignées à un rang donné.

Puis : **Suivi des reads**, **Export** (table, FASTA, phyloseq, script R reproductible,
session) et **Exploration** (voir ci-dessous).

## Réplicats PCR

Ajoute aux métadonnées une colonne indiquant l'échantillon biologique de chaque PCR (par ex.
`echantillon`) et choisis-la dans l'étape Contaminants. Pour chaque échantillon, la composition
de chaque PCR est comparée au centroïde de ses réplicats (Bray-Curtis) ; une PCR est jugée ratée
quand cette dissimilarité dépasse la variation habituelle entre réplicats. Le seuil est fixé
automatiquement là où se croisent les distributions de dissimilarité intra- et inter-échantillons
(comme dans metabaR), ou à la main. La PCR la plus éloignée est retirée en premier, puis les
centroïdes sont recalculés : une seule PCR ratée n'entraîne pas ses bons réplicats ; avec deux
réplicats seulement, les deux sont signalés car on ne peut pas savoir lequel est mauvais. Les
réplicats sont ensuite fusionnés. Les identifiants d'ASV restent les mêmes tout du long (numéro
dans la table sans chimères).

## Communauté mock

Choisis le ou les échantillons mock et charge les séquences attendues (FASTA, une par espèce ;
elles peuvent être plus longues que l'amplicon et dans n'importe quel sens). Chaque ASV est
alignée sur la séquence attendue la plus proche et classée : exacte, variante d'erreur (jusqu'à N
différences) ou inattendue.

## Outil : base de référence

À partir d'un FASTA, avec la taxonomie dans les en-têtes (avec ou sans identifiant au début,
séparateurs `;`, `|` ou `,`) ou dans un tableau à part (une colonne par rang, ou une seule
colonne de taxonomie). Nettoyage : majuscules, U→T, retrait des gaps, des caractères invalides,
des séquences trop courtes et des préfixes comme `k__`, troncature aux rangs vides ou
« unknown ». Doublons exacts retirés ; une même séquence avec des taxonomies différentes est
ramenée à leur ancêtre commun (ou retirée). Les noms rattachés à plusieurs parents sont signalés.
Fichiers produits (dans `refs/` du dossier de sortie) : un fichier `assignTaxonomy`, un fichier
`addSpecies` s'il y a un rang Species, et en option un classificateur IDTAXA (`LearnTaxa`). Un
bouton les envoie dans l'onglet Taxonomie. Note : `assignTaxonomy` n'a pas d'étape
d'entraînement séparée (il apprend à partir du FASTA à chaque lancement) ; IDTAXA s'entraîne une
fois et fait en général moins de sur-classification.

## BLAST

L'onglet **BLAST** envoie les ASV sélectionnées au BLAST du NCBI (interface « Common URL API »),
en une seule recherche de 100 séquences au maximum. Réglages : base (`core_nt` ou `nt`),
megablast ou blastn, nombre de hits, restriction Entrez optionnelle (ex.
`Bacillariophyta[Organism]`), exclusion des séquences « uncultured » / environnementales, et ton
adresse e-mail (demandée par le NCBI). La recherche tourne en arrière-plan (de 30 s à quelques
minutes). Pour chaque ASV : meilleur hit avec identité et couverture, espèces ex aequo (même
score) et assignation DADA2 pour comparaison ; chaque hit renvoie vers sa fiche GenBank. Les
résultats sont gardés dans la session, ajoutés à l'export de la table d'ASV et téléchargeables
en CSV. Les règles d'usage du NCBI sont respectées (au moins 10 s entre deux requêtes, au plus
une vérification du statut par minute).

## Exploration

L'objet phyloseq (table finale + taxonomie + métadonnées) est construit automatiquement à
l'ouverture de l'onglet ; les ASV y portent les mêmes identifiants que dans les exports.

- **Abondances relatives** : barres empilées ou heatmap, au niveau des ASV (sans taxonomie,
  avec option d'annoter chaque ASV par son genre) ou agrégées à n'importe quel rang. Les taxons
  les moins abondants sont regroupés dans « Autres », les non assignés dans « Non assigné ».
  Échantillons moyennés par une variable de métadonnées, facettes, ordre par nom, par
  similarité (Bray-Curtis, UPGMA) ou par taxon dominant. Export PDF, PNG et CSV.
- **Raréfaction** : courbes de raréfaction (richesse attendue, formule exacte) avec la profondeur
  choisie, et diversité moyenne sur des raréfactions répétées (Observed, Shannon, Simpson avec
  écarts-types ; Bray-Curtis moyen pour la NMDS).
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

## Appli nanopore

```r
shiny::runApp("chemin/vers/Dada2_shinyapp/nanopore")   # ou ouvrir nanopore/app.R > Run App
```

Entrée : les FASTQ basecallés (MinKNOW / Dorado), un par code-barres ou un sous-dossier par
code-barres (`fastq_pass/barcode01/…`). Étapes :

1. **Qualité** : longueur et qualité moyenne des reads par échantillon (N50, histogrammes).
2. **Filtre et amorces** : filtre de longueur et de qualité moyenne ; amorces (IUPAC, mélanges,
   préréglages) cherchées aux deux extrémités avec une tolérance d'erreur, reads remis dans le
   même sens et découpés.
3. **Consensus**, par échantillon : regroupement des reads (`DECIPHER::Clusterize`), séparation
   des groupes qui cachent deux variantes proches (positions où deux allèles co-varient),
   alignement (`AlignSeqs`) et consensus « poli » :
   - votes pondérés par la qualité et équilibrés entre brins (les erreurs nanopore dépendent
     souvent du brin) ;
   - correction de la longueur des homopolymères (longueur la plus probable d'après tous les reads) ;
   - positions couvertes par moins de la moitié des reads (restes d'amorce ou d'adaptateur) ignorées.
   Deux modes : *metabarcoding* (tous les groupes d'au moins N reads, pour garder les taxons rares)
   et *barcoding* (groupes majoritaires seulement). En option, expérimental : **Dorado polish**
   (`dorado aligner` + `dorado polish`) si dorado est installé (nécessite Rsamtools).
4. **Table** : les consensus quasi identiques sont fusionnés entre échantillons ; chimères retirées
   avec DADA2 (méthode « pooled », avec une différence tolérée, car le consensus d'une chimère rare
   repose sur peu de reads).
5. **Taxonomie** et filtre taxonomique, évaluation du **mock**, abondances relatives, Krona,
   exports (table, FASTA, phyloseq, script R reproductible).

Pas besoin de Python, de Linux ni de WSL (contrairement à NGSpeciesID / Medaka). Sur des données
simulées proches du R10 (≈ 2 % d'erreurs, surtout des indels dans les homopolymères ; 5 espèces
dont une rare et deux à 2 % l'une de l'autre ; reads dans les deux sens, bruit et chimères), toutes
les vraies séquences ont été retrouvées sans aucune différence, sur plusieurs tirages aléatoires,
et tous les consensus parasites étaient des chimères retirées à l'étape Table. C'est une
vérification sur données simulées : valide sur ta propre communauté mock avant de t'y fier.

Un jeu de test simulé (avec les séquences attendues) est dans `nanopore/exemple/`.

Dépendances : `DECIPHER`, `ShortRead`, `Biostrings`, `dada2` (Bioconductor) ; optionnels
`Rsamtools` (Dorado polish) et `phyloseq`.

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
*Nature Methods*), metabaR pour l'approche des réplicats (Zinger et al. 2021, *Methods in Ecology
and Evolution*), cutadapt (Martin 2011, *EMBnet.journal*), decontam (Davis et al. 2018,
*Microbiome*), phyloseq (McMurdie & Holmes 2013, *PLoS ONE*), Krona (Ondov et al. 2011,
*BMC Bioinformatics*), IDTAXA (Murali et al. 2018, *Microbiome*) le cas échéant, ainsi que la
base de référence utilisée.

## Licence

MIT, voir [LICENSE](LICENSE).
