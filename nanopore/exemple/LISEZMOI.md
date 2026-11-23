# Jeu d'exemple simulé

Reads simulés proches du nanopore R10 (≈ 2 % d'erreurs, surtout des indels dans les
homopolymères), amplicon rbcL de 400 pb (amorces Diat_rbcL_708F / R3), reads dans les deux sens,
3 % de bruit et 2 % de chimères.

- `reads/barcode01.fastq.gz` : mélange de 5 espèces (dont spD rare, 3 %, et spE à 2 % de spA)
- `reads/barcode02.fastq.gz` : mélange de 3 espèces
- `reads/barcode03.fastq.gz` : un seul spécimen (spA), pour le mode barcoding
- `mock_attendu.fa` : les 5 vraies séquences (onglet « Communauté mock »)
- `reference_test.fa` : petite base au format DADA2 (onglet Taxonomie, « Fichier local »)
- `metadata.csv` : métadonnées d'exemple
- `simulation.py` : le script qui a produit ces données

Essai rapide : dossier des FASTQ = `exemple/reads`, décocher « Un sous-dossier par
code-barres », préréglage « rbcL diatomées », filtre longueur 300-700, puis « Tout lancer ».
Résultat attendu : 5 séquences, toutes identiques aux séquences de `mock_attendu.fa`.
