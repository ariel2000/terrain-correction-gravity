# Atlas visuel des résultats

Cet atlas relie les figures aux questions scientifiques du projet. Les cartes
ne remplacent pas les fichiers numériques : elles servent à vérifier la
cohérence spatiale, repérer les anomalies et comparer les méthodes.

## 1. Zone d'étude et relief hybride

![Zone d'étude et MNT hybride](../outputs/reference/figures/study_area_hybrid_dem.jpg)

L'emprise des stations couvre le Sud-Ouest et l'Ouest du Cameroun, entre
8,502–10,9933° E et 3,5006–6,5000° N. La carte détaillée révèle les édifices
volcaniques, les hauts plateaux, les vallées et les secteurs littoraux que la
surface Surfer lissait fortement. Les stations sont colorées par leur altitude
mesurée, tandis que le fond représente le MNT hybride.

## 2. Comparaison Surfer–hybride

![Comparaison spatiale des MNT](../outputs/reference/figures/dem_spatial_comparison.jpg)

Les deux premiers panneaux utilisent la même échelle altimétrique. Le troisième
représente `MNT hybride − MNT Surfer`. Les différences organisées le long des
crêtes et vallées montrent que l'écart n'est pas un simple décalage vertical :
il provient aussi de la morphologie plus détaillée conservée par FABDEM.

Cette comparaison n'est pas une validation indépendante de Surfer, puisque la
grille Surfer a été interpolée à partir des mêmes stations.

## 3. Validation croisée du MNT hybride

![Validation du MNT hybride](../outputs/reference/figures/hybrid_dem_validation.jpg)

La validation retire des blocs spatiaux de 25 km avant de prédire les résidus.
Le MNT hybride ramène le biais régional de +1,67 m à −0,05 m. Les diminutions de
MAE, RMSE et P95 sont modestes : elles ne doivent pas être exagérées. Le choix du
modèle hybride est surtout justifié par la suppression du biais sans destruction
des formes topographiques satellitaires.

## 4. Comparaison des méthodes gravimétriques

![Comparaison des méthodes](../outputs/reference/figures/method_pairwise_comparison_200km.jpg)

Chaque point représente une station à 200 km. Nagy-prismes est placé en abscisse
comme correction conventionnelle de référence. La masse ponctuelle s'en approche
fortement dans le champ lointain. Hammer est systématiquement plus faible en
raison de la moyenne dans les compartiments. Vannes conserve son signe et sa
définition d'effet direct résiduel ; il ne doit donc pas être interprété comme
une correction positive interchangeable avec Nagy.

## 5. Influence du rayon

![Convergence radiale](../outputs/reference/figures/radial_convergence_detailed.jpg)

Le panneau gauche montre les valeurs cumulées moyennes et P95. Le panneau droit
montre le P95 de l'incrément produit par chaque couronne supplémentaire. Les
courbes restent nettement au-dessus de 0,01 mGal entre 100 et 200 km : la
convergence au seuil choisi n'est pas atteinte pour la majorité des stations.

## 6. Contrôle qualité spatial

![Carte du contrôle qualité](../outputs/reference/figures/station_altitude_qc_map.jpg)

Les 24 croix rouges sont les coordonnées uniques exclues du conditionnement du
MNT, principalement parce que `|altitude station − FABDEM| > 150 m`. Elles sont
conservées dans les calculs gravimétriques et signalées dans les résultats. Leur
concentration sur certains reliefs justifie un contrôle indépendant des
altitudes et du référentiel vertical avant interprétation géologique.

## Régénération

Après téléchargement des dalles indiquées dans `data/external/README.md` :

```powershell
python scripts/make_visual_atlas.py `
  --hybrid-dir data/processed/hybrid_dem/tiles `
  --fabdem-dir data/external/fabdem_v1_2
```

Le script utilise les résultats tabulaires déjà calculés et écrit les figures
dans `outputs/reference/figures/`.
