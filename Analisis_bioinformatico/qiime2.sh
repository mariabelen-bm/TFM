######################################################################################################################
### ANÁLISIS BIOINFORMÁTICO DE MUESTRAS DE MOCO DE INTESTINO DE JUVENILES DE TRUCHA ARCOIRIS (ONCORHYNCHUS MYKISS) ###
######################################################################################################################

# TRABAJO DE FIN DE MÁSTER
# MÁSTER UNIVERSITARIO EN BIOINFORMÁTICA
# FACULTAD DE BIOLOGÍA - UNIVERSIDAD DE MURCIA

### AUTOR: María Belén Barquero Martínez
### AÑO: 2022

######################################################################################################################
######################################################################################################################

#################
### OBJETIVO ###
#################

# Determinar si se producen cambios en las poblaciones de microorganismos presentes en el moco del intestino
# de juveniles de trucha arcoiris tras hacerles ingerir dietas experimentales suplementadas con metionina al 1 (control),
# 2 y 3 % respectivamente (Grupos MET1, MET2 y MET3). 
# De producirse cambios sustanciales en el moco del intestino de estos especímenes, identificar si dichos
# cambios traen consigo una posible mejora en la respuesta ante patógenos en los individuos que determine así que 
# la suplementación con metionina en las dietas mejora la respuesta inmunitaria de esta especie.

######################################################################################################################
######################################################################################################################

#######################
### CÓDIGO EMPLEADO ###
#######################

## Inicialización del programa Qiime2 en local mediante ambiente conda

conda activate qiime2-2022.2 

## 1 - IMPORTACIÓN DE LAS MUESTRAS
# El archivo pe-33-manifest contiene la id de la muestra, así como la ruta hasta sus secuencias. 
# En este caso contamos con dos archivos por muestras: R1 (forward) y su R2 (reverse) al ser PairEnd.
# Se indica en este archivo la ubicación exacta en la que se encuentran ambas secuencias .fastq por muestra. 

# Mediante "tools import" se importan los datos de la secuenciación de las muestras en formato "artefacto" con extensión .qz.
# Dicha extensión .qz es con la que trabaja el programa QIIME 2.

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path pe-33-manifest.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

## 1.1 - Visualización de las muestras importadas y demultiplexadas

# Al tener el archivo pe-33-manifest obtendremos un archivo -qza demultiplexado (paired-end-demux.qza), ya que las muestras ya están 
# asociadas a su R1 y R2, por lo que no es necesario realizar un proceso de demultiplexión.

# Se puede visualizar un resumen de los resultados demultiplexados. Esto permite determinar cuántas secuencias
# se obtuvieron por muestra y obtener un resumen de la distribución de las calidades de las secuencias 
# en cada posición de los datos de la secuencia.

qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv

## 2 - ELIMINACIÓN DE LOS PRIMERS Y AJUSTE DE PARÁMETROS
# Se usa el comando cutadapt para eliminar los primers de las secuencias y ajustar algunos parámetros.
# p-front-f es el primer de la secuencia R1 (forward), mientras que p-front-r el de la R2 (reverse).
# match-adapter-wildcards sirve para interpretar los comodines presentes en los adaptadores de las secuencias.

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f CCTACGGGNGGCNGCAG \
  --p-front-r GACTACNNGGGTATCTAATCC \
  --p-indels \
  --p-match-adapter-wildcards \
  --o-trimmed-sequences paired-end-trimmed.qza

## 2.1 - Visualización de las secuencias tras la eliminación de los primers
# Con p-n modificamos el número de secuencias que deben ser seleccionadas al azar para la generación de los gráficos de calidad. 
# En este caso se selecciona que este número sea de 10000.

qiime demux summarize \
  --i-data paired-end-trimmed.qza \
  --p-n 10000 \
  --o-visualization paired-end-trimmed.qzv

# Si se compara la visualización del archivo paired-end-demux.qzv con paired-end-trimmed.qzv se podrá ver que se han eliminado correctamente los primers.

## 3 - ELIMINACIÓN DE RUIDO DE LOS DATOS
# Para el proceso de eliminación de ruido o denoising se empleará el pipeline dada2, que detecta y corrige 
# los datos de secuencias de amplicones de Illumina.
# En este proceso de eliminación de ruido o denoising la longitud para truncar las secuencias tras prueba y error 
# se ha ajustado a 275 bases para el caso de las secuencias forward y 225 para el caso de las reverse, ambas a la derecha.
# Con p-max-ee-f conseguimos descartar las lecturas de avance con un número de errores superior al valor asignado, en este caso 4.

# Se obtendrá un artefacto de salidas con las estadísticas resultantes al aplicar dicho proceso de eliminación de ruido (stats-dada2.qza) y el resultado
# de las secuencias representativas tras el proceso (rep-seqs-dada2.qza).

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-trimmed.qza \
  --p-trunc-len-f 275 \
  --p-trunc-len-r 225 \
  --p-max-ee-f 4 \
  --p-max-ee-r 4 \
  --o-table table-dada2.qza \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats stats-dada2.qza

## 3.1 - Visualización de las secuencias tras el proceso de denoising con dada2
# Se van a visualizar las estadísticas tras el filtrado realizado (stats-dada2.qzv) y las secuencias representativas 
# rep-seqs-dada2.qzv (FeatureData)

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs-dada2.qzv 

# Exportación de tabla obtenida y de las secuencias
# Se pasa de formato .biom a .tsv para su visualización

qiime tools export \
  --input-path table-dada2.qza \
  --output-path asv_table

biom convert -i asv_table/feature-table.biom -o asv_table/asv-table.tsv --to-tsv 

qiime tools export \
  --input-path rep-seqs-dada2.qza \
  --output-path asvs

# Para visualización de números de secuenciación por muestra (FeatureTable)

qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-metadata-dada2.qzv \
  --m-sample-metadata-file metadata.tsv

# Con la visualización de table-metadata-dada2.qzv se puede obtener una estimación de la profundidad de muestreo uniforme óptima (rarefacción) para el análisis.
# La mayoría de métricas calculadas en el posterior análisis de diversidad son sensibles a este parámetro, por lo que su buena estimación es primordial
# para obtener unos buenos resultados.
# Para seleccionar el mejor valor de este parámetro hay que escoger el valor que sea lo más alto posible, y con ello retener más secuencias por muestra, y que
# a su vez excluya el menor número de muestras, ya a que a mayor profundidad, más probabilidades de que el recuento total de algunas muestras sea menor al valor
# escogido y que por tanto, se eliminen al seleccionar dicho valor de profundidad.

# Tras analizar el archivo table-metadata-dada2.qzv se obtiene que la mejor profundidad es de 32607 = Retención de 684.747 (39,52%) rasgos en 21 (100,00%) muestras 
# a la profundidad de muestreo especificada.
# Como se ha mencionado, se utilizará dicha profundidad posteriormente para el cálculo de las métricas de diversidad alfa y beta en el análisis taxonómico.

## 4 - ANÁLISIS TAXONÓMICO
# A continuación se explorará la composición taxonómica de las muestras y dichos resultados se relacionarán con los metadatos de nuestras muestras.
# El objetivo de dicho análisis taxonómico es determinar si la composición de las poblaciones bacterianas en la moco de MET2 y MET3 varía
# con respecto al grupo control o MET1.

# Precisamos usar una base de datos con la que comparar la composición bacteriana de nuestras muestras. En nuestro caso se usó la de Silva 138-99, que es
# la más actualizada hasta el momento en cuanto a bases de datos 16S https://www.arb-silva.de/documentation/release-138/.

# A esta base de datos se le asocian los primers usados para secuenciar nuestras muestras con el objetivo de reducir la base de datos al amplicón de interés 
# que queremos analizar. Tras eso se obtendrá el clasificador Naive Bayes previamente entrenado.

# Se crea un directorio en el que alojar la base de datos y el clasificador:

mkdir silva
cd silva

# Descarga de los archivos .qza para qiime2 de la base de datos Silva:

wget https://data.qiime2.org/2022.4/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2022.4/common/silva-138-99-tax.qza

# Extracción de lecturas de referencia (amplicón)

qiime feature-classifier extract-reads \
  --i-sequences $HOME/analisis_bioinformatico/silva/silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCNGCAG \
  --p-r-primer GACTACNNGGGTATCTAATCC \
  --o-reads $HOME/analisis_bioinformatico/silva/silva-138-99-extracts.qza

# Entrenamiento de un clasificador Naive Bayes

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $HOME/analisis_bioinformatico/silva/silva-138-99-extracts.qza \
  --i-reference-taxonomy $HOME/analisis_bioinformatico/silva/silva-138-99-tax.qza \
  --o-classifier $HOME/analisis_bioinformatico/silva/silva-138-99-classifier.qza

# Clasificación de las secuencias representativas con feature-classifier classify-sklearn 
# Permite la clasificación taxonómica de las secuencias empleando el clasificador Silva entrenado previamente:

qiime feature-classifier classify-sklearn \
  --i-classifier $HOME/analisis_bioinformatico/silva/silva-138-99-classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification silva_tax_sklearn.qza

## 4.1 - Visualización de composición taxonómica de cada muestra
# Se empleará para ello taxa barplot para crear el histograma resultante de la asignación taxonómica de las muestras usando table-dada2.qza, el archivo de 
# clasificación taxonómica silva_tax_sklearn.qza y el archivo de metadatos con la información de las muestras dependiendo de la dieta experimental suministrada.
# Se obtendrá un archivo .qzv barplot_taxonomoy.qzv en el que se puede visualizar dicha asignación por niveles taxonómicos (del 1 al 7).

qiime taxa barplot \
  --i-table table-dada2.qza \
  --i-taxonomy silva_tax_sklearn.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization barplot_taxonomoy.qzv

## 4.2 Exportación de tablas de asignacion taxonomica
# Se produce un archivo tabulated-combined-metadata.qzv con la asignación de cada secuencia a una id de referencia según el análisis taxonómico.

qiime tools export \
  --input-path silva_tax_sklearn.qza \
  --output-path asv_tax_dir

mv asv_tax_dir/taxonomy.tsv asv_tax_dir/silva_taxonomy.tsv

qiime metadata tabulate \
  --m-input-file rep-seqs-dada2.qza \
  --m-input-file silva_tax_sklearn.qza \
  --o-visualization tabulated-combined-metadata.qzv

## 5 - ALINEAMIENTO Y FILOGENIA
# Se generará un árbol filogenético que relaciona características entre sí.
# Para generar un árbol filogenético usamos el comando align-to-tree-mafft-fasttree. 
# El programa mafft produce una alineación de secuencia múltiple en nuestro FeatureData[Sequence] (rep-seqs-dada2.qza) para crear 
# un FeatureData[AlignedSequence] (aligned-rep-seqs.qza). 
# Después, el programa masks filtra el alineamiento para eliminar las posiciones que son muy variables (masked-aligned-rep-seqs.qza), y el programa FastTree
# generar un árbol filogenético a partir del alineamiento enmascarado. Este programa FastTree crea un árbol sin raíces (unrooted-tree.qza), 
# por lo que el último paso de este comando consta de aplicar un enraizamiento de punto medio (rooted-tree) para colocar la raíz del árbol en el punto medio de la 
# distancia más larga entre las puntas del árbol filogenético (rooted-tree.qza).

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

## 6 - ANÁLISIS DE DIVERSIDAD
## 6.1 - Rarefacción alfa
# Se usa para explorar y visualizar la variación de los índices alpha de acuerdo a la profundidad de secuenciación. 
# Esta visualización calcula una o más métricas de diversidad alpha en múltiples profundidades de muestreo, en pasos entre 1 
# (controlado con --p-min-depth si se aplica) y el valor proporcionado como --p-max-depth.
# En este caso aplicamos una profundidad máxima de 32607, que es la que hemos estimado anteriormente como la mejor para este análisis.

qiime diversity alpha-rarefaction \
  --i-table table-dada2.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 32607 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

## 6.2 - Índices de diversidad alfa y beta
# El análisis de diversidad admite el cálculo de métricas de diversidad alpha y beta, la aplicación de pruebas estadísticas relacionadas 
# y la generación de visualizaciones de los resultados de los índices de diversidad. 

# Primero aplicamos el comando core-metrics-phylogenetic a una profundidad especificada usando como input el árbol filogenético rooted-tree.qza 
# y la table-dada2.qza (FeatureTable[Frequency]), tabla que aporta información sobre la longitud de las de secuencias, estadísticas sobre la longitud 
# de las secuencias y de las secuencias de las que disponemmos. Con esto se generarán todos los índices alpha y beta a la vez y se alojarán en el 
# directorio creado (--output-dir metrics_32607).
# La profundidad escogida para calcular las métricas de diversidad alpha y beta la calculamos anteriormente mediante prueba y error con la visualización 
# del archivo table-metadata-dada2.qzv en el apartado 3.1 del análisis.

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-dada2.qza \
  --p-sampling-depth 32607 \
  --m-metadata-file metadata.tsv \
  --output-dir metrics_32607

# Tras calcular los archivos .qza de todos los índices de diversidad alpha y beta se procede a preparar sus archivos de visualización .qzv:

## Índices de diversidad alpha
# Los índices de diversidad alpha son:
# - Índice de diversidad de Shannon: una medida cuantitativa de la riqueza de la comunidad.
# - Características observadas: una medida cualitativa de la riqueza de la comunidad.
# - Diversidad filogenética de Faith: una medida cualitativa de la riqueza de la comunidad que incorpora las relaciones filogenéticas entre los rasgos.
# - Evennes (igualdad o igualdad de Pielou): una medida de la igualdad de la comunidad.

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_32607/faith_pd_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization metrics_32607/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_32607/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization metrics_32607/shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_32607/evenness_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization metrics_32607/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity metrics_32607/observed_features_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization metrics_32607/observed-features-group-significance.qzv

## Índices de diversidad beta
# Los índices de diversidad beta son:
# - Distancia de Jaccard: una medida cualitativa de la disimilitud de las comunidades.
# - Distancia Bray-Curtis: una medida cuantitativa de la disimilitud de las comunidades.
# - Distancia UniFrac no ponderada: una medida cualitativa de la disimilitud de la comunidad que incorpora las relaciones filogenéticas entre las características.
# - Distancia UniFrac ponderada: una medida cuantitativa de la disimilitud de la comunidad que incorpora las relaciones filogenéticas entre los rasgos.

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_32607/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --p-pairwise \
  --o-visualization metrics_32607/unweighted-unifrac-treatment1-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_32607/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --p-pairwise \
  --o-visualization metrics_32607/weighted-unifrac-treatment1-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_32607/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --p-pairwise \
  --o-visualization metrics_32607/bray_curtis-treatment1-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix metrics_32607/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --p-pairwise \
  --o-visualization metrics_32607/jaccard-treatment1-significance.qzv

## 7 - ANÁLISIS DE ABUNDANCIA DIFERENCIAL
# El análisis de abundancia diferencial se usa para identificar características que son diferencialmente abundantes 
# (es decir, presentes en diferentes abundancias) en los grupos de muestra. En nuestro caso nos sirve para determinar si hay microorganismos más abundantes en alguna
# en las muestras y si tiene relación con la concentración de metionina suministrada. Este análisis diferencial se puede realizar usando q2-composition o q2-gneiss.

# Se usó ANCOM para hacer el análisis de abundancia diferencial entre muestras. Primero se hizo dicha prueba de abundancia diferencial 
# a un nivel taxonómico específico usando la tabla de datos table-dada2.qza con taxa collapse. Se añade un pseudoconteo para evitar errores
# en el análisis estadístico y finalmente se aplica el plugin q2-composition ancom para realizar dicho análisis de abundancia, seleccionando
# el archivo de metadatos e indicando la columna sobre la que se quiere analizar. En nuestro caso interesa analizar la columna de los tratamientos o
# dietas experimentales. Este análisis de diversidad se realizó a los siete niveles taxonómicos que admitía el programa, obteniendo un archivo .qzv que permite
# visualizar las características más abundantes por grupo o dieta experimental.

## 7.1 - Abundancia diferencial con q2-composition (ANCOM)

# Nivel 2 taxonomía (Phylum)

qiime taxa collapse \
  --i-table table-dada2.qza \
  --i-taxonomy silva_tax_sklearn.qza \
  --p-level 2 \
  --o-collapsed-table table-level2.qza

qiime composition add-pseudocount \
  --i-table table-level2.qza \
  --o-composition-table comp-table-level2.qza

qiime composition ancom \
  --i-table comp-table-level2.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization ancon-level2.qzv \

# Nivel 3 taxonomía (Clase)

qiime taxa collapse \
  --i-table table-dada2.qza \
  --i-taxonomy silva_tax_sklearn.qza \
  --p-level 3 \
  --o-collapsed-table table-level3.qza

qiime composition add-pseudocount \
  --i-table table-level3.qza \
  --o-composition-table comp-table-level3.qza

qiime composition ancom \
  --i-table comp-table-level3.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization ancon-level3.qzv \

# Nivel 4 taxonomía (Orden)

qiime taxa collapse \
  --i-table table-dada2.qza \
  --i-taxonomy silva_tax_sklearn.qza \
  --p-level 4 \
  --o-collapsed-table table-level4.qza

qiime composition add-pseudocount \
  --i-table table-level4.qza \
  --o-composition-table comp-table-level4.qza

qiime composition ancom \
  --i-table comp-table-level4.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization ancon-level4.qzv \

# Nivel 5 taxonomía (Familia)

qiime taxa collapse \
  --i-table table-dada2.qza \
  --i-taxonomy silva_tax_sklearn.qza \
  --p-level 5 \
  --o-collapsed-table table-level5.qza

qiime composition add-pseudocount \
  --i-table table-level5.qza \
  --o-composition-table comp-table-level5.qza

qiime composition ancom \
  --i-table comp-table-level5.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization ancon-level5.qzv \

# Nivel 6 taxonomía (Género)

qiime taxa collapse \
  --i-table table-dada2.qza \
  --i-taxonomy silva_tax_sklearn.qza \
  --p-level 6 \
  --o-collapsed-table table-level6.qza

qiime composition add-pseudocount \
  --i-table table-level6.qza \
  --o-composition-table comp-table-level6.qza

qiime composition ancom \
  --i-table comp-table-level6.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization ancon-level6.qzv \

# Nivel 7 taxonomía (Especie)

qiime taxa collapse \
  --i-table table-dada2.qza \
  --i-taxonomy silva_tax_sklearn.qza \
  --p-level 7 \
  --o-collapsed-table table-level7.qza

qiime composition add-pseudocount \
  --i-table table-level7.qza \
  --o-composition-table comp-table-level7.qza

qiime composition ancom \
  --i-table comp-table-level7.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization ancon-level7.qzv \

## 7.2 - Abundancia diferencial con q2-gneiss
# Tras aplicar ANCOM, se usó q2-gneiss para generar mapas de calor que relacionaran los tratamientos o dietas experimentales y que arrojaran conclusiones sobre 
# sus características más abundantes.

# Primero se usó correlation-clustering para definir particiones de microorganismos que comúnmente coexisten entre sí utilizando 
# el clustering jerárquico de Ward: si dos microbios están muy correlacionados, el valor de la métrica se reduce a cero. 
# El cluster jerárquico de Ward utiliza esta métrica de distancia para agrupar los grupos de microorganismos que están correlacionados entre sí. 
# Al final, el árbol que se obtiene resaltará la estructura de alto nivel e identificará cualquier bloque dentro de los datos (hierarchy.qza).

# Tras ello con ayuda de feature-table filter-samples y empleando la table-dada2.qza y el archivo de metadatos crearemos tres archivos, en los que filtraremos
# los tratamientos o dietas experimentales, obteniendo los archivos treatment1-table.qza, treatment2-table.qza y treatment3-table.qza, y usando estos datos, con
# ayuda del comando gneiss dendrogram-heatmap y usando estas tablas filtradas con tratamientos y el archivo hierarchy.qza se crearán los mapas de calor:

qiime gneiss correlation-clustering \
  --i-table table-dada2.qza \
  --p-pseudocount 1 \
  --o-clustering hierarchy.qza

qiime feature-table filter-samples \
  --i-table table-dada2.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[treatment]='MET1'" \
  --o-filtered-table treatment1-table.qza

qiime feature-table filter-samples \
  --i-table table-dada2.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[treatment]='MET2'" \
  --o-filtered-table treatment2-table.qza

qiime feature-table filter-samples \
  --i-table table-dada2.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[treatment]='MET3'" \
  --o-filtered-table treatment3-table.qza

qiime gneiss dendrogram-heatmap \
  --i-table treatment1-table.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --p-color-map Blues \
  --o-visualization tree_heatmap_tr1.qzv

qiime gneiss dendrogram-heatmap \
  --i-table treatment2-table.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --p-color-map Blues \
  --o-visualization tree_heatmap_tr2.qzv

qiime gneiss dendrogram-heatmap \
  --i-table treatment3-table.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --p-color-map Blues \
  --o-visualization tree_heatmap_tr3.qzv

## 8 -  ANÁLISIS FUNCIONAL Y EXPRESIÓN DIFERENCIAL
## 8.1 - Asignación funcional
# Se usa el pipeline picrust2 (picrust2 full-pipeline) para hacer asignación funcional en base a la composición de baterias indentificadas 
# en cada muestra así como para para predecir el contenido metagenómico. Se emplea para ello la tabla table-dada2.qza  y las secuencias
# representativas rep-seqs-dada2.qza. Se obtienen como resultados el archivo ko_results.qza, que es el metagenoma previsto para los ortólogos KEGG, 
# el archivo ec_results.qza, que corresponde al metagenoma previsto para números CE y el archivo metacyc_results.qza, que corresponde a la 
# predicción de las abundancias de las vías MetaCyc.

qiime picrust2 full-pipeline \
  --i-table table-dada2.qza \
  --i-seq rep-seqs-dada2.qza \
  --o-ko-metagenome ko_results.qza \
  --o-ec-metagenome ec_results.qza \
  --o-pathway-abundance metacyc_results.qza

qiime metadata tabulate \
  --m-input-file ec_results.qza \
  --o-visualization ec_results.qzv

qiime metadata tabulate \
  --m-input-file ko_results.qza \
  --o-visualization ko_results.qzv

qiime metadata tabulate \
  --m-input-file metacyc_results.qza \
  --o-visualization metacyc_results.qzv

## 8.2 - Análisis de expresión diferencial
# A continuación se va usar de nuevo q2-composition ancom para realizar un estudio de expresión diferencial que pueda determinar
# la presencia o no de alguna ruta metabólica enriquecida por la presencia de la concentración de metionina en los tratamientos.
# Por otro lado, el análisis de expresión diferencial se usa para determinar la presencia de rutas metabólicas enriquecidas en muestras.
# Para ello utilizamos las tables de los resultados de las asignaciones funcionales para implementarlas en q2-composition ancom, obteniendo los archivos
# de visualización del análisis diferencial ancom_ec-treatment.qzv, ancom_KO-treatment.qzv y ancom_metacyc-treatment.qzv.

qiime composition add-pseudocount \
  --i-table ec_results.qza \
  --o-composition-table ec_comp-table.qza

qiime composition ancom \
  --i-table ec_comp-table.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization ancom_ec-treatment.qzv

qiime composition add-pseudocount \
  --i-table ko_results.qza \
  --o-composition-table ko_comp-table.qza

qiime composition ancom \
  --i-table ko_comp-table.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization ancom_KO-treatment.qzv 

qiime composition add-pseudocount \
  --i-table metacyc_results.qza \
  --o-composition-table metacyc_comp-table.qza

qiime composition ancom \
  --i-table metacyc_comp-table.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization ancom_metacyc-treatment.qzv 











