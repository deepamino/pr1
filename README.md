<h1 align="center">Apino</h1>

<p align="center">La API de DeepAmino</p>

<br>
<p align="justify">
APINO es un proyecto open-source diseñado para ofrecer una API intuitiva y eficiente que permite la extracción de información de secuencias de nucleótidos desde el NCBI (National Center for Biotechnology Information). Nuestro objetivo es facilitar el acceso a información genómica relevante mediante una interfaz simplificada y accesible, optimizando así el trabajo de investigadores, bioinformáticos y entusiastas de la biotecnología.
</p>

## Objetivo del Proyecto

<p align="justify">
El propósito de APINO es brindar una solución sencilla y de fácil manejo para la consulta de datos genómicos específicos mediante identificadores universales. Utilizando los identificadores únicos asignados por el NCBI para cada secuencia, los usuarios pueden extraer información esencial sobre secuencias de nucleótidos de manera rápida y eficiente. Nos enfocamos principalmente en secuencias relacionadas con genes, lo cual facilita la integración de APINO en estudios de genética, proyectos de investigación y aplicaciones clínicas.
</p>

## Características Principales

- **Extracción de datos mediante API** <p align="justify">Con APINO, cualquier usuario puede acceder a la información genética relevante simplemente proporcionando un identificador de secuencia.</p>

- **Compatibilidad con múltiples tipos de secuencias** <p align="justify">Aunque el enfoque principal es en genes, APINO permite la consulta de diversos tipos de secuencias nucleotídicas presentes en el NCBI.</p>

- **Fácil integración y escalabilidad** <p align="justify">Diseñado para ser fácilmente integrable en otras aplicaciones y adaptable a las necesidades del usuario.</p>

## Ejemplo de uso

La API ofrece dos métodos principales en las rutas `/info` y `/sequence`, los cuales son ambos métodos GET que reciben como `query param` el identificador de la secuencia. A continuación, se muestra un ejemplo de uso de ambas:


- Para obtener información de una determinada secuencia:
  
```bash
GET http://<domain>:8080/info?id=NM_001301717
```

```json
{
   "id":"NM_001301717",
   "locus":"NM_001301717 2191 bp mRNA linear PRI 01-OCT-2024",
   "definition":"Homo sapiens C-C motif chemokine receptor 7 (CCR7), transcript",
   "source":"Homo sapiens (human)",
   "organism":"Homo sapiens Eukaryota; Metazoa; Chordata",
   "references":[
      {
         "authors":"Gao,J., Wang,Z., Lin,S., Tian,Y., Wu,H., Li,Z. and Liu,F.",
         "title":"CCR7/DUSP1 signaling Axis mediates iCAF to regulates head and neck squamous cell carcinoma growth ",
         "Journal":"Cell Signal 122, 111305 (2024)"
      },
      {
         "authors":"Cardona,C.I., Rodriguez,A., Torres,V.C., Sanchez,A., Torres,A.,",
         "title":"C-C Chemokine Receptor 7 Promotes T-Cell Acute Lymphoblastic Leukemia Invasion of the Central Nervous System via beta2-Integrins ",
         "Journal":"Int J Mol Sci 25 (17), 9649 (2024)"
      },
      {
         "authors":"Geraldo,L.H., Garcia,C., Xu,Y., Leser,F.S., Grimaldi,I., de Camargo",
         "title":"CCL21-CCR7 signaling promotes microglia/macrophage recruitment and chemotherapy resistance in glioblastoma ",
         "Journal":"Cell Mol Life Sci 80 (7), 179 (2023)"
      },
      {
         "authors":"Xu,D., Liu,X., Ke,S., Guo,Y., Zhu,C. and Cao,H.",
         "title":"CCL19/CCR7 drives regulatory T cell migration and indicates poor prognosis in gastric cancer ",
         "Journal":"BMC Cancer 23 (1), 464 (2023)"
      },
      {
         "authors":"Moschovakis,G.L. and Forster,R.",
         "title":"Multifaceted activities of CCR7 regulate T-cell homeostasis in health and disease ",
         "Journal":"Eur J Immunol 42 (8), 1949-1955 (2012)"
      },
      {
         "authors":"Yoshida,R., Imai,T., Hieshima,K., Kusuda,J., Baba,M., Kitaura,M.,",
         "title":"Molecular cloning of a novel human CC chemokine EBI1-ligand chemokine that is a specific functional ligand for EBI1, CCR7 ",
         "Journal":"J Biol Chem 272 (21), 13803-13809 (1997)"
      },
      {
         "authors":"Burgstahler,R., Kempkes,B., Steube,K. and Lipp,M.",
         "title":"Expression of the chemokine receptor BLR2/EBI1 is specifically transactivated by Epstein-Barr virus nuclear antigen 2 ",
         "Journal":"Biochem Biophys Res Commun 215 (2), 737-743 (1995)"
      },
      {
         "authors":"Schweickart,V.L., Raport,C.J., Godiska,R., Byers,M.G., Eddy,R.L.",
         "title":"Cloning of human and mouse EBI1, a lymphoid-specific G-protein-coupled receptor encoded on human chromosome 17q12-q21.2 ",
         "Journal":"Genomics 23 (3), 643-650 (1994)"
      },
      {
         "authors":"Birkenbach,M., Josefsen,K., Yalamanchili,R., Lenoir,G. and Kieff,E.",
         "title":"Epstein-Barr virus-induced genes: first lymphocyte-specific G protein-coupled peptide receptors ",
         "Journal":"J Virol 67 (4), 2209-2220 (1993)"
      },
      {
         "authors":"Itoh,H., Toyama,R., Kozasa,T., Tsukamoto,T., Matsuoka,M. and",
         "title":"Presence of three distinct molecular species of Gi protein alpha subunit. Structure of rat cDNAs and human genomic DNAs ",
         "Journal":"J Biol Chem 263 (14), 6656-6664 (1988)"
      }
   ]
}
```

- Para obtener la secuancia de un identificador:

```bash
GET http://<domain>:8080/sequence?id=NM_001301717
```

```json
{
      "sequences": 
             {
                    "NM_001301717":"CTCTAGATGAGTCAGTGGAGGG..."
             }
}
```



<br>

<p align="center">© DeepAmino 2024</p>
