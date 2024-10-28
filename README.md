# Práctica 1: Operatividad con Biopython

**Participantes:**
- Ricardo Juan Cárdenes Pérez
- Susana Suárez Mendoza

Esta práctica consiste en la elaboración de un programa en Python que tome como entrada un fichero FASTA, con varias secuencias de ADN, que haga las siguientes operaciones:

- Calcular el porcentaje de nucleótidos **GC** en cada secuencia y lo muestre por pantalla
- Producir las cadenas de **aminoácidos** asociadas a **cada secuencia** de ADN
- Escribir estas secuencias de aminoácidos en otro **fichero FASTA** con los mismo identificadores que las secuencias de ADN

## 1. Implementación

Todo lo mencionado previamente ha sido implementado, incluyendo una interfaz de usuario accesible a través de una página web desarrollada con la librería `Flask` de Python. A continuación, se describen los diferentes módulos y componentes que conforman el proyecto.
- `/uploads`: este directorio almacena temporalmente los archivos FASTA que el usuario carga en la página web para su traducción de nucleótidos a aminoácidos. Una vez completada la traducción, los archivos FASTA son eliminados. Además, contiene cuatro archivos que no se eliminan, los cuales corresponden a las cuatro proteínas más frecuentes en el apartado de "Most Searched".
- `/results`: este directorio almacena los archivos FASTA que han sido traducidos a secuencias de aminoácidos, permitiendo al usuario descargarlos para su uso. Una vez que el usuario realiza la descarga, los archivos son eliminados permanentemente de este directorio.
- `/process`: este módulo implementa la lógica necesaria para traducir secuencias de nucleótidos a aminoácidos, calcular el porcentaje de contenido de GC y gestionar la obtención de información de proteínas mediante una API desarrollada en un módulo complementario en Java.
  - `/translate_factory.py`:La clase `TranslatorFactory` es una fábrica que, a partir de una clave específica (`key`), devuelve una instancia del traductor correspondiente, como `BiopythonTranslator`, y lanza un error si la clave no es válida.
  - `/biopython_translator.py`: La clase `BiopythonTranslator` facilita la lectura de archivos FASTA, traduce secuencias de nucleótidos a proteínas, calcula el contenido de GC, escribe secuencias de proteínas en formato FASTA, y permite eliminar archivos de manera opcional, gestionando errores en cada etapa del proceso.
  - `/protein_ncbi.py`:La clase `ProteinNCBI` representa una proteína con sus atributos principales, que incluyen un identificador (`id`), una descripción, la secuencia de aminoácidos (`protein_sequences`), y el contenido de GC (`gc`).
  - `/bio_object_info.py`: La clase `BioObjectInfo` encapsula información biológica clave, como definición, locus, fuente, organismo, comentario, gen, sinónimos de gen y referencias, y proporciona métodos `set_` para establecer cada uno de estos atributos.
  - `/api_manager.py`: La clase `ApiManager` facilita la interacción con una API externa para obtener información biológica y secuencias de proteínas mediante peticiones HTTP, procesando la respuesta en objetos de las clases `BioObjectInfo` y `ProteinNCBI`, y traduce secuencias de nucleótidos a proteínas con cálculo de GC usando `BiopythonTranslator`.

<div align="center">
    <img src="images_readme/uml_process.png" alt="Diagrama Process" width=700 />
      <p><strong>Figura 1.</strong> Diagrama de clases para el módulo "process".</p> 
  </div>

