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
