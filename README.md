# 🎭Mineria-de-datos-bioinformaticos
En este proyecto se pone en práctica la aplicación de las herramientas de minería de datos bioinformáticos:
- Manejo de API para realizar consultas a bases de datos biológicas y químicas.
- Uso de bibliotecas especializadas en la manipulación de estos datos.

Utilizaremos distintos archivos:
- heteromoleculas.csv, que contiene los datos de nombre e ID de cada molécula.
- Archivos .cif correspondientes a los ID de cada molécula extraídos de la API de Protein Data Bank.

👩🏽‍💻 ¿QUÉ CONTIENE ESTE CÓDIGO?

1️⃣ Obtención de información de PROTEÍNAS y COMPUESTOS BIOLÓGICOS:
  - **Funcion** de creación de la lista con ID's desde la **API** y **función** de descarga como .cif
2️⃣ PROGRAMA PRINCIPAL: extraigo el ID 1tup y determino su identificador en UniProt a través de su API. 
3️⃣ PANDAS DATAFRAME: creación del df con las columnas necesarias.
4️⃣ Investigar: ID Inchy, InchyKey, ID iupac...almacenarlo como compuesto.
5️⃣ Biopython y rdkit: aplicación de MMCIFParser. Obtención de sus SMILES.
