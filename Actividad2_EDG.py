#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 14:40:14 2024

@author: emmettdiez
"""

#%%

# PRIMERA PARTE: OBTENER INFORMACIÓN DE PROTEÍNAS Y COMPUESTOS QUÍMICOS.
# Importo las librerías necesarias para esta práctica.
import requests
from time import sleep
import pandas as pd
import os
import pubchempy as pcp
from Bio import PDB
from Bio.PDB import  MMCIF2Dict
import libchebipy
from rdkit import Chem
from rdkit.Chem import AllChem

# Defino la ruta del Escritorio, donde almaceno los archivos.
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")

# Creo la lista con dichoso ID.
pdb_ids = ['1tup', '2xyz', '3def', '4ogq', '5jkl', '6mno', '7pqr', '8stu']

# Compruebo que la ruta es correcta y existe.
if not os.path.exists(desktop_path):
    os.makedirs(desktop_path)
    print(f"Directorio creado: {desktop_path}")
else:
    print("El directorio existe.")
#%%

# FUNCIONES NECESARIAS
# Creo una función para recibir los ID de PdBk y otra para descargarlos como .CIF

def download_pdb_file(pdb_id, output_path):
    "Descargo el archivo PDB usando libreria requests"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_path, 'wb') as file:
            file.write(response.content) # Uso content porque abro el archivo con "wb" que es formato binario. Sino, no puedo leer un string.
            print(f"Archivo PDBx/mmCIF descargado y guardado en '{output_path}'.")
    else:
            raise ValueError(f"Error al descargar el archivo PDBx/mmCIF: {pdb_id}")

def download_pdbx_mmCIF_file(pdb_id, output_path):
    """Descarga el archivo CIF usando requests."""
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            file.write(response.content)
        print(f"Archivo PDBx/mmCIF descargado y guardado en '{output_path}'.")
    else:
        raise ValueError(f"Error al descargar el archivo PDBx/mmCIF: {pdb_id}")


# Con la siguiente función, averiguo el ID de 1tup.pdb desde Uniprot.
# Esto se consigue gracias al URL de la API de UniProt para mapear IDs.

def obtener_uniprot_id_desde_pdb(pdb_id):
    url = "https://rest.uniprot.org/idmapping/run"
    # Parámetros para el mapeo de PDB a UniProt
    datos = {
        'from': 'PDB',  # Cambiado a 'PDB' para el mapeo correcto
        'to': 'UniProtKB',
        'ids': pdb_id
    }
    
    # Hacemos la solicitud POST para obtener el job_id
    respuesta = requests.post(url, data=datos)
    if respuesta.status_code == 200:
        # Obtenemos el job_id para hacer seguimiento del estado
        job_id = respuesta.json()['jobId']

        # Obtenemos el resultado del mapeo
        url_datos = f"https://rest.uniprot.org/idmapping/results/{job_id}"
        resultado_final = requests.get(url_datos).json()
        
        if 'results' in resultado_final and resultado_final['results']:
            # Retornamos el primer UniProt ID mapeado
            uniprot_id = resultado_final['results'][0]['to']
            
            return uniprot_id
        else:
            print(f"No se encontraron resultados para el PDB ID {pdb_id}.")
            return None
    else:
        print(f"No se pudo obtener el UniProt ID para el PDB ID {pdb_id}.")
        return None

# Creamos una funcion para obtener la información de UniProt:

def obtener_informacion_uniprot(uniprotid):
    if not uniprot_id:  # Verifica que el ID no sea None o vacío
        return None
        # URL de la API de UniProt
    url = f"https://www.uniprot.org/uniprot/{uniprotid}.json"
        
        # Hacemos la solicitud a la API
    respuesta = requests.get(url)
        
    if respuesta.status_code == 200:
            datos = respuesta.json()

            # Extraer información relevante
            fecha_publicacion = datos['entryAudit']['firstPublicDate']
            fecha_modificacion = datos['entryAudit']['lastAnnotationUpdateDate']
            revisado = "Swiss-Prot" if datos['entryType'] else "Trembl"
            nombre_gen = datos['genes'][0]['geneName']['value']
            sinonimos = ', '.join([syn['value'] for syn in datos['genes'][0].get('synonyms', [])])
            organismo = datos['organism']['scientificName']
            nombre_proteina = datos['proteinDescription']['recommendedName']['fullName']['value']
            secuencia_aa = datos['sequence']['value']
            pdb_ids = ', '.join([ref['id'] for ref in datos['uniProtKBCrossReferences'] if ref['database'] == 'PDB'])
            
            # Mostrar la información extraída
            print(f"UniProt ID: {uniprot_id}")
            print(f"Fecha de publicación: {fecha_publicacion}")
            print(f"Fecha de modificación: {fecha_modificacion}")
            print(f"Revisado: {revisado}")
            print(f"Nombre del gen: {nombre_gen}")
            print(f"Sinónimos del gen: {sinonimos}")
            print(f"Organismo: {organismo}")
            print(f"Nombre completo de la proteína: {nombre_proteina}")
            print(f"Secuencia de aminoácidos: {secuencia_aa}")
            print(f"PDB IDs asociados: {pdb_ids}")
            
            return pd.Series([fecha_publicacion, fecha_modificacion, revisado, nombre_gen, sinonimos, organismo, nombre_proteina, secuencia_aa, pdb_ids])
            
    else:
            print(f"No se pudo obtener la información de UniProt para el ID {uniprot_id}.")


#%%
# PROGRAMA PRINCIPAL

#PASO 1: Descargamos en Desktop los archivos PDB y CIF con las funciones anteriormente creadas.
for pdb_id in pdb_ids:
    
    # Usando desktop_path para asegurarte de que estás en el Escritorio
    pdb_file_path = os.path.join(desktop_path, f"{pdb_id}.pdb")
    pdbx_mmCIF_file_path = os.path.join(desktop_path, f"{pdb_id}.cif")

    # Descargo los archivos
    download_pdb_file(pdb_id, pdb_file_path)
    download_pdbx_mmCIF_file(pdb_id, pdbx_mmCIF_file_path)

print('Ha terminado la lectura de datos')    

# PASO 2: Creamos un dataframe con los ID
df = pd.DataFrame({'Id_pdb': pdb_ids})

uniprot_id = []

# Obtenemos los IDs de Uniprot

for id_ in pdb_ids:
    print('Estoy procesando el ID:', id_)
    uniprot_id.append(obtener_uniprot_id_desde_pdb(id_))

data = {'Uniprot_id':[], 'fecha_publicacion':[], 'fecha_modificacion':[], 'revisado':[], 'nombre_gen':[], 'sinonimos':[], 'organismo':[], 'nombre_proteina':[], 'secuencia_aa':[], 'pdb_id':[]}# Modificamos el dataframe con las siguientes columnas:
for x in uniprot_id:
    salida = obtener_informacion_uniprot(x)
    data['Uniprot_id'].append(x)
    data['fecha_publicacion'].append(salida[0])
    data['fecha_modificacion'].append(salida[1])

pass
    
#%%

# APARTADO C: investigar si existe cofactor del 1tup y ver dónde se identifica.
# Entramos a la API de PubChem para obtener el ID como CID
# URL de la API de PubChem para obtener información de un compuesto por CID

# Definir el identificador PDB de la proteína de la que queremos obtener el CID.
# Consulto la web de Protein Data Bank para obtener el cofactor ZN, y lo corroboro en Uniprot.Puede verse en el documento adjunto.
cofactor_id = 'ZN'

# Construir la URL para la API de PubChem
url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cofactor_id}/cids/JSON'

# Hacer la solicitud GET a la API de PubChem.
response = requests.get(url)

# Verificar si la respuesta fue exitosa
if response.status_code == 200:
    # Obtener el JSON de la respuesta
    data = response.json()
    
    # Imprimir el CID(s)
    if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
        print(f'El CID para {cofactor_id} es: {data["IdentifierList"]["CID"][0]}')
        cid = data["IdentifierList"]["CID"][0]
    else:
        print(f'No se encontró un CID para {cofactor_id}')
else:
    print(f"Error al obtener el CID: {response.status_code}")

url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON'

response = requests.get(url)
if response.status_code == 200:
    compound_info = response.json()
    print(compound_info)
else:
    print(f"Error: {response.status_code}")


# Convertir el diccionario en df

cids = compound_info.get('PC_Compounds', [{}])[0].get('id', {}).get('id', [])
print(cids['cid'])  # Suponiendo que el primer CID es el relevante
propiedades = compound_info.get('PC_Compounds', [{}])[0].get('props', {})

info_cofactor = {} # Creo un diccionario vacío para que me guarde la info de mi cofactor que le pido debajo.

for i in propiedades:
        
    var = i['urn']['label']
    value = i['value'].get('sval')
        
    if var == 'InChIKey':
        print(f'InCHIKey: {value}')
        info_cofactor ['InChIKey'] = [f'{value}']
    if var == 'InChI':
        print(f'InChI: {value}')
        info_cofactor ['InChI'] = [f'{value}']
    if var == 'Molecular Weight':
        print(f'Molecular Weight: {value}')
        info_cofactor ['Molecular Weight'] = [float(f'{value}')]
    if var == 'IUPAC Name':
        print(f'IUPAC Name: {value}')
        info_cofactor ['IUPAC Name'] = [f'{value}']

df_cofactor = pd.DataFrame(info_cofactor, index = None)

print(df_cofactor)

#%%
# PARTE 2: MANIPULACIÓN DE DATOS BIOLÓGICOS

# Parseamos el CIF 4ogq:
parser = PDB.MMCIFParser(QUIET=True)
pdbx_mmCIF_file_path = '/Users/emmettdiez/Desktop/4ogq.cif'
structure = parser.get_structure("protein", pdbx_mmCIF_file_path)
heteromoleculas = []

# Iteramos a través de los modelos, cadenas y residuos:
for model in structure:
    for chain in model:
        for residue in chain:
            # Excluir el agua ('HOH')
            print(residue.id[0], residue.resname)
            if residue.id[0] != ' ' and residue.resname != 'HOH':
                heteromoleculas.append(residue.resname)

print("Heteromoléculas encontradas (sin incluir agua):")
print(heteromoleculas)

mmcif_dict = MMCIF2Dict.MMCIF2Dict(pdbx_mmCIF_file_path) #Convierte 4ogq en un diccionario


# Extraigo los nombres y primeras letras de las heteromoléculas:
if mmcif_dict:
    nonpoly_names = mmcif_dict.get('_pdbx_entity_nonpoly.name', [])
    nonpoly_ids = mmcif_dict.get('_pdbx_entity_nonpoly.comp_id', [])
else:
    nonpoly_names = []
    nonpoly_ids = []

# Crear un DataFrame con los nombres y sus identificadores de tres letras:
df_heteromoleculas = pd.DataFrame({
    'Nombre': nonpoly_names,
    'ID de tres letras': nonpoly_ids
})

# Mostrar el DataFrame resultante
print("\nInformación de heteromoléculas:")
print(df_heteromoleculas)

# Guardo el archivo en formato CSV
df_heteromoleculas.to_csv('heteromoleculas.csv', index=False)

#%%

# Función para obtener el SMILES de una molécula a partir de su nombre en PubChem
def get_smiles_from_pubchem(Nombre):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{Nombre}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            return response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
        except (KeyError, IndexError):
            return None
    return None


df_heteromoleculas['SMILES'] = df_heteromoleculas['Nombre'].apply(get_smiles_from_pubchem)

print(df_heteromoleculas)

#%%

# Crear archivo SDF y calcular el peso molecular
file_path = "/Users/emmettdiez/Desktop/heteromolecules.sdf"
with Chem.SDWriter(file_path) as sdf_writer:
# with Chem.SDWriter("heteromolecules.sdf") as sdf_writer:
    for index, row in df_heteromoleculas.iterrows():
        try:
            mol = Chem.MolFromSmiles(row['SMILES'])
            if mol:
                mol.SetProp("Molecular_weight", str(AllChem.CalcExactMolWt(mol)))
                sdf_writer.write(mol)
        except:
            pass
            #print(row['SMILES'])
            
#%%    

# Verificación: cargar las moléculas del archivo SDF
supplier = Chem.SDMolSupplier(r"/Users/emmettdiez/Desktop/heteromolecules.sdf")
for mol in supplier:
    if mol is not None:
        print(Chem.MolToMolBlock(mol))
        print('--------------------------')
        print(mol.GetProp("Molecular_weight"))

