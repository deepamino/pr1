import requests
from .ncbi_model.bio_object_info import BioObjectInfo
from .ncbi_model.protein_ncbi import ProteinNCBI
from .ncbi_model.biopython_translator import BiopythonTranslator

class ApiManager:

    _url='http://localhost:8080/info'
    _url_sequence='http://localhost:8080/sequence'

    def request_info(self, id):
        params = {'id': id}
        response = requests.get(ApiManager._url, params=params)

        if response.status_code == 200:
            response_json = response.json()
            
            info = self._get_info(response_json)
            return info
        
        else:
            return None
    
    def get_sequences(self, id):
        params = {'id': id}
        response = requests.get(ApiManager._url_sequence, params=params)
        translator = BiopythonTranslator()

        if response.status_code == 200:
            response_json = response.json()
            sequence = response_json['sequences'][id]
            gc = translator.calculate_gc_content(sequence)
            aminoacids = translator.translate_to_protein(sequence)
            
            return ProteinNCBI( id=id, description="", protein_sequences=aminoacids, gc=gc)
        
        else:
            return None
    
    def _get_info(self, json):
        info = BioObjectInfo()

        for key in json:
            if hasattr(info, key):
                if key == 'geneSynonyms' or key == 'organism':
                    json[key] = json[key].replace(';', ',')
                setattr(info, key, json[key])

        return info
    
    
    
