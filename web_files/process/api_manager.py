import requests
from .protein_info import ProteinInfo

class ApiManager:

    _url='http://localhost:8080/info'

    def reques_info(self, id):
        params = {'id': id}
        response = requests.get(ApiManager._url, params=params)

        if response.status_code == 200:
            response_json = response.json()
            
            info = self._get_info(response_json)
            return info
        
        else:
            return None
    
    def _get_info(self, json):
        info = ProteinInfo()

        for key in json:
            if hasattr(info, key):
                setattr(info, key, json[key])

        return info
    
    
