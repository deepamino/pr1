from .ncbi_model.biopython_translator import BiopythonTranslator

class TranslatorFactory:

    __translators = {
        'Biopython': BiopythonTranslator
    }

    @staticmethod
    def initialize_collector(key):
        translator_class = TranslatorFactory.__translators.get(key)
        if translator_class is None:
            raise ValueError('Invalid translator key')
        return translator_class()