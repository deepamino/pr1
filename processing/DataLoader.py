from abc import ABC, abstractmethod
from Bio import Entrez, SeqIO


class DataLoader(ABC):
    @abstractmethod
    def collect(self, *args):
        pass


class ApiDataLoader(DataLoader):

    def __init__(self):
        self.email = 'example@gmail.com'
        self.type = 'nuccore'

    def collect(self, retmax=10):
        Entrez.email = self.email
        
        if not isinstance(retmax, int) or retmax <= 0:
            raise ValueError("retmax must be a positive integer.")

        chromosome_ids = self.search_chromosome("Homo sapiens[ORGN] AND chromosome", retmax=retmax)
        sequences = self.fetch_sequences(chromosome_ids)

        for seq_record in sequences:
            yield (seq_record.id, seq_record.seq)

    def fetch_sequences(self, id_list):
        ids = ",".join(id_list)
        handle = Entrez.efetch(db=self.type, id=ids, rettype="fasta", retmode="text")
        sequences = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        return sequences

    def search_chromosome(self, term, retmax=10):
        handle = Entrez.esearch(db=self.type, term=term, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]
        print(f"Found {len(id_list)} results.")
        return id_list
    

class FastaLoader(DataLoader):
    def collect(self, file):
        dna_seq = SeqIO.parse(file, 'fasta')
        for seq in dna_seq:
            yield seq


class LoaderFactory:
    _collectors = {
        'fasta': FastaLoader,
        'APIncbi': ApiDataLoader
    }

    @staticmethod
    def initialize_loader(key):
        loader_class = LoaderFactory._collectors.get(key)
        if loader_class:
            return loader_class()
        else:
            raise ValueError(f"Invalid collector key: {key}")