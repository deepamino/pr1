from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from process.protein import Protein

class BiopythonTranslator:

    def read_fasta(self, fasta_file):
        sequences = []
        try:
            with open(fasta_file, "r") as fasta:
                for record in SeqIO.parse(fasta, "fasta"):
                    sequences.append(record)
        except FileNotFoundError:
            raise Exception(f"Archivo {fasta_file} no encontrado.")
        
        return sequences
    
    def translate_to_protein(self, nucleotide_seq):
        try:
            seq_obj = Seq(nucleotide_seq)
            protein_seq = seq_obj.translate()
        except Exception as e:
            raise Exception(f"Error en la traducción de nucleótidos a aminoácidos: {str(e)}")
        
        return protein_seq
    
    def write_protein_fasta(self, protein_sequences, output_file):
        try:
            with open(output_file, "w") as output_handle:
                SeqIO.write(protein_sequences, output_handle, "fasta")
        except Exception as e:
            raise Exception(f"Error al escribir el archivo FASTA: {str(e)}")
        
    def calculate_gc_content(self, nucleotide_seq):
        gc_pairs_count = nucleotide_seq.count('GC')
        gc_content = (gc_pairs_count * 2 / len(nucleotide_seq)) * 100
        return round(gc_content,2)
        
    def process_fasta(self, fasta_file, output_fasta):
        nucleotide_sequences = self.read_fasta(fasta_file)
        protein_records = [] 
        proteins_results = []

        for record in nucleotide_sequences:
            gc_content = self.calculate_gc_content(record.seq)

            protein_seq = self.translate_to_protein(record.seq)
            
            protein_result = Protein(id=record.id, description=record.description, protein_sequences=protein_seq, gc=gc_content)
            protein_record = SeqRecord(Seq(protein_seq), id=record.id, description=record.description)

            proteins_results.append(protein_result)
            protein_records.append(protein_record)

        self.write_protein_fasta(protein_records, output_fasta)

        return proteins_results


