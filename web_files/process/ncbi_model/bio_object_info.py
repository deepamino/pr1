
class BioObjectInfo:

    definition = ""
    locus = ""
    source = ""
    organism = ""
    comment = ""
    gene = ""
    geneSynonyms = ""
    references = ""

    def set_definition(self, definition):
        self.definition = definition
    
    def set_locus(self, locus):
        self.locus = locus
    
    def set_source(self, source):
        self.source = source
    
    def set_organism(self, organism):
        self.organism = organism

    def set_comment(self, comment):
        self.comment = comment
    
    def set_gene(self, gene):
        self.gene = gene
    
    def set_geneSynonyms(self, geneSynonyms):
        self.geneSynonyms = geneSynonyms
    
    def set_references(self, references):
        self.references = references