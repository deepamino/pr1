
class ProteinInfo:

    description = ""
    locus = ""
    source = ""
    organism = ""
    comment = ""
    gene = ""
    gene_synonyms = ""
    articles = ""

    def set_description(self, description):
        self.description = description
    
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
    
    def set_gene_synonyms(self, gene_synonyms):
        self.gene_synonyms = gene_synonyms
    
    def set_articles(self, articles):
        self.articles = articles