import pandas as pd

"""
The TSV format presents the match data in columns as follows:

Protein Accession (e.g. P51587)
Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
Sequence Length (e.g. 3418)
Analysis (e.g. Pfam / PRINTS / Gene3D)
Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
Signature Description (e.g. BRCA2 repeat profile)
Start location
Stop location
Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
Status - is the status of the match (T: true)
Date - is the date of the run
(InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
(InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
(GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
(Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
"""

class Interpro:
    """ Interpro object (now only use for maker) """
    def __init__(self, strain):
        self.strain = strain # like UMSG1
        self.table = pd.read_table("data/" + strain + ".tsv", sep="\t", names = ["proAcc", "md5", "seqLen", "anaMethod", "sigAcc", "sigDes", "start", "end", "score", "status", "date", "intAnn", "intDes", "GO", "pathway"])
        self.ttable = self.table[self.table.score != "-"] # Phobius, MobiDBLite, Coils, ProSitePatterns do not have scores, so ommitted
        self.ttable = self.ttable[self.ttable.score.astype("float") < 0.000001] # only score < 1e-6
        self.ttable = self.ttable[self.ttable.GO.notnull()]  # GO is not NA
        
    def generate_go_report_for_wego(self):
        proAcc2GOs = dict()
        for index, row in self.ttable.iterrows():
            proAcc = row['proAcc']
            GOs = row['GO']
            for GO in GOs.split("|"):
                if not proAcc in proAcc2GOs.keys():
                    proAcc2GOs[proAcc] = set()
                proAcc2GOs[proAcc].add(GO)

        with open("output/" + self.strain + ".go.wego.txt", "w") as f:        
            for proAcc in sorted(proAcc2GOs.keys()):
                f.write(proAcc + "\t" + "\t".join(sorted(list(proAcc2GOs[proAcc]))) + "\n")
    
if __name__ == '__main__':
    for strain in ['Bgh2', 'BghR1', 'Bgt', 'Ene', 'UCSC1', 'UMSG1', 'UMSG2', 'UMSG3']:
        i = Interpro(strain)
        i.generate_go_report_for_wego() 


