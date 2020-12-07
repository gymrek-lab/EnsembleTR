"""
ReliabilityMetrics object includes computation of metrics indicating the reliability of a record
It includes locus-level and allele-level scores
"""

from trtools.utils.utils import GetHardyWeinbergBinomialTest

class ReliabilityMetrics:
    def __init__(self, record_object, faminfo=None, popinfo=None):
        self.record_object = record_object
        self.faminfo = faminfo
        self.popinfo = popinfo
        self.locus_stats = {}
        self.allele_stats = {}
        self.ComputeLocusStats()
        self.ComputeAlleleStats()

    def ComputeLocusStats(self):
        # Call rate
        self.locus_stats["callrate"] = self.record_object.record.call_rate

        # HWE by population
        samplelists = []
        if self.popinfo is None:
            samplelists.append(None)
        else:
            for pop in self.popinfo:
                samplelists.append(self.popinfo[pop])
        hwepvals = []
        for sl in samplelists:
            allele_freqs = self.record_object.hm_record.GetAlleleFreqs(samplelist=sl)
            genotype_counts = self.record_object.hm_record.GetGenotypeCounts(samplelist=sl)
            hwepvals.append(GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts))
        self.locus_stats["hwe-all"] = hwepvals # for each population
        self.locus_stats["hwe-worst"] = min(hwepvals) # worst of any population

        # Mendelian inheritance - TODO (if we have family info)

    def ComputeAlleleStats(self):
        pass # TODO

    def __str__(self):
        strout = "Locus stats for %s:\n"%self.record_object.vcf_type
        for stat in self.locus_stats:
            strout += "    %s=%s\n"%(stat, self.locus_stats[stat])
        return strout
        
