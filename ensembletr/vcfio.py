"""
This file contains classes for reading/writing 
VCF files.
"""

import trtools.utils.common as common
import trtools.utils.utils as utils
import trtools.utils.mergeutils as mergeutils
import trtools.utils.tr_harmonizer as trh
import cyvcf2
import vcf

from . import recordcluster as recordcluster

class VCFWrapper:
    """
    Simple class to keep track of VCF files and
    associated attributes

    Parameters
    ----------
    reader : cyvcf2.VCF
       VCF Reader
    vcftype : trh.TRRecordHarmonizer.vcftype
       Type of the VCF file (e.g. Hipstr, GangSTR, etc.)

    Attributes
    ----------
    reader : cyvcf2.VCF
       VCF Reader
    vcftype : trh.TRRecordHarmonizer.vcftype
       Type of the VCF file (e.g. Hipstr, GangSTR, etc.)
    """
    def __init__(self, reader, vcftype):
        self.vcfreader = reader
        self.vcftype = vcftype

class Writer:
    """
    Class to write the merged VCF file

    Parameters
    ----------
    out_path : str
          Name of the output VCF file
    samples : list of str
          IDs of samples to be included in the output
    command : str
          Command used to invoke this tool

    Attributes
    ----------
    vcf_writer : file
          Writeable file object to write VCF file to
    """
    
    def __init__(self, out_path, samples, command):
        self.vcf_writer = open(out_path, "w")
        self.vcf_writer.write('##fileformat=VCFv4.1\n')
        self.vcf_writer.write('##command=%s\n'%command)
        self.vcf_writer.write('##INFO=<ID=START,Number=1,Type=Integer,Description="First position in all alleles">\n')
        self.vcf_writer.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Last position in all alleles">\n')
        self.vcf_writer.write('##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Length of motif (repeat unit)">\n')
        self.vcf_writer.write('##INFO=<ID=RU,Number=1,Type=String,Description="Motif (repeat unit)">\n')
        self.vcf_writer.write('##INFO=<ID=METHODS,Number=1,Type=String,Description="Methods that attempted to genotype this locus (AdVNTR, EH, HipSTR, GangSTR)">\n')
        self.vcf_writer.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        self.vcf_writer.write('##FORMAT=<ID=SRC,Number=1,Type=String,Description="Source(s) of the merged call">\n')
        self.vcf_writer.write('##FORMAT=<ID=CERT,Number=1,Type=String,Description="Set to True if we are certain in the merged call">\n')
        self.vcf_writer.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')

    def WriteRecord(self, rcres):
        r"""
        Write a VCF record for the record cluster

        Parameters
        ----------

        rcres : recordcluster.RecordResolver
            resolver of record cluster 
        """
        if not rcres.resolved:
            common.WARNING("Warning: attempting to write record for unresolved record cluster")
            return

        FORMAT = ['GT', 'NCOPY', 'SRC','CERT','INPUTS']
        INFO_DICT = {'START': rcres.record_cluster.first_pos,
                     'END': rcres.record_cluster.last_end,
                     'PERIOD': len(rcres.record_cluster.canonical_motif),
                     'RU': rcres.record_cluster.canonical_motif,
                     'METHODS': "|".join([str(int(item)) for item in rcres.record_cluster.vcf_types])}
        
        res_pre_alleles = rcres.res_pas
        res_cert = rcres.res_cer
#        res_pre_alleles, res_cert = self.ResolveAllSampleCalls()
        out_rec = OutVCFRecord(res_pre_alleles, rcres.record_cluster)
        ####

        SAMPLE_DATA=[]
        for sample in rcres.record_cluster.samples:
            SAMPLE_DATA.append(':'.join(
                [out_rec.sample_to_GT[sample],
                out_rec.sample_to_NCOPY[sample],
                out_rec.sample_to_SRC[sample],
                str(res_cert[sample]),
                out_rec.sample_to_INPUTS[sample]
                ]
                ))
        
        ALTS = out_rec.alts
        if len(ALTS) == 0: ALTS.append(".")
        INFO = get_info_string(INFO_DICT)
        record_template = rcres.record_cluster.record_objs[0].cyvcf2_record
        # TODO remove record template and get information from pre allele
        self.vcf_writer.write('\t'.join([str(record_template.CHROM), 
            str(record_template.POS), 
            '.',
            out_rec.ref,
            ','.join(out_rec.alts),
            '.',
            '.',
            INFO,
            ':'.join(FORMAT),
            '\t'.join(SAMPLE_DATA)]) + '\n')

    def Close(self):
        r"""
        Close the writer file object
        """
        self.vcf_writer.close()

# TODO edit 
def get_info_string(data):
    out_recs = []
    for key in data:
        out_recs.append(str(key) + '=' + str(data[key]))
    return ';'.join(out_recs)

class OutVCFRecord:
    def __init__(self, res_pas, rc):
        self.pre_alleles = res_pas
        self.record_cluster = rc
        self.ref = ""
        self.alts = []
        self.sample_to_GT = {}
        self.sample_to_NCOPY = {}
        self.sample_to_SRC = {}
        self.sample_to_INPUTS = rc.PrintRawCalls()
        # Add INFO fields as well

        for sample in self.pre_alleles:
            GT_list = []
            NCOPY_list  = []
            SRC_list = []
            for pa in self.pre_alleles[sample]:
                if self.ref == "":
                    self.ref = pa.reference_sequence
                # TODO Fix assertion
                #assert self.ref == pa.ref_seq

                if pa.allele_sequence != self.ref:
                    # we have alt allele
                    if pa.allele_sequence not in self.alts:
                        self.alts.append(pa.allele_sequence)
                    GT_list.append(str(self.alts.index(pa.allele_sequence) + 1))
                    NCOPY_list.append(str(pa.allele_ncopy))
                else:
                    # ref allele
                    GT_list.append('0')
                    NCOPY_list.append(str(pa.reference_ncopy))
                for caller in pa.support:
                    if caller.name not in SRC_list:
                        SRC_list.append(caller.name)
            if len(GT_list) == 0:
                GT_list = ['.']
            if len(SRC_list) == 0:
                SRC_list = ['.']
            if len(NCOPY_list) == 0:
                NCOPY_list = ['.']
            self.sample_to_GT[sample] = '/'.join(GT_list)
            self.sample_to_NCOPY[sample] = ','.join(NCOPY_list)
            self.sample_to_SRC[sample] = ','.join(SRC_list)

class Readers:
    """
    Class to keep track of VCF readers being merged

    Parameters
    ----------
    vcfpaths : list of str
       List of paths to each of the input VCF files
    ref_genome : pyfaidx.Fasta
       Reference genome

    Attributes
    ----------
    ref_genome : pyfaidx.Fasta
       Reference genome
    vcfwrappers : list of VCFWrapper
       VCF wrappers for each input VCF 
    samples : list of str
       Samples shared by input VCF files
    """
    def __init__(self, vcfpaths, ref_genome):
        self.ref_genome = ref_genome
        self.vcfwrappers = []
        self.samples = []

        # first pass, determine shared samples across all vcf files
        for invcf in vcfpaths:
            vcffile = cyvcf2.VCF(invcf)
            hm = trh.TRRecordHarmonizer(vcffile)
            if len(self.samples) == 0:
                self.samples = vcffile.samples
            else:
                # setting sample list to overlap of sample lists
                self.samples = list(set(self.samples).intersection(set(vcffile.samples)))
        # Second pass, only load the shared samples
        for invcf in vcfpaths:
            vcffile = cyvcf2.VCF(invcf, samples = self.samples)
            hm = trh.TRRecordHarmonizer(vcffile)
            self.vcfwrappers.append(VCFWrapper(vcffile, hm.vcftype))
  
        # Get chroms and check if valid
        self.chroms = []
        for wrapp in self.vcfwrappers:
            if len(self.chroms) == 0:
                self.chroms = utils.GetContigs(wrapp.vcfreader)
            else:
                self.chroms = list(set(self.chroms) | set(utils.GetContigs(wrapp.vcfreader)))

        # Load current records
        self.current_tr_records = [trh.HarmonizeRecord(wrapper.vcftype, next(wrapper.vcfreader))
                                   for wrapper in self.vcfwrappers]
        if not self.areChromsValid():
            raise ValueError('Invalid CHROM detected in record.')
        self.done = all([item.vcfrecord is None for item in self.current_tr_records])
        self.is_min_pos_list = mergeutils.GetMinRecords(self.getCurrentRecordVCFRecs(),
                                                        self.chroms)
        self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos = \
            self.getCurrentRange()
        self.is_overlap_min = self.getOverlapMinRecords()

    def areChromsValid(self):
        r"""
        Check if chromosomes of current records are valid

        Returns
        -------

        is_valid : bool
           Equal to true if all are valid, else False
        """
        is_valid = True
        for r, wrapper in zip(self.current_tr_records, self.vcfwrappers):
            if r is None or r.vcfrecord is None:
                continue
            if r.vcfrecord.CHROM not in self.chroms:
                common.WARNING((
                                   "Error: found a record in file {} with "
                                    "chromosome '{}' which was not found in the contig list "
                                    "({})".format(wrapper.vcftype.name, r.vcfrecord.CHROM,
                                                  ", ".join(self.chroms))))
                is_valid = False
        return is_valid

    def getCurrentRecordVCFRecs(self):
        r"""
        Get list of the current records for each VCF reader

        Returns
        -------
        ret : list of trh.TRRecord
           List of VCF records for each reader
        """
        ret = []
        for item in self.current_tr_records:
            if item is None or item.vcfrecord is None:
                ret.append(None)
            else:
                ret.append(item.vcfrecord)
        return ret

    def getCurrentRange(self):
        r"""
        Get range (chrom, start, end) of current records

        Returns
        -------

        chrom : str
            Chromosome of current range
        start_pos : int
            Start position of the list of records
        end_pos : int
            Largest end position of the group of records
        """
        start_pos = None
        end_pos = -1
        chrom = None
        for i in range(len(self.is_min_pos_list)):
            if self.is_min_pos_list[i]:
                chrom = self.current_tr_records[i].vcfrecord.CHROM
                start_pos = self.current_tr_records[i].vcfrecord.POS
                end = start_pos + len(trh.HarmonizeRecord(self.vcfwrappers[i].vcftype,
                                                          self.current_tr_records[i].vcfrecord).ref_allele)
                if end > end_pos:
                    end_pos = end
        return chrom, start_pos, end_pos

    def getOverlapMinRecords(self):
        r"""
        Check of each record overlaps the min record

        Returns
        -------
        is_overlap_min : list of bool
           One entry for each VCF reader. Set to True
           if the current record for that reader overlaps
           the min record
        """
        # Set is_overlap_min for anything overlapping
        is_overlap_min = []
        for i in range(len(self.current_tr_records)):
            if self.current_tr_records[i] is None:
                is_overlap_min.append(False)
                continue
            if self.current_tr_records[i].vcfrecord.CHROM != self.cur_range_chrom:
                is_overlap_min.append(False)
                continue
            start_pos = self.current_tr_records[i].vcfrecord.POS
            end = start_pos + len(self.current_tr_records[i].vcfrecord.REF)
            if (self.cur_range_start_pos <= start_pos <= self.cur_range_end_pos) or \
                    (self.cur_range_start_pos <= end <= self.cur_range_end_pos):
                is_overlap_min.append(True)
            else:
                is_overlap_min.append(False)
        return is_overlap_min

    def getMergableCalls(self):
        r"""
        Determine which calls are mergeable
        Make one list of record clusters for each canonical motif

        Returns
        -------
        ov_region : recordcluster.OverlappingRegion
           contains VCF records in an overlapping region
        """
        record_cluster_list = []
        for i in range(len(self.current_tr_records)):
            if self.is_overlap_min[i] and self.current_tr_records[i] is not None:
                curr_ro = recordcluster.RecordObj(self.current_tr_records[i].vcfrecord, self.vcfwrappers[i].vcftype)
                canon_motif = utils.GetCanonicalMotif(curr_ro.hm_record.motif)
                added = False
                for rc in record_cluster_list:
                    if rc.canonical_motif == canon_motif:
                        rc.AppendRecordObject(curr_ro)
                        added = True
                if not added:
                    record_cluster_list.append(recordcluster.RecordCluster([curr_ro], self.ref_genome, \
                                                                           canon_motif, self.samples))
        ov_region = recordcluster.OverlappingRegion(record_cluster_list)
        return ov_region

    def goToNext(self):
        r"""
        Get next records for each reader
        """
        prev_records = self.current_tr_records
        new_records = []
        for idx, rec in enumerate(prev_records):
            if mergeutils.GetMinRecords(self.getCurrentRecordVCFRecs(), self.chroms)[idx]:
                try:
                    new_records.append(
                        trh.HarmonizeRecord(self.vcfwrappers[idx].vcftype,
                                            next(self.vcfwrappers[idx].vcfreader)))
                except StopIteration:
                    new_records.append(None)
            else:
                if prev_records[idx] is None:
                    new_records.append(None)
                else:
                    new_records.append(
                        trh.HarmonizeRecord(self.vcfwrappers[idx].vcftype,
                                            prev_records[idx].vcfrecord)
                    )
        self.current_tr_records = new_records

        # Update
        if not self.areChromsValid():
            raise ValueError('Invalid CHROM detected in record.')
        self.done = all([item is None or item.vcfrecord is None for item in self.current_tr_records])
        self.is_min_pos_list = mergeutils.GetMinRecords(self.getCurrentRecordVCFRecs(), self.chroms)
        self.cur_range_chrom, self.cur_range_start_pos, self.cur_range_end_pos = \
            self.getCurrentRange()
        self.is_overlap_min = self.getOverlapMinRecords()
