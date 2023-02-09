import cyvcf2
from cyvcf2 import VCF, Writer
from itertools import islice, count
from collections import defaultdict
import numpy as np
import math
import sys
import time


file_name = sys.argv[1] # Input vcf file
vcf = VCF(file_name, strict_gt=True)
samples = vcf.samples


def merge(merging_list):
    # extracting the coordinates of maximum mutual string \
    # calling trim function to trim alleles \
    # update format fields for the merging list.
    max_start = 0
    min_end = 1000000000000
    for record in merging_list:
        if int(record.POS) > max_start:
            max_start = int(record.POS)
        if (int(record.POS) + len(record.REF)) - 1 < min_end:
            min_end = (int(record.POS) + len(record.REF)) - 1
    alleles, allele_map, ref_allele = trim(merging_list, max_start, min_end)
    updated_format = update_format(alleles, allele_map, merging_list, ref_allele)
    assert(max_start > 0)
    assert(min_end > 0)
    assert(min_end < 1000000000000)
    return updated_format, alleles, ref_allele, max_start


def trim(merging_list, start, end):
    # trimming alleles according to given coordinates. return the new list \
    # of alleles,  mapping of initial alleles to their \
    # corrected allele, corrected ref allele and corrected POS.

    alleles = set()
    refs = set()
    allele_map = {}
    for i in range(len(merging_list)):
        record = merging_list[i]
        start_diff = start - record.POS 
        end_diff = record.POS + len(record.REF) - 1 - end
        trimmed_ref = record.REF[start_diff:len(record.REF) - end_diff]
        refs.add(trimmed_ref)
        assert(start_diff >= 0)
        assert(end_diff >= 0)
        for allele in record.ALT:
            trimmed_allele = allele[start_diff:len(allele) - end_diff]
            alleles.add(trimmed_allele)
            allele_map[(allele,i)] = trimmed_allele
    assert(len(refs) == 1)  # All the ref alleles should be same after trimming.
    ref_allele = list(refs)[0]
    alleles = list(alleles)
    for allele in alleles:
           if allele == ref_allele:
                alleles.remove(allele)
                break
    alleles.sort()
    alleles.insert(0,ref_allele)
    assert(len(ref_allele) > 0)
    return alleles, allele_map, ref_allele


def update_format(alleles, allele_map, merging_list, ref_allele):

    #  update format field for the merging records. GT and GB is updated \
    #  according to the new list of alleles and mapping of old alleles. \
    #  other values are copied.
    updated_format = {}
    for sample in samples:
        updated_format[sample] = {}
    for j in range(len(merging_list)):
        record = merging_list[j]
        format_data = {}
        for format_field in record.FORMAT[1:]:
            format_data[format_field] = record.format(format_field)
        for i in range(len(samples)):
            
            # updating GT and GB
            gt = record.genotypes[i]
            sample = samples[i]
            
            if gt[0] == -1:
                continue
            gt1 = gt[0]
            gt2 = gt[1]
            phasing = gt[2]
            new_gt1 = new_gt2 = new_gb1 = new_gb2 = 0
            if gt1 != 0:
                new_allele = allele_map[(record.ALT[gt1-1],j)]
                if new_allele == '':  # empty sequence means that allele got fully trimmed
                    continue
                new_gt1 = alleles.index(new_allele)
                new_gb1 = len(new_allele) - len(ref_allele)
            if gt2 != 0:
                new_allele = allele_map[(record.ALT[gt2-1],j)]
                if new_allele == '':
                    continue
                new_gt2 = alleles.index(new_allele)
                new_gb2 = len(new_allele) - len(ref_allele)
            if phasing:
                updated_format[sample]['GT'] = str(new_gt1) + "|" + str(new_gt2)
            else:
                updated_format[sample]['GT'] = str(new_gt1) + "/" + str(new_gt2)
            updated_format[sample]['GB'] = str(new_gb1) + "|" + str(new_gb2)

            # updating other format fields
            for format_field in record.FORMAT[2:]:
                if format_data[format_field][i] == -2147483648:
                    updated_format[samples[i]][format_field] = "."
                    continue
                if format_data[format_field][i] != ".":
                    if type(format_data[format_field][i]) == np.ndarray:
                        if math.isnan(format_data[format_field][i][0]):
                            updated_format[samples[i]][format_field] = "."
                        else:
                            updated_format[samples[i]][format_field] = format_data[format_field][i][0]
                    else:
                        updated_format[samples[i]][format_field] = format_data[format_field][i]
    
    # check if any of samples had no call at all
    for sample in updated_format:
        if not updated_format[sample]:
            for format_field in merging_list[0].FORMAT:
                updated_format[sample][format_field] = "."
    return updated_format


def get_record_str(record):
    # convert unchanged vcf record to string to be written
    ALL_SAMPLE_DATA = []
    format_data = {}
    for format_field in record.FORMAT[1:]:
        format_data[format_field] = record.format(format_field)
    for i in range(len(samples)):
        sample_data = []
        sample_GT = record.genotypes[i]
        if sample_GT[0] == -1:
            sample_data.append(".")
        else:
            if sample_GT[2]:
                sample_data.append(str(sample_GT[0]) + "|" + str(sample_GT[1]))
            else:
                sample_data.append(str(sample_GT[0]) + "/" + str(sample_GT[1]))

        for format_field in record.FORMAT[1:]:
            if format_data[format_field][i] == -2147483648:
                sample_data.append(".")
                continue
            if type(format_data[format_field][i]) == np.ndarray:
                if math.isnan(format_data[format_field][i][0]):
                    sample_data.append(".")
                else:
                    sample_data.append(str(format_data[format_field][i][0]))
            else:
                sample_data.append(format_data[format_field][i])
        sample_data = ":".join([data for data in sample_data])
        ALL_SAMPLE_DATA.append(sample_data)
    INFO = {'START': str(record.POS), 'END': str(record.POS + len(record.REF) - 1), 'PERIOD':str(record.INFO['PERIOD'])}
    output = '\t'.join([record.CHROM, str(record.POS), record.ID,
            record.REF, ",".join(record.ALT), ".", ".", ";".join(["%s=%s"%(key, INFO[key]) for key in INFO]),
            ':'.join(record.FORMAT),
            '\t'.join(ALL_SAMPLE_DATA)]) + '\n'
    return output


def get_updated_record_str(updated_format, alleles, ref_allele, record, pos):
    # convert a merged record to string to be written.
    ALL_SAMPLE_DATA = []
    allele_string = ",".join(alleles[1:])
    if len(alleles) == 1 or (len(alleles) == 2 and alleles[1] == ""):
        allele_string = "."
    for sample in updated_format:
        sample_data = updated_format[sample]
        sample_data = ":".join([str(sample_data[key]) for key in sample_data])
        ALL_SAMPLE_DATA.append(sample_data)
    INFO = {'START': str(pos), 'END': str(pos + len(ref_allele) - 1), 'PERIOD':str(record.INFO['PERIOD'])}
    return '\t'.join([record.CHROM, str(pos), record.ID,
            ref_allele, allele_string, ".", ".",
            ";".join(["%s=%s"%(key, INFO[key]) for key in INFO]), ':'.join(record.FORMAT),
            '\t'.join(ALL_SAMPLE_DATA)]) + '\n'

def vcf_window(it, n):
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def main():
    corrected_addr = sys.argv[2]
    vcf_writer = open(corrected_addr, "w")
    vcf_writer.write(vcf.raw_header)

    window = vcf_window(vcf, 80)  # Looking up to 80 records further to find conflicts.
    skip = 0
    for records in window:
        if skip > 0:
            skip -= 1
            continue
        merging_list = []
        current_record = records[0]
        current_id = current_record.ID
        merging_list.append(current_record)
        next_index = 1
        while True:
            if next_index >= len(records):
                break
            next_record = records[next_index]
            if next_record.ID == current_id:
                merging_list.append(next_record)
            else:
                break
            next_index += 1
        if len(merging_list) == 1:
            vcf_writer.write(get_record_str(merging_list[0]))
        else:
            print([record.ID for record in merging_list])
            updated_format, alleles, ref_allele, pos = merge(merging_list)
            vcf_writer.write(get_updated_record_str(updated_format, alleles, ref_allele, merging_list[0], pos))
            skip = len(merging_list) - 1

    # last window is not fully iterated
    index = skip + 1
    while True:
        if index >= len(records):
            break
        merging_list = []
        current_record = records[index]
        current_id = current_record.ID
        merging_list.append(current_record)
        next_index = index + 1
        while True:
            if next_index >= len(records):
                break
            next_record = records[next_index]
            if next_record.ID == current_id:
                merging_list.append(next_record)
            else:
                break
            next_index += 1
        if len(merging_list) == 1:
            vcf_writer.write(get_record_str(merging_list[0]))
            index += 1
        else:
            print([record.ID for record in merging_list])
            updated_format, alleles, ref_allele, pos = merge(merging_list)
            vcf_writer.write(get_updated_record_str(updated_format, alleles, ref_allele, merging_list[0], pos))
            index += len(merging_list)

    vcf_writer.close()

main() 
