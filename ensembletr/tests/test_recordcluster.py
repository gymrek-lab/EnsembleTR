from .. import recordcluster
import trtools.utils.tr_harmonizer as trh

import os
import pytest
import cyvcf2

def test_RecordObj(vcfdir):
	vcffile = os.path.join(vcfdir, "gangstr-chr20.vcf.gz")
	reader = cyvcf2.VCF(vcffile)
	record = next(reader)
	ro_gangstr = recordcluster.RecordObj(record, trh.VcfTypes.gangstr)

	assert(ro_gangstr.canonical_motif == "A")
