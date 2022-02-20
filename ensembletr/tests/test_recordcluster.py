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
	assert(ro_gangstr.GetCalledAlleles() == set([0, 1, 2, 3]))
	# Test sample 0 - hom ref
	assert(ro_gangstr.GetROSampleCall(0)[0] == 0)
	assert(ro_gangstr.GetROSampleCall(0)[1] == 0)
	assert(ro_gangstr.GetSampleString(0) == "gangstr=15.0,15.0")
	# Test sample 9 - het
	assert(ro_gangstr.GetROSampleCall(9)[0] == 1)
	assert(ro_gangstr.GetROSampleCall(9)[1] == 0)
	assert(ro_gangstr.GetSampleString(9) == "gangstr=14.0,15.0")
