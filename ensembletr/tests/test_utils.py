from ensembletr import utils

def test_EHScore():
	score = utils.GetEHScore("14-14/14-14", "14/14", 2)
	assert(score == 1)

def test_motif_equality():
	assert(utils.MotifEquality("AGRG", 'AGGG') == True)
	assert (utils.MotifEquality("AGG", 'GGR') == True)
	assert (utils.MotifEquality("GBG", 'GGG') == True)
	assert (utils.MotifEquality("GHG", 'GGG') == False)
	assert (utils.MotifEquality("NHG", 'GAC') == True)
	assert (utils.MotifEquality("NHG", 'HAC') == True)
	assert (utils.MotifEquality("MHG", 'AGG') == False)

def test_all_possible_seqs():
	assert(utils.GetAllPossibleSequences("R") == ['A', 'G'])
	assert (utils.GetAllPossibleSequences("M") == ['A', 'C'])
	assert (utils.GetAllPossibleSequences("TM") == ['TA', 'TC'])
	assert (sorted(utils.GetAllPossibleSequences("BG")) == sorted(['CG', 'GG','TG']))
	assert (len(utils.GetAllPossibleSequences("BGRR")) == 12)
	assert (len(utils.GetAllPossibleSequences("BGBB")) == 27)
	assert ('CGGA' in utils.GetAllPossibleSequences("BGDH"))
	assert (utils.GetAllPossibleSequences('AGTT') == ['AGTT'])

test_all_possible_seqs()