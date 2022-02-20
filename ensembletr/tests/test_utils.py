from .. import utils
import pytest

def test_EHScore():
	score = utils.GetEHScore("14-14/14-14", "14/14", 2)
	assert(score == 1)