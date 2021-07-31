from math import nan
from trtools.statSTR import GetHWEP
import numpy as np

class LocusReliability:
    def __init__(self, hm_record):
        self.hm_record = hm_record
        # TODO do something more appropriate to average over all samples
        self.hwep = min(GetHWEP(hm_record))


    # Score between zero and 1 (1 more confident, 0 less confident)
    def GetLocusScore(self):
        if self.hwep == 0 or self.hwep == 1:
            return self.hwep
        if self.hwep is np.nan:
            return 0
        log_hwep = np.log(self.hwep)
        neg_log = -1 * log_hwep
        one_over_neg_log = 1 / neg_log
        multiplier = 3
        ret_value = multiplier * one_over_neg_log

        if ret_value >= 1:
            return 1
        elif ret_value < 0:
            raise ValueError("Negative locus score!", ret_value)
        else:
            return ret_value
        

