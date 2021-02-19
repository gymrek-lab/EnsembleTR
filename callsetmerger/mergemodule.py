#!/usr/bin/env python3

"""
Merge Module includes code and classes for merging record clusters.
Work in progress

"""
from callsetmerger.recordcluster import ClusterGraph




# class MergeOverlappingRegion:
#     def __init__(self, overlapping_region):
#         self.ov_region = overlapping_region
#         self.rc_merge_objects = []
#         for rc in overlapping_region.RecordClusters:
#             self.rc_merge_objects.append(RecordClusterMerger(rc))

#     def GetVCFLines(self):
#         ret = []
#         for rc_merge_obj in self.rc_merge_objects:
#             line = rc_merge_obj.GetVCFLine()
#             if line is not None:
#                 ret.append(line)
#         return ret



class RecordClusterMerger:
    def __init__(self, rc, samples):
        self.record_cluster = rc
        self.graph = ClusterGraph(rc)
        self.samples = samples
        for sample in self.samples:
            samp_gt = self.record_cluster.GetSampleCall(sample)
            resolved_gt = self.graph.ResolveGenotypes(samp_gt)
            print(sample, resolved_gt)

    def GetVCFLine(self):
        line = ''
        # for sample in samples:
            # Get call for sample
            # if len(self.graph.vcf_types) == 1:
            #     # Only one callser
            #     pass
            # else:
            #     # More than one caller
            #     if self.graph.GetConfusionScore() == 1.0:
            #         # Simple graph with no confusion
            #         pass
            #     else:
            #         # More complex graph with confusion
            #         pass
        return line


# For each sample, we have 1 object that we pass the graph to -> call for that sample
