#!/usr/bin/env python
from __future__ import print_function

import re
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
import logging

log = logging.getLogger(__name__)


NAMESPACE = "DRAGEN FastQC"


class DragenFastqcMetrics(BaseMultiqcModule):
    def add_fastqc_metrics(self):
        data_by_sample = dict()
        for f in self.find_log_files("dragen/fastqc_metrics"):
            s_name, data = parse_fastqc_file(f)
            s_name = self.clean_s_name(s_name, f)
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="stats")
            data_by_sample[s_name] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            return set()

        # Write data to file
        self.write_data_file(data_by_sample, "fastqc")

        headers = OrderedDict()
        headers["avg_gc_content_percent"] = {
            "title": "% GC",
            "description": "Average % GC Content",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "Set1",
            "format": "{:,.0f}",
        }
        self.general_stats_addcols(data_by_sample, headers, namespace=NAMESPACE)
        return data_by_sample.keys()


def parse_fastqc_file(f):
    """
    PREFIX.fastqc_metrics.csv
    POSITIONAL BASE MEAN QUALITY,Read1,ReadPos 4 A Average Quality,36.59
    """

    s_name = re.search(r"(.*)\.fastqc_metrics.csv", f["fn"]).group(1)

    data = defaultdict(dict)
    gcreads = 0
    totalreads = 0
    for line in f["f"].splitlines():
        if(line.split(",")[0] != "READ GC CONTENT"):
            continue

        qualtype, read, metric, stat = line.split(",")

        try:
            reads = float(stat)      
            percent = float(float(metric.split('%')[0]) / 100)
            gcreads += (reads * percent)
            totalreads += reads
        except ValueError:
            pass

    data["avg_gc_content_percent"] = gcreads / totalreads * 100
    return s_name, data
