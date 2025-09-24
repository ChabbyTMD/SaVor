include: "rules/common.smk"
include: "rules/fastq_ingestion.smk"
include: "rules/alignment.smk"
include: "rules/lumpy.smk"
include: "rules/delly.smk"
include: "rules/wham.smk"
include: "rules/svcallprocess.smk"
include: "rules/benchmark.smk"

import sys
import os
from pathlib import Path
import pandas as pd

configfile: "config/config.yaml"

wildcard_constraints:
    window=r"\d+"

samples = parse_sample_sheet(config)
REFGENOME = samples['refGenome'].unique().tolist()

rule all:
    input:
        svArcher_output,