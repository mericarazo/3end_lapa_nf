#!/usr/bin/env python

from lapa.result import LapaResult

dir = 'lapa_out'
result = LapaResult(dir)

counts = result.attribute('count')
counts.to_csv(f'{dir}/lapa_counts.tsv', index=True, sep="\t")

usage = result.attribute('usage')
usage.to_csv(f"{dir}/lapa_usage.tsv", index=True, sep="\t")