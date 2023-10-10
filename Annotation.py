import subprocess
from pathlib import Path

import polars as pl


def is_gtf_header_line(line):
    return line.startswith("#")


class Annotation:
    def __init__(self, gtf_file):
        self.gtf_file = gtf_file
        self.df = None

        if self.gtf_file_exists():
            self.gtf_to_df()
            self.filter_only_genes()
            self.sort_df_by_coordinates()

    def gtf_file_exists(self):
        return self.gtf_file.exists()

    def gtf_count_header_lines(self):
        no_header_lines = 0
        with open(self.gtf_file, "r") as gtf:
            while is_gtf_header_line(gtf.readline()):
                no_header_lines += 1
        return no_header_lines

    def gtf_to_df(self):
        df = pl.read_csv(
            self.gtf_file,
            has_header=False,
            separator="\t",
            skip_rows=self.gtf_count_header_lines(),
            new_columns=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
            dtypes={
                "seqname": str,
                "source": str,
                "feature": str,
                "start": pl.Int64,
                "end": pl.Int64,
                "score": str,
                "strand": str,
                "frame": str,
                "attribute": str,
            },
        )
        self.df = df
        return self

    def sort_df_by_coordinates(self):
        self.df = self.df.sort(["seqname", "start", "end"])
        return self

    def filter_only_genes(self):
        self.df = self.df.filter(pl.col("feature") == "gene")
        return self
