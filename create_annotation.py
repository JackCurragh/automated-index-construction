'''
This is script version 2 for annotation file creation

Inputs are: 

- GTF/GFF file
- FASTA file
- psuedo UTR length (default is 0)
- output file name

'''

import argparse
import sqlite3

import gffpandas.gffpandas as gffpd
import gzip
import numpy as np
import os
import pandas as pd
import tempfile

from Bio import SeqIO


def parse_fasta(fasta_path: str) -> dict:
    """
    Read in the transcriptome fasta file at the provided path
    and return a dictionary

    Inputs:
        fasta_path: Path to the transcriptome fasta file

    Outputs:
        transcript_dict: Dictionary containing the
                         transcript information
    """
    transcript_dict = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        transcript_dict[record.id.split('.')[0]] = record

    return transcript_dict


def is_gzipped(file_path: str) -> bool:
    """
    Checks whether the file is gzipped or not

    Inputs:
        file_path: Path to the file to be checked

    Outputs:
        True if gzipped, otherwise False
    """
    try:
        with open(file_path, 'rb') as f:
            # Read the first two bytes of the file
            header = f.read(2)

        # Check if the file starts with the gzip magic number (0x1f 0x8b)
        return header == b'\x1f\x8b'

    except IOError:
        # File not found or unable to open
        return False


def extract_transcript_id(attr_str):
    for attr in attr_str.split(";"):
        # Ensembl GFF3 support
        if attr.startswith("Parent=transcript:") \
                or attr.startswith("ID=transcript:"):
            return attr.split(":")[1]
        # Gencode GFF3 support
        elif attr.startswith("transcript_id="):
            return attr.split("=")[1]
        # Ensembl GTF support
        elif attr.startswith(" transcript_id "):
            return attr.split(" ")[2].replace('"', "")
    return np.nan


def parse_gff(gff_path: str) -> pd.DataFrame:
    """
    Read in the gff file at the provided path and return a dataframe

    Inputs:
        gff_path: Path to the gff file

    Outputs:
        gff_df: Dataframe containing the gff information
    """
    if is_gzipped(gff_path):
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_filepath = temp_file.name
            with gzip.open(gff_path, 'rt') as f:
                for line in f:
                    temp_file.write(line.encode())

        gff_df = gffpd.read_gff3(temp_filepath).df

        os.remove(temp_filepath)

    else:
        gff_df = gffpd.read_gff3(gff_path).df

    gff_df.loc[:, "transcript_id"] = gff_df["attributes"].apply(
        extract_transcript_id
        )

    cds_df = gff_df[gff_df["type"] == "CDS"]
    coding_tx_ids = cds_df["transcript_id"].unique()

    # subset GFF DataFrame to only include transcripts in the transcript_list
    gff_df = gff_df[gff_df["transcript_id"].isin(coding_tx_ids)]

    # Sort GFF DataFrame by transcript ID
    gff_df = gff_df.sort_values("transcript_id")

    return gff_df, coding_tx_ids


def parse_attributes(df: pd.DataFrame) -> dict:
    '''
    Parse the attributes into a dictionary form df

    Inputs:
        attributes: Attributes df

    Outputs:
        attributes_dict: Dictionary containing the attributes
    '''
    attributes_dict = {}
    for idx, row in df.iterrows():
        for attribute in row['attributes'].split(';'):
            if attribute != '':
                if attribute.split('=')[0] not in attributes_dict:
                    attributes_dict[attribute.split('=')[0]] = [attribute.split('=')[1]]
                else:
                    if attribute.split('=')[1] not in attributes_dict[attribute.split('=')[0]]:
                        attributes_dict[attribute.split('=')[0]].append(attribute.split('=')[1])

    return attributes_dict


def transcripts_table(gff_df: pd.DataFrame, transcripts: dict) -> pd.DataFrame:
    '''
    Handle the population of the transcripts table

    Transcript schema is as follows:

    Inputs:
        gff_df: GFF dataframe
        transcripts: Dictionary containing the transcript information

    Outputs:
        transcripts_df: Dataframe containing the transcripts table
    '''
    transcripts_df = pd.DataFrame(
        columns=['transcript', 'gene', 'length', 'cds_start',
                 'cds_stop', 'sequence', 'strand', 'stop_list',
                 'start_list', 'exon_junctions', 'tran_type',
                 'gene_type', 'principal', 'version', 'gc',
                 'five_gc', 'cds_gc', 'three_gc', 'chrom']
                 )

    # group by transcript id
    gff_df = gff_df.groupby('transcript_id')

    # iterate over each transcript
    for transcript_id, transcript_df in gff_df:

        attributes = parse_attributes(transcript_df)
        transcripts_df_entry = {
            "transcript": transcript_id.split('.')[0],
            "gene": attributes['gene_name'][0],
            "length": transcripts[transcript_id.split('.')[0]].__len__(),
            # "cds_start": transcript_df['CDS'].iloc[0].start,
        }

        print(transcript_id, transcript_df)
        print(transcripts[transcript_id.split('.')[0]])

        print()
        for elem in transcripts_df_entry:
            print(elem, transcripts_df_entry[elem])
        break


def main(args):
    '''
    Run Sqlite creation
    '''
    print("Parsing FASTA file...")
    transcripts = parse_fasta(args.fasta)
    print("Done")
    print("Parsing GFF/GTF file...")
    gff_df, coding_tx_ids = parse_gff(args.gtf)
    print("Done")

    transcripts_table(gff_df, transcripts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create annotation file for Trips-Viz")
    parser.add_argument('-g', '--gtf', help='GTF/GFF file', required=True)
    parser.add_argument('-f', '--fasta', help='FASTA file', required=True)
    parser.add_argument('-p', '--pseudo', help='Pseudo UTR length', default=0)
    parser.add_argument('-o', '--output', help='Output file name', required=True)

    args = parser.parse_args()
    main(args)

'''
python create_annotation.py \
        -g organism_sqlites/homo_sapiens_gencode25/gencode.v25.annotation.gff3 \
        -f /home/jack/projects/index-construction/organism_sqlites/homo_sapiens_gencode25/gencode.v25.transcripts.fa \
        -o test.sqlite
'''