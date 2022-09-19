import subprocess
import pandas as pd
import os

from Bio import SeqIO


def build_index(command, fasta_path, outfile_base):
    '''
    build the index. 
    command can be specified for different programs (bowtie, star, bwa)
    '''
    subprocess.run([command, fasta_path, outfile_base])


def unzip(zipped_path):
    '''
    decompress the file at the given path
    '''
    subprocess.run(['gzip', '-d', zipped_path])



def get_fasta_file(directory_path, fasta_url):
    '''
    wget the fasta file from the url and place it at the directory path
    '''
    subprocess.run(['wget', fasta_url, '-P', directory_path])


def process_nc(fasta_path, output_path):
    '''
    split the nc fasta file into a tRNA and rRNA.
    '''

    rRNAs = []
    tRNAs = []
    uncompressed_fasta_path = '.'.join(fasta_path.split('.')[:-1])
    if not os.path.isfile(uncompressed_fasta_path):
        unzip(uncompressed_fasta_path)

    for seq_record in SeqIO.parse(uncompressed_fasta_path, "fasta"):
        if 'tRNA' in seq_record.description:
            tRNAs.append(seq_record)
        
        if 'rRNA' in seq_record.description:
            rRNAs.append(seq_record)

    with open(output_path + '/' + "rRNA.fasta", "w") as output_handle:
        SeqIO.write(rRNAs, output_handle, "fasta")

    with open(output_path + '/' + "tRNA.fasta", "w") as output_handle:
        SeqIO.write(tRNAs, output_handle, "fasta")




def bowtie_main(index_metadata):
    '''
    prepare the files and directories for bowtie index construction and execute
    '''
    metadata = pd.read_csv(index_metadata, sep="\t")
    index_types = list(metadata.index_type.unique())[::-1]
    for variety in index_types:

        if not os.path.isdir('indices/bowtie/' + variety):
            os.mkdir('indices/bowtie/' + variety)

        variety_df = metadata.query(f'index_type == "{variety}"')
        for row in variety_df.iterrows():
            index = row[1]['index_name']
            fasta_url = row[1]['fasta_url']
            output_path = 'indices/bowtie/' + variety + '/' + index

            if not os.path.isdir(output_path):
                os.mkdir(output_path)

            fasta_path = output_path + '/' + fasta_url.split('/')[-1]

            if fasta_url != 'None':
                if not os.path.isfile(fasta_path):
                    get_fasta_file(output_path, fasta_url)
            else: 
                continue
            basename = index + "_" + variety
            if variety == 'ncRNA':
                process_nc(fasta_path, output_path)

                for nc_type in ['rRNA', 'tRNA']:
                    fasta_path = output_path + '/' + nc_type + '.fasta'
                    basename = index + "_" + nc_type
                    build_index('bowtie-build', fasta_path, output_path + '/' + basename )

            else:
                build_index('bowtie-build', fasta_path, output_path + '/' + basename )
            # return('hi')

# build_index('bowtie-build', '/home/jack/projects/tools_for_Galaxy/test-data/Yeast_rRNA.fasta.fasta', 'test/test')
# bowtie_main('bowtie_index_construction.tsv')

process_nc('indices/bowtie/ncRNA/mus_musculus_m28/Mus_musculus.GRCm39.ncrna.fa.gz', 'indices/bowtie/ncRNA/mus_musculus_m28')