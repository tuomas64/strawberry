import vcf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Function to map heterozygous alleles to IUPAC ambiguity codes
def get_ambiguity_code(alleles):
    iupac_codes = {
        ('A', 'G'): 'R', ('G', 'A'): 'R',
        ('C', 'T'): 'Y', ('T', 'C'): 'Y',
        ('G', 'C'): 'S', ('C', 'G'): 'S',
        ('A', 'T'): 'W', ('T', 'A'): 'W',
        ('G', 'T'): 'K', ('T', 'G'): 'K',
        ('A', 'C'): 'M', ('C', 'A'): 'M',
        ('A', 'A'): 'A', ('T', 'T'): 'T',
        ('C', 'C'): 'C', ('G', 'G'): 'G'
    }
    return iupac_codes.get(alleles, 'N')  # Default to 'N' for unknown or complex alleles

def vcf_to_fasta_with_ambiguity(vcf_file, output_fasta):
    # Open the VCF file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    
    # Create a dictionary to hold sequences for each sample
    sequences = {sample: [] for sample in vcf_reader.samples}

    # Parse through each record (variant site)
    for record in vcf_reader:
        ref = record.REF
        alt = record.ALT[0]  # Assuming only one ALT allele; expand if needed
        
        for sample in record.samples:
            genotype = sample['GT']
            
            if genotype == '0/0':  # Homozygous reference
                sequences[sample.sample].append(ref)
            elif genotype == '1/1':  # Homozygous variant
                sequences[sample.sample].append(alt)
            elif genotype == '0/1' or genotype == '1/0':  # Heterozygous
                # Get IUPAC ambiguity code
                ambiguity_code = get_ambiguity_code((ref, alt))
                sequences[sample.sample].append(ambiguity_code)
            else:
                sequences[sample.sample].append('N')  # Handle missing data or other genotypes

    # Write the sequences to a FASTA file
    fasta_records = []
    for sample, seq_list in sequences.items():
        # Combine the sequence list into a single string
        sequence_str = ''.join(seq_list)
        # Create a SeqRecord for each sample
        fasta_records.append(SeqRecord(Seq(sequence_str), id=sample, description=""))
    
    # Write the records to a FASTA file
    SeqIO.write(fasta_records, output_fasta, 'fasta')

# Example usage
vcf_file = "input_file.vcf"
output_fasta = "output_file.fasta"
vcf_to_fasta_with_ambiguity(vcf_file, output_fasta)
