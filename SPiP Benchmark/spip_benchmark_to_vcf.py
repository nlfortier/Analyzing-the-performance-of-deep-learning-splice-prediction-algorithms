import os
import re
import sys

import pygautil

def findLatestTsf(directory, filePrefix, substring):
    latestFile = None
    latestTime = 0

    for fileName in os.listdir(directory):
        # Check if the file starts with the required prefix
        if fileName.startswith(filePrefix) and substring in fileName:
            filePath = os.path.join(directory, fileName)
            
            if not os.path.isfile(filePath):
                continue
            
            modTime = os.path.getmtime(filePath)
            
            # Update latest_file if this file is more recent
            if modTime > latestTime:
                latestTime = modTime
                latestFile = filePath
    
    return latestFile

DEPENDENCY_DIR = "/mnt/datarepo/raw/annotations/Homo sapiens, GRCh_38/Genomenon/Dependencies/"
GENE_ANNOTATION_DIR = "/mnt/datarepo/mirror/annotations/Genes and Regulation/"
SYSTEM_DATA_PATH = "/home/fortier/data/VarSeq_Sys_Data/"
REFSEQ_GENES_38_TSF = findLatestTsf(GENE_ANNOTATION_DIR, "RefSeqGenes", "GRCh_38")
REFERENCE_SEQUENCE_38 = DEPENDENCY_DIR + "ReferenceSequenceV2-NCBI_GRCh_38_Homo_sapiens.tsf"
REF_READER_38 = pygautil.openForRead(REFERENCE_SEQUENCE_38)

pygautil.setCoordSysId("GRCh_38,Chromosome,Homo sapiens")
pygautil.setReferenceSequenceSource(REFERENCE_SEQUENCE_38)
pygautil.setTranscriptSource(REFSEQ_GENES_38_TSF)
pygautil.setSystemDataPath(SYSTEM_DATA_PATH)

VCF_HEADER="""##fileformat=VCFv4.1
##source=VarSeqV2.5.0
##fileDate=20240108
##reference=GRCh_38,Chromosome,Homo sapiens
##contig=<ID=1,length=248956422>
##contig=<ID=2,length=242193529>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##contig=<ID=MT,length=16569>
##INFO=<ID=GeneName,Number=1,Type=String,Description="Gene name">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""

# GRCh38 RefSeq (NC_) accession numbers for human chromosomes
_CHR_TO_NC_GRCH38 = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10",
    "MT": "NC_012920.1",
    "M": "NC_012920.1",  # alternative mitochondrial designation
}


GENE_NAME_MAP = None
def load_gene_coordinates(tsfFile):
    gene_map = {}
    refseqTrack = pygautil.openForRead(tsfFile)
    coverageSpace = refseqTrack.coverageSpace()
    for chr in coverageSpace:
        refseqTrack.readChr(chr[0], ['Start', 'Stop', 'Gene Name'])
        while refseqTrack.hasNext():
            data = refseqTrack.next()
            geneName = data[-1]

            if geneName not in gene_map:
                gene_map[geneName] = [chr[0], data[0], data[1]]

            start = data[0]
            stop = data[1]
            if start < gene_map[geneName][1]:
                gene_map[geneName][1] = start

            if stop > gene_map[geneName][2]:
                gene_map[geneName][2] = stop

    return gene_map

def get_gene_coordinates(geneName):
    global GENE_NAME_MAP
    if GENE_NAME_MAP is None:
        GENE_NAME_MAP = load_gene_coordinates(REFSEQ_GENES_38_TSF)

    return GENE_NAME_MAP.get(geneName, None)

def chromosome_to_nc(chr_name: str):
    if chr_name is None:
        return None

    if chr_name.startswith("chr"):
        chr_name = chr_name.replace("chr", "")

    key = str(chr_name).strip().upper()
    if key.startswith("CHR"):
        key = key[3:].lstrip()
    return _CHR_TO_NC_GRCH38.get(key)


# Genomic sort order for chromosomes (1-22, X, Y, MT); unknown chrs sort last
_CHR_ORDER = {str(i): i for i in range(1, 23)}
_CHR_ORDER.update({"X": 23, "Y": 24, "MT": 25, "M": 25})


def _chr_sort_key(chr_name: str):
    """Return a comparable key for genomic chromosome order."""
    s = str(chr_name).strip().upper()
    if s.startswith("CHR"):
        s = s[3:].lstrip()
    return _CHR_ORDER.get(s, 999)


def sort_variants_by_genome(vcf_lines: list) -> list:
    """
    Sort a list of VCF line strings in genomic order: by chromosome (1-22, X, Y, MT),
    then by POS numerically. Modifies the list in place and returns it.
    Each line is tab-separated with CHROM in column 0 and POS in column 1.
    """
    def key(line):
        parts = line.split("\t", 2)
        chrom = parts[0] if len(parts) > 0 else ""
        pos = int(parts[1], 10) if len(parts) > 1 else 0
        return (_chr_sort_key(chrom), pos)
    vcf_lines.sort(key=key)
    return vcf_lines


_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str):
    return seq.translate(_COMPLEMENT)[::-1]


def _normalize_range_in_rest(rest: str) -> str:
    """Swap start and stop in g.start_stop so that start <= stop."""
    m = re.search(r"(\d+)_(\d+)", rest)
    if m and int(m.group(1)) > int(m.group(2)):
        start, stop = m.group(2), m.group(1)
        return rest[: m.start(1)] + start + "_" + stop + rest[m.end(2) :]
    return rest


def revcomp_alternate_in_rest(rest: str):
    # substitution ref>alt: reverse complement both ref and alt
    m = re.search(r"([ACGT]+)>([ACGT]+)$", rest, re.IGNORECASE)
    if m:
        return rest[: m.start(1)] + reverse_complement(m.group(1)) + ">" + reverse_complement(m.group(2))
    # delinsSEQ
    m = re.search(r"(delins)([ACGT]+)$", rest, re.IGNORECASE)
    if m:
        return rest[: m.start(2)] + reverse_complement(m.group(2))
    # insSEQ (negative lookbehind so we don't match "ins" in "delins")
    m = re.search(r"(?<!del)(ins)([ACGT]+)$", rest, re.IGNORECASE)
    if m:
        return rest[: m.start(2)] + "ins" + reverse_complement(m.group(2))
    return rest


def _strip_operator_trailing_digits(rest: str) -> str:
    """Remove invalid trailing digits after operators (e.g. del5 -> del, ins3 -> ins)."""
    # Match operator + digits only at end; order matters (delins before del/ins)
    for pattern in (r"(delins)(\d+)$", r"(del)(\d+)$", r"(?<!del)(ins)(\d+)$", r"(dup)(\d+)$", r"(inv)(\d+)$"):
        rest = re.sub(pattern, r"\1", rest, flags=re.IGNORECASE)
    return rest


def fix_hgvs(hgvs: str, chr_name: str, strand: str = "+"):
    if "delin" in hgvs and "delins" not in hgvs:
        return hgvs.replace("delin", "delins")
    
    if "c. " in hgvs:
        return hgvs.replace("c. ", "c.")

    parts = hgvs.strip().split(":")
    rest = "".join(parts[1:]) if len(parts) > 1 else ""
    rest = _strip_operator_trailing_digits(rest)
    # For g. notation only: use NC ref, normalize range, optional revcomp
    if rest.startswith("g."):
        nc_id = chromosome_to_nc(chr_name)
        if nc_id is None:
            return None
        rest = _normalize_range_in_rest(rest)
        if strand == "-" and rest:
            rest = revcomp_alternate_in_rest(rest)
        ref = nc_id
    else:
        ref = parts[0] if parts else ""
    return f"{ref}:{rest}" if rest else ref

def converted_coordinates_are_valid(coordinates):
    if len(coordinates) > 0 and len(coordinates[0]) > 0 and type(coordinates[0]) == list:
        start = coordinates[0][1]
        stop = coordinates[0][2]
        refalt = coordinates[0][3]
        return start >= 0 and stop >= 0 and len(refalt) < 100 and 'does not match' not in coordinates[0][-1]
    
    return False

def coordinates_mismatch_reference(converted_coordinates):
    return 'does not match' in converted_coordinates[0][-1]

def extract_number(text):
    return ''.join([c for c in text if c.isdigit()])

def extract_hgvs_refalt(hgvs):
    rest = "".join([c for c in hgvs.split(".")[-1] if not c.isdigit()])
    return rest.split(">")

def liftover_hgvs_coordinates(hgvs: str, chr_name: str, coordinates):
    position = coordinates[0][1]
    liftover = pygautil.liftover(chr_name.replace("chr", ""), position, position, "", "GRCh_38")

    if len(liftover) > 1:
        old_start_1_based = position + 1
        new_start_1_based = liftover[1] + 1
        hgvs = hgvs.replace(str(old_start_1_based), str(new_start_1_based))
        coordinates = pygautil.convertPdot(hgvs)

        if coordinates_mismatch_reference(coordinates) and ">" in hgvs:
            # probably mismatch due to reverse compliment issue
            ref, alt = extract_hgvs_refalt(hgvs)
            refc = reverse_complement(ref)
            altc = reverse_complement(alt)
            hgvs = hgvs.replace(f"{ref}>{alt}", f"{refc}>{altc}")
            coordinates = pygautil.convertPdot(hgvs)

    if converted_coordinates_are_valid(coordinates):
        return coordinates

    return []

def overlaps_gene(position: int, gene: str):
    gene_coordinates = get_gene_coordinates(gene)
    if gene_coordinates:
        return position >= gene_coordinates[1] and position <= gene_coordinates[2]

    return True

def convert_hgvs(hgvs: str, chr_name: str, position: int, nt_change: str, gene: str, strand: str):
    # First try to use the provided position and nt_change if simple substitution 
    if position > 0 and ">" in nt_change:
        if not overlaps_gene(position, gene):
            # position doesn't even overlap the gene, must be GRCh37
            liftover = pygautil.liftover(chr_name.replace("chr", ""), position, position, "", "GRCh_38")
            if len(liftover) > 1:
                position = liftover[1]

        nc_chr = chromosome_to_nc(chr_name)
        if strand == "-":
            ref, alt = nt_change.split(">")
            refc = reverse_complement(ref)
            altc = reverse_complement(alt)
            nt_change = f"{refc}>{altc}"
        gdot = f"{nc_chr}:g.{position}{nt_change}"

        coordinates = pygautil.convertPdot(gdot)
        if coordinates_mismatch_reference(coordinates):
            coordinates = liftover_hgvs_coordinates(gdot, chr_name, coordinates)
        
        if converted_coordinates_are_valid(coordinates):
            return coordinates

    # if that failed, try to convert the provided hgvs
    if '.' not in hgvs and len(hgvs) > 30:
        return []
    
    coordinates = pygautil.convertPdot(hgvs)
    if converted_coordinates_are_valid(coordinates):
        return coordinates
    
    hgvs = fix_hgvs(hgvs, chr_name, strand)
    coordinates = pygautil.convertPdot(hgvs)

    if coordinates_mismatch_reference(coordinates):
        # must be in GRCh37 so liftover to GRCh38
        coordinates = liftover_hgvs_coordinates(hgvs, chr_name, coordinates)


    if converted_coordinates_are_valid(coordinates):
        return coordinates

    return []

def get_reference_sequence_bases(chr, start, stop):
    global REF_READER_38
    ref_track = REF_READER_38

    if chr.startswith("chr"):
        chr = chr.replace("chr", "")

    ref_seq = ''
    ref_track.readRegion(chr, start, stop, ['Data'])
    ref_seq = ""
    while ref_track.hasNext():
        values = ref_track.next()
        ref_seq += values[0]
    return ref_seq

def get_vcf_output_line(variant):
    chrom = variant['chr']
    if chrom.startswith("chr"):
        chrom = chrom.replace("chr", "")

    values = [
        chrom,
        str(variant['start']+1),
        variant['id'],
        variant['ref'],
        variant['alt'],
        '.',
        'PASS',
        'GeneName='+variant["GeneName"],
    ]

    return "\t".join([str(v) for v in values])

def convert_to_vcf_line(variant):
    alts = [variant["alt"]]
    isIns = variant["ref"] == "-" or variant["ref"] == ""
    isDel = False
    if variant["alt"] == "-" or variant["alt"] == "":
        isDel = True

    # construct vcf variant and write to vcf
    vcfVariant = dict(variant)
    if isDel:
        delVariant = dict(variant)
        delVariant['start'] -= 1
        refBase = get_reference_sequence_bases(delVariant['chr'], delVariant['start'], delVariant['start']+1)
        delVariant['ref'] = refBase + vcfVariant['ref']
        delVariant['alt'] = refBase
        return get_vcf_output_line(delVariant)
    if isIns:
        vcfVariant['start'] -= 1
        refBase = get_reference_sequence_bases(vcfVariant['chr'], vcfVariant['start'], vcfVariant['start']+1)
        alts = [refBase + a for a in alts]
        vcfVariant['ref'] = refBase

    if len(alts) > 0:
        vcfVariant['alt'] = ",".join(alts)
        return get_vcf_output_line(vcfVariant)

def convert_csv_to_vcf(csv, vcf):
    global VCF_HEADER
    vcf_file = open(vcf, "w")
    vcf_file.write(VCF_HEADER + "\n")
    fields = []
    variant_count = 0
    variants_to_write = []
    dropped_variants = []
    
    print("Loading csv...")
    for line in open(csv, "r"):
        if line.startswith("#") or not line.strip():
            continue

        values = [v.strip() for v in line.split("\t")]
        if len(fields) == 0:
            fields = values
            continue
            
        csv_record = {fields[i]:v for i,v in enumerate(values)}

        if "chr" not in csv_record:
            print("Invalid csv line: ")
            print(line)
            continue

        hgvs = csv_record["varID"]
        chr_name = csv_record["chr"]
        position = int(csv_record["gNomen"])
        nt_change = csv_record["ntChange"]
        strand = csv_record["strand"]
        gene = csv_record["gene"]
        converted = convert_hgvs(hgvs, chr_name, position, nt_change, gene, strand)
        if len(converted) < 1:
            print(f"Failed to convert variant: {hgvs}")
            dropped_variants.append(converted)
            continue

        start = converted[0][1]
        ref = converted[0][3].split("/")[0]
        alt = converted[0][3].split("/")[1]
        variant = {
            "chr": chr_name,
            "start": start,
            "ref": ref,
            "alt": alt,
            "id": csv_record["varID"],
            "GeneName": csv_record["gene"]
        }
        vcf_line = convert_to_vcf_line(variant)
        variants_to_write.append(vcf_line)

        variant_count += 1
        if variant_count % 1000 == 0:
            print(f"    Loaded {variant_count} variants...")

    print("Sorting variants...")
    variants_to_write = sort_variants_by_genome(variants_to_write)
    
    print("Converting to vcf...")
    variant_count = 0
    for vcf_line in variants_to_write:
        vcf_file.write(vcf_line + "\n")

        variant_count += 1
        if variant_count % 1000 == 0:
            print(f"    Wrote {variant_count} variants...")
    
    print(f"Conversion complete. Wrote {variant_count} variants. Dropped {len(dropped_variants)}.")

def test():
    result = convert_hgvs("NM_024426.4:c.996G>A", "11", 32438041, "G>A", "WT1", "-")
    print(result)

    """
    result = convert_hgvs("NM_020865:g.154018887G>C", "3", "-")
    print(result)

    result = convert_hgvs("NM_000132:g.154130395C>T", "X", "-")
    print(result)

    result = fix_hgvs("NM_000132:g.154855744_154855744:delinsGG", "X", "-")
    print(result)

    result = convert_hgvs("NM_183061:g.112278365:T>C", "3", "-")
    print(result)

    result = convert_hgvs("NM_000059.3:c.470_475del5", "13", "+")
    print(result)

    result = convert_hgvs("NM_007294.3:c.167_168delinGA", "17", "-")
    print(result)

    result = convert_hgvs("NM_000169:c. 639+919G>A", "X", "-")
    print(result)
    
    result = convert_hgvs("NM_016124:c.336-2del", "1", "+")
    print(result)
    """

def list_observations(csv):
    fields = []
    observations = set()
    interpretations = set()
    print("Loading csv...")
    variant_count = 0
    for line in open(csv, "r"):
        if line.startswith("#") or not line.strip():
            continue

        values = [v.strip() for v in line.split("\t")]
        if len(fields) == 0:
            fields = values
            continue
            
        csv_record = {fields[i]:v for i,v in enumerate(values)}
        observations.add(csv_record["observation"])
        interpretations.add(csv_record["Interpretation"])

        variant_count += 1
        if variant_count % 1000 == 0:
            print(f"    Wrote {variant_count} variants...")

    observations = sorted(list(observations))
    #print("Observations")
    #for observation in observations:
    #    print(observation)

    print("Interpretations")
    interpretations = sorted(list(interpretations))
    for interpretation in interpretations:
        print(interpretation)

def example():
    convert_csv_to_vcf("example.csv", "example.vcf")

def main():
    convert_csv_to_vcf("HUMU-43-2308-s003.csv", "HUMU-43-2308-s003-converted.vcf")

main()
