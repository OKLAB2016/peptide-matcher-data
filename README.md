# peptide-matcher-data

This repository contains the workflow used to generate fasta datasets compatible with [Peptide matcher](https://github.com/OKLAB2016/peptide-matcher/) or potentially other software capable of parsing the data. To download the datasets, go to [releases](https://github.com/OKLAB2016/peptide-matcher-data/releases).

## Workflow description

The workflow takes [alphafold](https://alphafold.ebi.ac.uk/) structures and the corresponding [uniprot](https://www.uniprot.org/) records and extracts structural information from the two to create a single dataset file that can be used for analyses of structural features in proteins with a single-residue resolution.

## Dependencies

The workflow is implemented in [snakemake](https://snakemake.readthedocs.io/) and depends on `dssp`, `python>=3.6`, `wget` and requires `pandas` and `biopython`.

## Format description

The format of the datasets is an extension of [fasta](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/):

```
>P0A8G3 Uronate isomerase confidence:49545c5f5e606062626263616262636363636363636363636262636263636362636262626162626263626262626262616261616161605e5c52505b5f6161626262626262626262616162605f5d60626262626262626363636363636363636362626262626363626263626262626262626262615e6062626263626363626263636362636262626262625f60626262626263636262626263626263636363636363636363636363636262626161615f5f606062626262636363636363626262626261605e60626262636362636363636363626262626162636363636363636363636363636363626262616262626363636262616262626363626262616263636363636362626262626162616262626363636363636363636363636363636363636262626062636363636363626262605f626160616060605f5f5e5e5d5e5f616262636262616262636363636363636363636262605e5f61626262636363636363636362636363636363636363626262625a58585a5e616263636363636363626263636363636363636363636363636363626262636262636363636362615f5e6062626262626363636363636363636363636363636363626262626262626362626363636363636363636363626262606062615e5b595949 secstruct:3-2T1-4T2-1S5H6I2T2-3E1-1S4-5H2T3-2S5H3T1-9H2T2-3G2T3S1-11H4G2T1S11H3T3-1S3-1S3T12H1T1S3G1S7H2T8E1-2T3-8H1-2T2-1S6E2-4H1-2T1S2T12H1T4-1S17H2T2-6E2S7-12H2T4-22H1T1-6E1-2E4-7H1-4S2-2E5-12H3T3-6E2S3G7H4G2-2T1S2T1S3E1-1S2-3G1-1S14H1S1-3G8-1S1-4T19H2T2S3-9H5I6H1T4- accessibility:6c473b250c18291f080e032701310a001a26000116350d00243d0f0e0301000000001f011b0e2a1400052f1d112404430b0008070600012b1304130c01051f0b010f4621420f29152e160117224204251c1a2705001e0a00001f01001e3b053f50240f020a150200040001010514090333033559152a000e35460e022f2e00092c1e00133e100031283d220201000705000030240c1f0226000000010f06090012212902351503122b04444a1b4b3b05310414050e0101000f000c2a00000e002329453001133b04033531000c4621054b0c34083b290510100019280000172408002a16002242151c0221020000000405181506100d0734182c31362504122701002716091d3c312b510b23383410050607000802100002140c00000909001633182603000000000001050122151b281f35440d1f453b1c142d2e0d0510021a2008342e4804203918001c1d00002716031b384441110203080000000104000c09280704251700002d05022e3808251f5b5521340e0901020100040805081405031314290100112201001c11021a41231a3a0520290002000404020d082e060101010101010209020104000000000404010027170b3b323f4110313f250c3631070026100016180106060e08021f3105002b05481556 transmembrane:470-
MTPFMTEDFLLDTEFARRLYHDYAKDQPIFDYHCHLPPQQIAEDYRFKNLYDIWLKGDHY
KWRAMRTNGVAERLCTGDASDREKFDAWAATVPHTIGNPLYHWTHLELRRPFGITGKLLS
PSTADEIWNECNELLAQDNFSARGIMQQMNVKMVGTTDDPIDSLEHHAEIAKDGSFTIKV
LPSWRPDKAFNIEQATFNDYMAKLGEVSDTDIRRFADLQTALTKRLDHFAAHGCKVSDHA
LDVVMFAEANEAELDSILARRLAGETLSEHEVAQFKTAVLVFLGAEYARRGWVQQYHIGA
LRNNNLRQFKLLGPDVGFDSINDRPMAEELSKLLSKQNEENLLPKTILYCLNPRDNEVLG
TMIGNFQGEGMPGKMQFGSGWWFNDQKDGMERQMTQLAQLGLLSRFVGMLTDSRSFLSYT
RHEYFRRILCQMIGRWVEAGEAPADINLLGEMVKNICFNNARDYFAIELN
```

In addition to the record identifier and description, the definition line contains additional data for each position of the protein:

| Field           | Description                        | Source    | Values                                            | Encoding               |
|:---------------:|:----------------------------------:|:---------:|:------------------------------------- -----------:|:----------------------:|
| `confidence`    | pLDDT confidence scores            | alphafold | Integers 1-100                                    | Hexadecimal numbers    |
| `secstruct`     | Secondary structure                | dssp      | dssp codes                                        | CIGAR-like compression |
| `accessibility` | Relative solvent accessibility     | dssp      | `ASA/max(ASA)` where max. values for abs. solvent accessibility are the theoretical maxima from [Tien et al 2013](https://doi.org/10.1371/journal.pone.0080635) | Hexadecimal numbers |
| `transmembrane` | Transmembrane regions              | uniprot   | T: transmembrane, S: signal peptide, -: otherwise | CIGAR-like compression |              |

Two encoding are used to compress the data to short strings:

* 2-digit hexadecimal numbers concatenated into a single string. It is obtained as `encoded = ''.join('%02x' % x for x in scores)` and can be decoded as: `scores = [ int(encoded[i:i+2], 16) for i in range(0, len(encoded), 2) ]`.
* The [CIGAR](http://samtools.github.io/hts-specs/SAMv1.pdf)-like compression. It is a simple format of the form `<int><char><int><char>...` where the integers signify the number of times the characters appear consecutively. The encoding is obtained as `encoded = ''.join(str(sum(1 for x in g)) + k for k, g in itertools.groupby(chars))` and can be decoded as `chars = list(itertools.chain(*[ [k] * int(n) for n, k in re.findall('(\d+)(.)', encoded) ]))`.

## How to run

Run as follows: `snakemake -c$CPUS --resources unisave=5 --config organism=$ORGANISM` where `$ORGANISM` is the organisms code in the form `UP000000589_10090_MOUSE` (uniprot proteome, taxid, short name) as appears on the alphafold's website.

The default configuration file is in `workflow/config.yml`, among the relevant values there is the alphafold version.
