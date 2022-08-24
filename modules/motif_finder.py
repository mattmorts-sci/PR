import pandas as pd
import seaborn as sns
from Bio import SeqIO
from modules.log import log
from datetime import datetime

def motif_finder(file, number, direction="n"):
    """
    Finds and counts motifs at either terminal end of the sequence
    Takes the number of residues in the motif and the terminal, either
    'n' or 'c'
    """

    seq_lst = []
    residues = []

    for seq_record in SeqIO.parse(f"{file}", "fasta"):
        sequence = str(seq_record.seq).upper()  # Puts the sequence into a variable
        seq_lst.append(sequence)  # Appends the sequence to seq_lst

    # Takes the specified number of amino acids and appends them to the list residues
    for sequence in seq_lst:
        if direction == "n":
            residues.append(sequence[:number])
        elif direction == "c":
            residues.append(sequence[-number:])
        else:
            print('Error: must specify "n" or "c"')
            break

    # Removes redundant amino acid combinations by creating a set then makes a dict of them
    # using the amino acid combinations as the keys and empty values of 0
    aa_set = set(residues)
    aa_dict = dict.fromkeys(aa_set, 0)

    # Goes back to the list of sequence and counts the number of times a particular amino acid
    # combination is found and updates the dict
    for sequence in seq_lst:
        if direction == "n":
            aa_triplet = sequence[:number]
            if aa_triplet in aa_dict:
                aa_dict[aa_triplet] += 1
        else:
            aa_triplet = sequence[-number:]
            if aa_triplet in aa_dict:
                aa_dict[aa_triplet] += 1
    # Creates a list of tuples using list comprehension
    tup_lst = [
        (k, v) for k, v in aa_dict.items()
    ]  

    # Sorts the tuple list by the key which is set to the second
    # element in the tuple
    tup_lst.sort(
        key=lambda x: x[1], reverse=True
    )  

    print(tup_lst[:10])

    main_aa = tup_lst[0][0]
    main_aa_count = tup_lst[0][1]

    count = 0
    with open(f"{file}", "r") as c:
        for line in c:
            if line.startswith(">"):
                count += 1

    a = main_aa_count / count

    print(
        f'The most common aa at {direction}-terminal is "{main_aa}" and its freq. is: '
    )
    print("{:.3f}".format(a))

    return aa_dict


def chart(file, number, direction, x_size=8, y_size=40):
    """Takes the dict from the motif_finder function and produces a chart of the most
    to least common motifs"""

    aa_dict = motif_finder(file, number, direction)

    # Converts the dict to a pandas dataframe with the columns amino_acid_combination and Count, resets index
    aa_df = pd.Series(aa_dict, name="Count")

    df2 = aa_df
    df2.index.name = "amino_acid_combination"
    df3 = df2.reset_index()

    # Visulise distribution of sequence lengths using seaborn
    sns.set(rc={"figure.figsize": (x_size, y_size)})

    sns.barplot(
        data=df3.sort_values("Count", ascending=False),
        y="amino_acid_combination",
        x="Count",
    )


def terminal_filter(file, data_source, terminus, *args):

    date = datetime.now().strftime("%y%m%d")

    codon = tuple([*args])

    temp_dict = {}

    for seq_record in SeqIO.parse(f"{file}", "fasta"):
        sequence = str(seq_record.seq).upper()  # Converts sequence to uppercase string
        temp_dict[
            seq_record.id
        ] = sequence  # Adds to dict with header as key and sequence as value

    with open(
        f"output/{date}_NR_{data_source}_{terminus}-terminalCodon.fasta", "w"
    ) as filt:
        for k, v in temp_dict.items():
            if terminus == "n":
                if v.startswith(codon):
                    filt.write(">" + k + "\n" + v + "\n")
            elif terminus == "c":
                if v.endswith(codon):
                    filt.write(">" + k + "\n" + v + "\n")

    log(
        f"output/{date}_NR_{data_source}_{terminus}-terminalCodon.fasta was created, filtering by {terminus}-terminal aa {codon} "
    )

    print(f"output/{date}_NR_{data_source}_{terminus}-terminalCodon.fasta")
    print(f"terminal aa included were: ")
    print(*codon, sep=", ")