# LB1 final project: Building a Profile Hidden Markov Model for a protein domain (Kunitz-type protease inhibitor example)

This repository contains the instructions for building a Hidden Markov Model profile for a protein domain, as well as an example model generated for the final project of the Laboratory of Bioinformatics course of the Master's Degree in Bioinformatics of the Università di Bologna, course 2022-2023.

The project consists of building a profile-HMM for the Kunitz domain, starting from available protein structures from RCSB PDB. Once the model is built, it will be used to annotate domains in UniProtKB/Swiss-Prot proteins.

For more details, see also the written report (`project_lb1b-Giacomo_Orsini.pdf`).

Written by Giacomo Orsini.

**Table of Contents**
-  [0. Software and databases](#0-software-and-databases)
- [1. Structure selection and fetching](#1-structure-selection-and-fetching)
    - [1.1. PDB Search](#11-pdb-search)
    - [1.2. Protein Alignment](#12-protein-alignment)
    - [1.3. (Optional) Trimming the MSA](#13-optional-trimming-the-msa)
- [2. HMM Generation](#2-hmm-generation)
    - [2.1. Train a HMM profile](#21-train-a-hmm-profile)
    - [2.2. Consistency test](#22-consistency-test)
- [3. Method Testing](#3-method-testing)
    - [3.1. Testing set](#31-testing-set)
    - [3.2. Cross-validation sets](#32-cross-validation-sets)
    - [3.3. Model evaluation](#33-model-evaluation)
    - [3.4. Conclusion](#34-conclusion)

## 0. Software and databases
For this project, you have to download the `HMMER` package ([website](http://hmmer.org/)). The example project was conducted using `HMMER v3.3.2` (Nov 2020) version. 

- All UniProt searches were performed using Release 2023_02.
- All PDB searches were conducted on 27/05/2023.
- All PDBeFold submissions were issued under PDBe Fold v2.59.

## 1. Structure selection and fetching
The first step to build a profile-HMM is to retrieve the available domain structures on [PDB database](https://www.rcsb.org/rcsb.org/) and generate a multiple sequence alignmnet with [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/). As said, we will work with the Kunitz domain as an example. 

### 1.1. PDB Search
First, access the [Advanced search](https://www.rcsb.org/search/advanced) section of the PDB website, and write the desired query. In my case, the parameters have been:
- `Pfam identifier OR SCOP2 lineage identifier OR CATH lineage identifier`: `PF00014`, `4003337` or `4.10.410.10` (select those structurally-resolved PDB structures that have been annotated as containing a Kunitz domain according to the main classification databases).
- `AND Refinement Resolution`: `<= 3`.
- `AND Polymer Entity Sequence Length`: `51 - 76` residues (size range of the Kunitz domains. This subquery also ensures the exclusion of protein chains that have more than one domain).
- `AND Polymer Entity Mutation Count`: `0` (wild-type versions of the protein, no mutants). Alternatively, as sometimes mutations are not reported, remove entries that have "mutant" in the title (`AND Structure Title has NOT any of words: mutant - variant`).

This query returned >100 structures. However, different PDB files may refer to the same protein. Therefore, you have to group  those structures that belong to the same protein. To do this, select in the advanced query the `Return Polymer Entities` that are `grouped by Sequence Identity 95%` and `displaying as Representatives`. 

After retrieving these entities, generate a Custom Tabular Report. You can choose different sets of information, although, for simplicity, these are the ones that you will need:
- `PDB ID`.
- `Auth Asym ID`: chain ID given by the author.
- `Accession code(s)`: UniProtIDs related to each PDB entry.

Download the tabular report in CSV format (`pdb_report.csv`). 

### 1.2. Protein Alignment
To conduct the MSA on PDBeFold, create a text file (`pdb_codes.txt`) that contains the PDB identifiers in the PDBeFold format (e.g 3abc:D, entry + chain containing the domain):
```
grep -v 'Identifier\|Entity\|,,' pdb_report.csv | cut -d ',' -f 2,3 | tr -d \" |tr "," ":" >pdb_codes.txt
```
- `grep -v`: select the lines in the input file that do <u>not</u> contain the subsequent terms. `-i` enables case insensitive matching.
- `cut`: cut out selected portions of each line of the input file, with the desired separator character `-d`, and output the selected fields `-f`.
- `tr`: substitute (translate) first character by second character. `-d` option deletes that character from the input file.
- `>`: redirect final output to the output file.

Then, access the PDBeFold website:
1. Click on `Launch PDBeFold`
2. Select, inside `Multiple Submission Form`, `List of PDB codes` as the source
3. Upload the file (`pdb_codes.txt`) 
4. Once all entries are displayed, click on `Update List` (with Jmol as the viewer)
5. Submit the query

Once the request has been completed (it might take a few minutes), download the MSA in FASTA format (`msa.fasta`).

&emsp;<u>Note</u>: To check that the MSA is correct, download the Multiple Alignment Results (`msa.dat`) or look at the results on the web page: all sequences must have an RMSD <2. If this is not the case, remove those entries (`grep -v`) from the PDB codes file (`pdb_codes_ref.txt`), re-run the MSA, and download the new file (`msa_ref.fasta`).

### 1.3. (Optional) Trimming the MSA
After retrieving an MSA, we can trim the alignment and check for eventual redundant entries. Poorly aligned or highly gapped regions can influence the resulting HMM profile, so trim these areas with the following commands:
```
grep . msa.fasta |awk '{if (substr($1,1,1)==">") {printf "\n%s ",$1} else {printf "%s",$0}}'|less
```
First, identify `pos1` and `pos2`, the initial and final position of the residues to include. Then proceed with the trimming:
```
grep . msa.fasta |awk '{if (substr($1,1,1)==">") {printf "\n%s ",$1} else {printf "%s",$0}}'|awk '{print $1; print toupper(substr($2,pos1,pos2))}'>msa_trimm.fasta
```
- `grep .`: remove empty lines.
- `awk`: it's a language to modify text files. Between "'{}'" are written the instructions used to clean the fasta file (for more details follow [this link](https://www.geeksforgeeks.org/awk-command-unixlinux-examples/)). We want to remove everything except the PDB ID and the chain, followed by the sequence, in uppercase, from `pos1` to `pos2`.

After trimming, check that there are no duplicated sequences within the MSA. If there are, remove one of the sequences (duplicate sequences can introduce redundancy and potentially bias the model training process). To do so, upload the trimmed file (`msa_trimm.fasta`) to the [Align](https://www.uniprot.org/aligng) tool at the UniProt website, and delete manually (using vim, vi, nano or any other text editor) the FASTA entry and sequence of one of the duplicated sequences.

## 2. HMM Generation
The second step consists in the generation of the HMM profile from the selected PDB entities. It's performance will be first verifyed on the training set.

### 2.1. Train a HMM profile

Once the (trimmed) MSA is ready, build the HMM profile (where `kunitz_mod.hmm` will be the output HMM profile and `msa_trimm.fasta` is the input MSA):
```
hmmbuild kunitz_mod.hmm msa_trimm.fasta
```

### 2.2. Consistency test
To verify that the trained HMM is able to recognize the proteins in the dataset (consistency test), perform a `hmmsearch` with the sequences of the proteins used to generate the model. 

The FASTA format sequences can be retrieved with the [Retrieve/ID Mapping](https://www.uniprot.org/id-mapping) tool at the UniProt website. The IDs to retrieve can be loaded from a text file. 

Therefore, create a text file (`training_set.idlist`) that contains the (unique) UniProt identifiers of the model proteins:
```
grep -v 'Identifier\|Entity\|,,' pdb_report.csv | cut -d ',' -f 4 | tr -d \" |sort |uniq >training_set.idlist
```
- `sort`: sort input file in ascending order
- `uniq`: filter out repeated lines in the input file (file must be sorted)

&emsp;<u>!Note!</u>: Remember to remove (if any) entries that scored RMSD <2 when performing the MSA, and the eventual duplicates removed in the trimming procedure (although the duplicates usually refer to the same Uniprot id). 

Map the IDs in the UniProt website to retrieve their sequences and download them, compressed, as FASTA (canonical) format (`training_set.fasta`). Finally, run `hmmsearch` (Usage: `hmmsearch [options] <hmmfile> <seqdb>`):
```
hmmsearch --max --noali -o training_set.out kunitz_mod.hmm training_set.fasta
```
- `--max`: turn off all the heuristics for cutting off distantly related proteins.
- `--noali`: don't include in the output the visual representation of the alignment between the model and each entry.
- `-o`: name of the output file.


Check the output file on the training set, evaluating the E-value, alignment quality, consistency, and coverage. These should all be very high, as we are running the model over the same set of sequences used to generate it. If some entries score much less than the others or they have a high E-value, you might need, for example, to exclude them from the MSA or retrieve a different set of structures.

## 3. Method Testing
The third step is to retrieve a suitable test set, consisting of known proteins that contain (positive set) or not (negative set) the Kunitz domain, search it against the trained model, and compute the scoring indexes to evaluate the HMM profile on cross-validation and testing sets.

### 3.1. Testing set
All FASTA sequences for the test set have been retrieved with the Advanced search option at the [UniProt website](https://www.uniprot.org/), restricting the search to UniProtKB/Swiss-Prot entries (reviewed proteins). We will first download two subsets. To retrieve those proteins that contain a Kunitz domain, add the subquery Pfam identifier `PF00014` in the Advanced Search. This search will contain the proteins that will make up the **positive set** (`/UniProt_Swiss-Prot_data/db_kunitz.fasta`). All the remaining proteins in UniProtKB/Swiss-Prot will be used as the **negative set** (`/UniProt_Swiss-Prot_data/db_nonkunitz.fasta`).

Once both sets are downloaded, create a file that contains the UniProtIDs for each set:
```
grep '^>' db_kunitz.fasta |cut -d '|' -f 2 >db_kunitz.idlist
grep '^>' db_nonkunitz.fasta |cut -d '|' -f 2 >db_nonkunitz.idlist
```

For a fair test on the model, the positive test set should exclude the training data. 
To remove those sequences that are part of the positive set and that have been used in the training test, create a file (`db_kunitz_toremove.idlist`) that doesn't include the UniProtIDs of the proteins used to build the model (`training_set.idlist`):
```
comm -13 <(sort training_set.idlist) <(sort /UniProt_Swiss-Prot_data/db_kunitz.idlist) >db_kunitz_toremove.idlist
```
- `comm`: compare two files and show the lines that are unique or in common. The first column shows lines unique to the first file, the second column shows the lines unique to the second file and the third column shows lines common in both. `-13` returns only the unique lines found in the second file.

Once the UniprotIDs of the proteins that have <u>not</u> been used in the training set are obtained, go to the [Retrieve/ID Mapping](https://www.uniprot.org/id-mapping) tool at the UniProt website to retrieve and download the FASTA file (`db_rmkunitz.fasta`).

Next, create a FASTA file that contains all of the sequences (`total_set.fasta`) of the positive and negative sets, which will be our final testing set:
```
cat db_rmkunitz.fasta ../UniProt_Swiss-Prot_data/uniprot_sp_nonkunitz.fasta >total_set.fasta
```

Run `hmmsearch`, with the file that was just generated (`total_set.fasta`) as a sequence database:
```
hmmsearch --max --noali -o total_set.search kunitz_mod.hmm total_set.fasta 
```
- `--noali` (optional)

Refine the output to just include the hits that contain the E-values, UniProtIDs, and other information (`final_pos` corresponds to the number of the last row containing this information):
```
head -n final_pos val_set.search |tail -n +18|grep -v inclusion >test_set.out
```
- `head -n [x]`: display the first x lines of the input file.
- `tail -n +[x]`: display the lines starting from line x of the input file, till the end.

&emsp;<u>Note</u>: to know the `final_pos`, you'll have to look at the file manually. To do that, I suggest using the vi text editor, as it shows the number of the row you are in. 

Now that the hits are refined, generate a file for each two sub-validating set () that will be used later for the evaluation of the performance in the format `UniProtID E-value class`, where `E-value` is the E-value for the highest scoring hit (domain) and the `class` is the actual ("real") annotation, obtained from UniProt: 1 for proteins that contain a Kunitz domain and 0 for proteins that don't contain a Kunitz domain.

```
awk '{split($9,i,"|"); print i[2],$4}' val_set.out |grep -f pos_set.idlist |awk '{print $1, $2, 1}' >test_set.classp
awk '{split($9,i,"|"); print i[2],$4}' val_set.out |grep -v -f pos_set.idlist |awk '{print $1, $2, 0}' >test_set.classn
```
- `grep -f`: find a pattern (could be another file) inside a file. The output will be the lines containg a word from the pattern. The option `-v` outputs, at the contrary, the lones that have no pattern.    

Finally, add all the proteins retrieved from UniProt as non-Kunitz to the class file. To do so, first, get the UniProtIDs that have been a hit against the model:
```
cut -d ' ' -f 1 test_set.classn >neg_set.hits
```

And lastly, retrieve those IDs that did not appear in the `hmmsearch` output, as the alignment was not computable, and add them to the class file (with a high E-value, eg `100`, to evidence that they did not match the model):
```
comm -13 <(sort neg_set.hits) <(sort ../UniProt_Swiss-Prot_data/db_nonkunitz.idlist) | awk '{print $0,100,0}' >>test_set.classn
```

### 3.2. Cross-validation sets
The evaluation procedure consists of a 2-fold cross-validation test: the positive and negative sets are split into two subsets (`set1.txt` and `set2.txt`), then the optimal classification (E-value) threshold is found for each subset (cross-validation set) and the performance are tested inverting the best thresholds on the two subsets(testing set). The average of the two optimal thresholds will also be tested on a total set containing set1 and set2 (`val_set.txt`). For this procedure, you will need a python script (`confusion_matrix.py`)

Before dividing the cross validation sets, shuffle the content of the files, and save it in a new file (`val_set.randomp` and `val_set.randomn`):

```
sort -R val_set.classp > val_set.randomp
sort -R val_set.classn > val_set.randomn
```
- `sort -R`: sort the file randomly

To divide the sets into two different subsets, count the number of entries/lines (`wc -l`) there are in each of the `val_set.randomp` and `val_set.randomn` files, add it to a variable (`n_p`, number of positives, `n_n`, number of negatives) and half these variables to obtain the number of entries in half of each set (`n_p_h`, half number of positives, `n_n_h`, half number of negatives):
```
n_p=$(wc -l val_set.randomp |awk '{print $1}')
n_p_h=$(($n_p/2))
n_n=$(wc -l val_set.randomn |awk '{print $1}')
n_n_h=$(($n_n/2))
```

Finally, select the first half of each of the `val_set.randomp` and `val_set.randomn` files, and add it to the first subset (`set1.txt`) and the second half of each file and add it to the second subset (`set2.txt`):
```
head -n $n_n_h val_set.randomn >set1.txt
head -n $n_p_h val_set.randomp >>set1.txt
tail -n +$(($n_n_h+1)) val_set.randomn >set2.txt
tail -n +$(($n_p_h+1)) val_set.randomp >>set2.txt
```
Finally, create the text file that will contain the full validation set (`val_set.txt`), to compute the performance on the whole dataset: 
```
cat val_set.randomp val_set.randomn > val_set.txt
```
### 3.3. Model evaluation
Compute the scoring indexes for evaluating the profile HMM on the validation sets. This is done in 5 sub-steps:
1. Calculate the optimised E-value on the first dataset (`set1.txt` = cross-validation set)
2. Use this E-value on the second dataset and evaluate its performances (`set2.txt` = testing set)
3. Calculate the optimised E-value on the second dataset (`set2.txt` = cross-validation set)
4. Use this E-value on the first dataset and evaluate its performances (`set1.txt` = testing set)
5. Compute the average E-value and get the performances on the whole dataset (`val_set.txt` = testing set)

To find the optimised E-values and their performances - Matthews correlation coefficient (MCC), accuracy (ACC), and confusion matrix (CM) - run the python script `confusion_matrix.py`, which takes as first argument the file to be evaluated and as second argument the E-value threshold to use (where `eval2` and `eval1` are the optimised E-values obtained from the cross-validation sets 1 and 2, respectively). The output is saved in a file (`.cvs` stores the results of the cross-validation sets and `.tv` of the testing set):
```
for i in `seq 1 20`; do python ../confusion_matrix.py set1.txt 1e-$i; done >set1.cvs
for i in `seq 1 20`; do python ../confusion_matrix.py set2.txt 1e-$i; done >set2.cvs
python ../confusion_matrix.py set1.txt eval2 > set1.tv
python ../confusion_matrix.py set2.txt eval1 > set2.tv
```

Finally, run the Python script with the average E-value (`evalavg`), and save the output on a file. This will be the E-value threshold at which your model works best:
```
python ../confusion_matrix.py val_set.txt evalavg >val_set.tv
```
### 3.4. Conclusion
You now have a complete, tested profile-HMM to annotate Kunitz-type domains. You can perform the same procedures for any other protein domain. 
 
