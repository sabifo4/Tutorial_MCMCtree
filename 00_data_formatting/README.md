# Data formatting

For this example, we assume that we have inferred both our molecular alignment and the corresponding phylogeny, and hence have the alignment and the tree files ready. In addition, we also assume that all the needed quality control checks have been carried out to make sure the best practice to infer both the sequence alignment and the corresponding phylogeny has been followed. Before proceeding with timetree inference, however, we will make sure that:

1. The alignment file is in PHYLIP format and easy to read (i.e., ideally one sequence per line).
2. The tree file is in Newick format.

## Alignment file

Inside the directory where this `README.md` file is (i.e., [`00_data_formatting`](README.md)), you can find a subdirectory called [`00_raw_data/alignment`](00_raw_data/alignment) with the sequence alignment that we will use in this tutorial.

If you open [the alignment file](00_raw_data/alignment/raw_aln.fa), you will see that each aligned sequence is not in a unique line, which sometimes makes it more difficult to parse the file. First, we will run [an in-house PERL script called `one_line_fasta.pl`](../src/one_line_fasta.pl) to convert the `raw_aln.fa` FASTA file in a "one-line" FASTA file:

```sh
# Run the next commands from the 
# `00_raw_data/alignment` directory
name=`ls *fa`
printf "Converting "$name" into a one-line FASTA file\n"
chmod 775 ../../../src/*
perl ../../../src/one_line_fasta.pl $name
onefa=$( echo $name | sed 's/\.fa/\_one\_line\.fa/' )
namefa=$( echo $name | sed 's/\.fa//' )
mv $onefa $namefa.fasta
```

After running the code snippet above, you will see that a new FASTA file called `raw_aln.fasta` has just been generate in [the `alignment` directory](00_raw_data/alignment/), and the format is exactly what we wanted: one sequence per line. Now, we just need to run [another in-house PERL script called `FASTAtoPHYL.pl`](../../src/FASTAtoPHYL.pl), which will convert the sequence file into a PHYLIP formatted file:

```sh
# You should still be inside `00_raw_data/alignment`
# If not, please move to this directory and run the
# following commands
aln_name=`ls *fasta`
a_noext=$( echo $aln_name | sed 's/\.fasta//' )
num=$( grep '>' $aln_name | wc -l )
len=$( sed -n '2,2p' $aln_name | sed 's/\r//' | sed 's/\n//' | wc -L )
perl ../../../src/FASTAtoPHYL.pl $aln_name $num $len 
# Create a directory for input data for `MCMCtree`
mkdir ../../01_inp_data
mv $a_noext.phy ../../01_inp_data/example_aln.phy
```

You will see a new directory called `01_inp_data` inside [`00_data_formatting`](README.md) directory. If you navigate to this newly created directory, you will find the alignment in PHYLIP format (i.e., the input file we need!). You will also find a log file called `log_lenseq.txt` inside [the `alignment` directory](00_raw_data/alignment/) where you can read how many taxa were parsed and the length of the sequences.

The alignment is now in the correct format, so we can parse the tree file!

## Tree file

After inferring your phylogeny, you will obtain a tree topology in Newick format with branch lengths such as what you can see in our input file: [`tree_ML.tree`](00_raw_data/trees/tree_ML.tree).

Your next goal is to calibrate this tree topology. First, however, we will need to get rid of the branch lengths!

```sh
# Run from `00_raw_data/trees`
cp tree_ML.tree ../../01_inp_data/tree_example_uncalib.tree
sed -i 's/:[0-9]*\.[0-9]*//g' ../../01_inp_data/tree_example_uncalib.tree
# NOTE: This regular expresion will work with that example
# file. You may have to use more complex regular expressions
# if you have `E-` or even bootstrap values that you need
# to get rid of!
#
# Add header
sed -i '1s/^/4 1\n/' ../../01_inp_data/tree_example_uncalib.tree
```

Now, we can calibrate this tree topology following the calibration files we have designed. In a nutshell, the tasks carried out by our [R in-house script](scripts/Include_calibrations.R) are the following:

* Read the input file with calibration information (i.e., [`calibrations.csv`](00_raw_data/calibs/calibrations.csv)), and the corresponding input tree file with the tree topology that shall be fixed when running PAML programs (i.e., the uncalibrated tree file we have just generated and saved in `01_inp_data`).
* Find the corresponding nodes for the MRCAs for each calibrated node so that the tag that will be used to identify the calibration for such node (i.e., first column in the calibration files) is incorporated as a node label.
* Identify nodes that are assigned more than one calibration (i.e., potential duplicates in the calibration file) and decide what to do with them. Then, repeat the previous step with the filtered calibration file, if generated, and output the resulting tree topology with node labels.
* Read the previously output file as a character vector (not as a `phylo` object!) to replace the node labels with the corresponding node calibrations in `MCMCtree` notation. Columns 4-7 (or, alternatively, 8 if you already know the `MCMCtree` notation) in the calibration files are used for this purpose.

Once you have run our [R in-house script `Include_calibrations.R`](scripts/Include_calibrations.R), you will see that the calibrated tree files with nodes labelled following `MCMCtree` notation are output in the `01_inp_data` directory (i.e., check file which name ends with `calib_MCMCtree.tree`). In addition, a file ending with `fordisplay_calib_MCMCtree.tree` will have been output in directory [`00_raw_data/trees`](00_raw_data/trees/), which can be opened with graphical viewers such as `FigTree` or `TreeViewer` ([Bianchini and Sánchez-Baracaldo, 2024](https://onlinelibrary.wiley.com/doi/10.1002/ece3.10873)) to see which nodes have been calibrated and how.

Now, we can move onto the next step: [we can calculate the Hessian, the branch lengths, and the gradient with `CODEML`!](../01_PAML/00_CODEML/README.md)
