# Data formatting

At this stage, we assume that the alignment and the tree files are ready and that all the needed quality control checks have been carried out to make sure the best practice to infer both the sequence alignment and the corresponding phylogeny has been followed. Before proceeding with timetree inference, however, we will make sure that:

1. The alignment file is in PHYLIP format and easy to read (i.e., we recommend having one sequence per line in the alignment file).
2. The tree file is in Newick format.

In this tutorial, we will use [an example dataset](00_raw_data) to see how we can convert an alignment from FASTA format into PHYLIP format and from a tree in NEXUS format to a tree in Newick format.

## Alignment file

Inside the directory where this `README.md` file is, you can find a subdirectory called [`00_raw_data`](00_raw_data) with the example alignment and tree files aforementioned.

If you open [the alignment file](00_raw_data/raw_aln.fa), you will see that each aligned sequence is not in a unique line, which sometimes makes it more difficult to parse the file. The first thing that we will do now is to run [an in-house PERL script called `one_line_fasta.pl`](scripts/one_line_fasta.pl) to convert the `raw_aln.fa` FASTA file in a FASTA file in which all sequences are written in one line:

```sh
# Run the next commands from the 
# `00_data` directory
cd 00_raw_data
name=`ls *fa`
printf "Converting "$name" into a one-line FASTA file\n"
../scripts/one_line_fasta.pl $name
onefa=$( echo $name | sed 's/\.fa/\_one\_line\.fa/' )
namefa=$( echo $name | sed 's/\.fa//' )
mv $onefa $namefa.fasta
```

After running the code snippet above, you will see that a new FASTA file called `raw_aln.fasta` has been generate in [the `00_raw_data` directory](00_raw_data) with the format we wanted: one sequence per line. Now, we just need to run [another in-house PERL script called `FASTAtoPHYL.pl`](scripts/FASTAtoPHYL.pl), which will convert this newly generated alignment file from FASTA into PHYLIP format:

```sh
# You should still be inside `00_data/00_raw_data`
# If not, please move to this directory and run the
# following commands
aln_name=`ls *fasta`
a_noext=$( echo $aln_name | sed 's/\.fasta//' )
num=$( grep '>' $aln_name | wc -l )
len=$( sed -n '2,2p' $aln_name | sed 's/\r//' | sed 's/\n//' | wc --m )
perl ../scripts/FASTAtoPHYL.pl $aln_name $num $len 
# Create a directory for input data for `MCMCtree`
mkdir ../01_inp_data
mv $a_noext.phy ../01_inp_data
```

You will now see a new directory called `01_inp_data` inside the `00_data`. If you navigate to this newly created `01_inp_data` directory, you will find the alignment in PHYLIP format (i.e., the input file we need!). You will also find a log file called `log_lenseq.txt` inside [the `00_raw_data` directory](00_raw_data) where you can read how many taxa were parsed and the length of the sequence.

The alignment is now in the correct format, so we can start to parse the tree file!

## Tree file

If the calibrated input tree is in NEXUS format, as in [this example dataset](00_raw_data/tree_ML.nexus), we can use bash scripting to easily convert this file into Newick format. To process the example tree file, please run the following commands:

```sh
# You should still be inside `00_data/00_raw_data`
# If not, please move to this directory and run the
# following commands
tname=`ls *nexus`
t_noext=$( echo $tname | sed 's/\.nexus//' )
# Create a tree file in PHYLIP format by 
# getting the number of taxa from the FASTA
# alignment
aln_name=`ls *fasta`
num=$( grep '>' $aln_name | wc -l )
printf $num" 1\n" > ../01_inp_data/$t_noext"_calib.tree"
# Extract the line that has the tree in Newick format
# from the NEXUS TREE, remove characters before the tree,
# remove branch lengths, and remove weird characters
grep 'tree ' $tname | sed 's/..*\[\&R\]\ //' | sed 's/\[\&label\=//g' | sed 's/\]//g' | sed 's/\:[0-9]\.[0-9]*//g' | sed 's/\:[0-9]//g' | sed 's/\:[0-9]*\.[0-9]*e-[0-9]*//g' | sed 's/\:[0-9]e-[0-9]*//g' | sed 's/E-[0-9]*//g' >> ../01_inp_data/$t_noext"_calib.tree"
# Now, use the calibrated tree to generate the uncalibrated
# tree, which will be used to run `BASEML`. In this case, only
# ST calibrations are available. Otherwise, the `sed` commands
# should be modified to include the rest of the calibrations
# that need to be deleted
cp ../01_inp_data/$t_noext"_calib.tree" ../01_inp_data/$t_noext"_uncalib.tree"
sed -i 's/'\"'//g' ../01_inp_data/$t_noext"_uncalib.tree"
sed -i 's/ST([0-9]*\.[0-9]*,[0-9]*\.[0-9]*,[0-9]*\.[0-9]*,[0-9]*\.[0-9]*)//g' ../01_inp_data/$t_noext"_uncalib.tree"
sed -i 's/ST([0-9]*\.[0-9]*,[0-9]*\.[0-9]*,\-[0-9]*\.[0-9]*,[0-9]*\.[0-9]*)//g' ../01_inp_data/$t_noext"_uncalib.tree"
sed -i 's/B([0-9]*\.[0-9]*,[0-9]*\.[0-9]*)//g' ../01_inp_data/$t_noext"_uncalib.tree"
# If we want to have single quotation marks for the calibrations
# in the calibrated tree instead of double quotation marks, 
# we should do the following
sed -i 's/'\"'/'\''/g' ../01_inp_data/$t_noext"_calib.tree" 
```

> **NOTE 1**: If you decide to later adapt this tutorial to analyse a different dataset, please check the `sed` commands used in the code snippet above to generate the uncalibrated trees. Perhaps, there is a branch length or a calibration that has not been yet removed in your tree file/s by any of the regex used above as you may have included different characters (e.g., the `sed` command has not targeted `-`, other type of calibrations such as low bounds, etc.). Before moving onto the next step, please open the tree file with your preferred text editor such as `FigTree` or `TreeViewer` (or your preferred graphical software to visualise the edited phylogeny) and make sure that (i) there are no issues and that (ii) the final tree has no branch lengths, calibrations, or other labels. In other words, you should make sure that the `*_uncalib.tree` file **only** includes the tree topology in Newick format. The topology + calibrations (no branch lengths!) should be available **only** in the `*_calib.tree` file (i.e., file with the calibrated tree).
> **NOTE 2**: We have covered a specific formatting example (i.e., from NEXUS to Newick) as it tends to be the most common for tree files. Nevertheless, other format conversions may be required. Please let us know if you would like to see other data conversions as we are working on other tutorials that could complement this one.

---

Now, you are ready to move onto the next step: [calculating the Hessian and the gradient!](../01_analyses/01_Hessian)
