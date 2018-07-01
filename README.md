RNAMake
=======
-----------
(c)  J.Yesselman, Stanford University, 2014-2018

**RNAMake** is a toolkit for designing and optimizing RNA 3D structure. It allows 
the alignment between RNA motifs. These motif are small modular pieces of RNA that are 
believed to fold independently, thus attaching them together with helix flanking both 
sides allows users of **RNAMake** to build large segments of RNA with a high success 
rate of forming the predicted structure _in vitro_.

Install
-------
```bash
git clone https://github.com/jyesselm/RNAMake.git
cd RNAMake
sudo pip install -r requirements.txt
```

Make sure to have a valid c++ compiler either g++ (> 4.6) or clang as well as python 2.7. In your .bashrc (.bash_profile if mac OSX) add:

```bash
# location to RNAMake directory used by RNAMake
export RNAMAKE=<RNAMake Path> 
# adding point python to RNAMake python code 
export PYTHONPATH=$PYTHONPATH:$RNAMAKE
```
[OPTIONAL] but useful

```bash
# location of RNAMake c++ executables 
export PATH=$PATH:$RNAMAKE/rnamake/lib/RNAMake/cmake/build/ 
# location of python executable scripts                         

note if you are using c shell or another non bash shell you will need to use the equivalent commands


Compile
------- 
requires `cmake` and `ninja`

cmake: https://github.com/Kitware/CMake <br>
ninja: https://github.com/ninja-build/ninja

make sure they are in your `$PATH`. Then run:

```bash
python compile.py 
```
executables are located in

```bash
$RNAMAKE/rnamake/lib/RNAMake/cmake/build/
```

Tests [OPTIONAL]
-----
To run unit tests for python:
```
$RNAMAKE/rnamake/unittests/run_unittests.sh
```

To run unit tests for C++ code:

```bash
cd $RNAMAKE/rnamake/lib/RNAMake/cmake/build/
python run_unittests.py
```


Applications
============

design_rna
-----------
Generates segments are RNA between two Watson-Crick basepairs can also perform sequence optimization for helical sequences. 

```
design_rna  [-pdb pdb_file.pdb ]
			[-start_bp "start_bp_name" ]
			[-end_bp "end_bp_name" ]
			[-mg motif_graph.mg ]
			[-out_file name_of_outfile.out ]
			[-score_file name_of_scorefile.scores ]
			[-designs num_of_designs ]
			[-seqs_per_design seq_per_designs ]
			[-search.accept_score accept_score]
			[-show_sections ]
			[-mc ]
			[-only_ideal]
			[-verbose ]
```
#### Required Arguments

Argument  | Description
------------- | -------------
-pdb		    | PDB file containing starting RNA, make sure that both the start basepair and end basepair are watson and crick base pairs that are end at RNA chains. Must supply both `-start_bp` and `-end_bp`. See figure below.
-start_bp			    | The Watson-Crick basepair to start building the RNA segment from. Example "A194-A252" the base pair between resiudes 192 and 252 both on chain A. 
-end_bp			 |	The Watson-Crick basepair to end the RNA segment. Same naming convention as `start_bp`. See examples below.
-mg			    | Supplies a motif graph file instead of a pdb, `start_bp` and `end_bp`; see documention for more info. 

#### Example of Basepair Ends that can be connected with `design_rna`

![basepair_end_examples](readme_resources/ggaa_tetraloop.png "Basepair End Example")

#Examples
<a href="#-">link</a> - (click to go to first anchor of this comment)

examples are located: /RNAMake/examples/cpp/design_rna


start.pdb:
![basepair_end_examples](readme_resources/ggaa_tetraloop.png "Basepair End Example")

<a id="-"></a>


Simplest use, generating one design

```
design_rna -pdb start.pdb -start_bp A222-A251 -end_bp A149-A154 -pdbs
> DESIGN RNA: loaded pdb from file: start.pdb
> DESIGN RNA: generated 1 design(s)! if you would like more please specify how many you would like with -designs #Num
```
saved solution in: design.0.pdb
![basepair_end_examples](readme_resources/solution_1.png "RNAMake Solution")

Getting more designs:

```
design_rna -pdb start.pdb -start_bp A222-A251 -end_bp A149-A154 -pdbs -designs 100
> DESIGN RNA: loaded pdb from file: start.pdb
> DESIGN RNA: generated 100 design(s)! if you would like more please specify how many you would like with -designs #Num
```

Controlling size of the solution:

```
design_rna -pdb start.pdb -start_bp A222-A251 -end_bp A149-A154 -pdbs -search.max_size 100
```

![basepair_end_examples](readme_resources/controlling_size.png "Controlling the size of solutions")

Ensuring solutions end with a helix: 

```
design_rna -pdb start.pdb -start_bp A222-A251 -end_bp A149-A154 -pdbs -search.helix_end
```
