[![Build Status](https://travis-ci.com/RNAMake/RNAMake.svg?token=Kxiycasibo9yqZt7eayf&branch=devel)](https://travis-ci.com/RNAMake/RNAMake)

RNAMake
=======
-----------
(c)  J.Yesselman, Stanford University, 2014-2018

**RNAMake** is a toolkit for designing and optimizing RNA 3D structure. It allows
the alignment between RNA motifs. These motif are small modular pieces of RNA that are
believed to fold independently, thus attaching them together with helix flanking both
sides allows users of **RNAMake** to build large segments of RNA with a high success
rate of forming the predicted structure _in vitro_.

## Development

Trying to make code changes to `RNAMake`? Check out the [development guide](CONTRIBUTING.md).


Install
-------
```shell
git clone https://github.com/jyesselm/RNAMake.git
cd RNAMake
```

Make sure to have a valid c++ compiler either g++ (> 4.6) or clang as well as python 3.0+. In your .bashrc (.bash_profile if mac OSX) add:

```shell
# location to RNAMake directory used by RNAMake
export RNAMAKE=<RNAMake Path>
```
[OPTIONAL] but useful

```shell
# location of RNAMake c++ executables
export PATH=$PATH:$RNAMAKE/rnamake/lib/RNAMake/cmake/build/
# location of python executable scripts
export PATH=$PATH:$RNAMAKE/rnamake/bin/
```

note if you are using c shell or another non bash shell you will need to use the equivalent commands


Compile
-------
requires `cmake` and `ninja`

cmake: https://github.com/Kitware/CMake <br>
ninja: https://github.com/ninja-build/ninja

make sure they are in your `$PATH`. Then run:

```shell
cd RNAMake/cmake/build_XXX
# where XXX is your compiler of choice, for clang
cd RNAMake/cmake/build_clang
python ../make_project.py -target OS {mac|linux}
cmake -G Ninja
ninja
```

Applications
============
1. <a href="#design_rna_scaffold">design_rna_scaffold</a>

<a id="design_rna_scaffold"></a>

design_rna_scaffold
-----------
Generates segments are RNA between two Watson-Crick basepairs can also perform sequence optimization for helical sequences.

```
DesignRNAScaffold
Usage: ./design_rna_scaffold [OPTIONS]

Options:
  -h,--help                                                                     Print this help message and exit

Core Inputs:
  --pdb TEXT:(FILE) AND (ends in .pdb)                                          path to a PDB file with input RNA structure
  --start_bp TEXT:format [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]                starting basepair to be used in structure format: [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]
  --end_bp TEXT:format [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]                  ending basepair to be used in structure format: [CHAIN ID][NT1 NUM]-[CHAIN ID][NT2 NUM]
  --designs INT:POSITIVE=1                                                      number of designs to create. Default is 1
  --log_level TEXT:{debug,error,fatal,info,verbose,warn}=info                   level for global logging
  --extra_pdbs TEXT                                                             , deliminted list of other pdbs used in building

I/O Options:
  --dump_intermediate_pdbs                                                      flag to dump intermediate pdbs
  --dump_pdbs                                                                   TODO
  --new_ensembles TEXT                                                          flag to include new structural ensembles
  --no_out_file=0                                                               if you only want the summary and not the actual structures
  --out_file TEXT=default.out                                                   output file that contains all information to rebuild solutions
  --score_file TEXT=default.scores                                              name of output file containining scoring information for design


Search Parameters:
  --ending_helix TEXT                                                           ending helix for design solution. Format = [TODO]
  --exhaustive_scorer TEXT=default                                              TODO
  --max_helix_length INT=99                                                     maximum number of basepairs in a solution helix
  --mc_scorer TEXT:{default,scaled_scorer}=default                              TODO
  --min_helix_length INT=4                                                      minimum number of basepairs in a solution helix
  --motif_path TEXT                                                             TBD
  --no_basepair_checks                                                          flag to disable basepair checks
  --scaled_score_d FLOAT=1                                                      TODO
  --scaled_score_r FLOAT=2                                                      TODO
  --search_cutoff FLOAT=7.5                                                     TODO
  --search_max_size INT:POSITIVE=999999                                         maximum number of steps for a design search
  --search_type TEXT:{exhaustive,mc,path_finding}=path_finding                  search type for traversing motif space
  --solution_filter TEXT:{NoFilter,RemoveDuplicateHelices}=RemoveDuplicateHelices
                                                                                TODO
  --starting_helix TEXT                                                         starting helix for design solution. Format = [TODO]
    --no_sterics                                                                turns off sterics checks againsts supplied RNA structure
  --only_tether_opt                                                             ignore supplied structure other than sterics

Sequence Optimization Parameters:
  --skip_sequence_optimization                                                  flag to skip sequence optimization of the design
  --sequences_per_design INT=1                                                  number of sequences to try per motif design

Thermo Fluc Parameters:
  --thermo_fluc                                                                 run thermo fluc procedure to estimate thermo fluc of helices

```
#### Required Arguments

Argument  | Description
------------- | -------------
--pdb		    | PDB file containing starting RNA, make sure that both the start basepair and end basepair are watson and crick base pairs that are end at RNA chains. Must supply both `-start_bp` and `-end_bp`. See figure below.
--start_bp			    | The Watson-Crick basepair to start building the RNA segment from. Example "A194-A252" the base pair between resiudes 192 and 252 both on chain A.
--end_bp			 |	The Watson-Crick basepair to end the RNA segment. Same naming convention as `start_bp`. See examples below.

#### Example of Basepair Ends that can be connected with `design_rna`

![basepair_end_examples](readme_resources/ggaa_tetraloop.png "Basepair End Example")

##Examples

examples are located: $RNAMAKE/examples/cpp/design_rna


start.pdb:
![basepair_end_examples](readme_resources/ggaa_tetraloop.png "Basepair End Example")

examples in $RNAMAKE/tests/design_ran_scaffold/

each directory has a cmd/COMMAND which gives a functioning command to explore the different options available.

