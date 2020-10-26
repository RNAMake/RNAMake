# build_sqlite_libraries status
So far have only focused on `/motif_libraries_new/` files. 
## Files where the python script's results disagree with the original 
+ bp_steps.db
+ new_bps_steps
+ unique_twoway

These seem to differ because the clustering algorithm being used is dependent on the order in which 
the motifs are seen. The way the code is structured, there is not a definite ordering of how the 
motifs are ordered when they are given to the clustering algorithm.

## Files where this application disagrees with the python script's output
+ flex_helices.db => this file is not created because no motifs can be made. This makes use of the 
`motif::Motif(String const &,structure::ResidueTypeSet const &)` ctor and seems to throw an exception 
from the constructor because many of the basepairs are `nullptr` (not sure if this is the actual cause
though). 
+ bp_steps.db, new_bp_steps.db, unique_twoway.db => These files similiarly disagree with the original 
files but also with the python script. Likely due to variance introduced by the clustering algorithm.
+ tcontact.db, nway.db and hairpin.db => see some but many fewer misses. Differences between the C++ and
python starts to creep up when the motif libraries are loaded in. Namely, the number of loaded motifs 
for each type differs and there are a number of x3dna-dssr segfaults when the code attempts to generate
into the new motifs. (These problems are all found in `BuildSqliteLibraries::_build_basic_libraries()`)

## Usage
build_sqlite_libraries has two subcommand run modes: BUILD and VALIDATE
To BUILD, use the following:

`./build_sqlite_libraries BUILD --in_dir /path/to/motifs/folder --out_dir /path/to/new_dir/ --all`
+ `--in_dir` Specifies the folder where the motif folders are.
+ `--out_dir` Specifies where the new db files will be put. Note that the code will setup the 
appropriate subfolders.

To VALIDATE, use the following:

`./build_sqlite_libraries VALIDATE [--orig /path/to/orig/db/files] --remake /path/to/new/db/files [--summary_file name.csv]`

+ `--orig` Refers to the directory where the "true" db files will be found. Note that it is optional and by default is 
set to the resources folder.
+ `--remake` Path to the new db files that have just been made. 
+ `--summary_file` Optional, creates a .csv file that includes the names of original motifs that were NOT found or aligned
within the remade set.

**NOTE: `VALIDATE` mode ONLY looks at files in `/motif_libraries_new/`.**

## Other
+ a lot of errors seem to tie back to there not being a 1:1 corresponding function for the python `motif.str_to_motif()`
function
