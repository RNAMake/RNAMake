RNAMake
=======

**RNAMake** is a toolkit for designing and optimizing RNA 3D structure. It allows 
the alignment between RNA motifs. These motif are small modular peices of RNA that are 
believed to fold independently, thus attaching them together with helix flanking both 
sides allows users of **RNAMake** to build large segments of RNA with a high success 
rate of forming the predicted structure in vitro.

Install
-------

```bash
git clone https://github.com/jyesselm/RNAMake.git
cd RNAMake
sudo pip install -r requirements.txt
```

To install call
```
python setup.py 
```

Compile [optional]
------- 

to compile, make sure you have `cmake` and `ninja` installed with their binaries set up in your `$PATH`. Then run:
```
python compile.py 
```


Tests
-----
To run unit tests for python:
```
RNAMake/rnamake/unittests/run_unittests.sh
```

To run unit tests for C++ code:
```
cd RNAMAke/rnamake/lib/RNAMake/cmake/build/
python run_unittests.py
```


