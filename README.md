# Eternabot 2.0

Eternabot 2.0 is an update algorithm for predicting the likelyhood of a RNA sequence will fold into a target secondary structure or base pairing arrangement.

### How to install

Requirements:

c++ compiler that supports c++14

clang 3.4 or greater (successfully compiled on clang 11)

gcc 5 or greater (successfully compiled on gcc 8)



Compiling tool requirements: 

cmake >= 3.00 

ninja >= 1.90

python >= 2.7



**Note all commands below use `$` before a command to indicate its at command line. You do not need to include this** 

```bash
$ cd cmake/build
$ python ../make_project.py 
$ cmake -G Ninja
$ ninja 
```



### How to run 

after successful compiling there should be an executable 'eternabot' in ```cmake/build``` 

you can add this directory to your path so you can run ```eternabot``` anywhere

```bash
# a standard run 
$ eternabot -seq "SSSSWNNNNWSSSS" -ss "(((((....)))))"
96.186 0.0434825 GCGGUGAAAACCGC (((((....))))) (((((....)))))
```

where:

 `-seq` is the sequence constraint string based on the IUPAC notation: https://en.wikipedia.org/wiki/Nucleic_acid_notation

`-ss` is the target secondary structure in dot bracket notation: http://eternawiki.org/wiki/index.php5/Dot-Bracket_Notation



the output currently contains 5 elements per solution (elements from example above)

96.186: the predicted eterna score out of 100. 100 would be the sequence is predicted to perfectly match the target secondary structure.

0.0434825: the difference in base pair probability between the predicted structure by vienna 1.8.5 and the target structure. 0 would be no difference at all. Generally numbers under 1 are as good as you can get.

GCGGUGAAAACCGC: the sequence solution for the target secondary structure and design constraints.

(((((....))))): The target structure

(((((....))))): the predicted structure, in this case the same. Predicted by vienna 1.8.5. 



There is also an outputed csv file, which by default is eternabot.csv.

```ba
$ cat eternabot.csv
opt_num,structure,opt_score,bp_diff_score,opt_sequence,longest_gc_stretch
0,(((((....))))),96.186,0.0434825,GCGGUGAAAACCGC,4
```

Which is a lot of the same information.



### other options 











