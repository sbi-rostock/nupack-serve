<div id="top"></div>

[![Docker Repository on Quay](https://quay.io/repository/bagnacan/nupack-serve/status "Docker Repository on Quay")](https://quay.io/repository/bagnacan/nupack-serve)

# nupack-serve

Nupack-serve is a wrapper for [NUPACK](http://www.nupack.org/): a software
suite for the analysis and design of nucleic acid structures, devices, and
systems [(Wolfe et al. 2017)](https://www.doi.org/10.1021/jacs.6b12693).  
This wrapper modifies some of NUPACK's functions in order to:
- Incoporate NUPACK within computational pipelines
- Enhance the machine-readability of results
- Restrict computation to CPU-bound operations

Nupack-serve was initially created to port NUPACK's *mfe*, *complexes* and
*concentrations* within the [Triplexer](https://github.com/sbi-rostock/triplexer)
pipeline. However, the C functions written to remove I/O bound operation and
provide provenance medatada and JSON output can be adopted by the remaining
NUPACK utilities.  

- [Motivation](#motivation)
- [Practical benefits](#practical-benefits)
  - [CPU-bound execution](#cpu-bound-execution)
  - [Human and machine readability](#human-and-machine-readability)
  - [Reproducibility](#reproducibility)
- [Installation requirements](#installation-requirements)
- [Run nupack-serve](#run-nupack-serve)
- [Usage example](#usage-example)



## Motivation

NUPACK's functionalities rely on the reading/writing of files on the hard
drive. This design principle is handy when results need to be overviewed by
users. However, when placed within the context of computational pipelines, I/O
operations can exert a negative impact on the performances of the underlying
system.  
For the same reason, when results are provided for human-only consultation, it
becomes harder to reuse them in an automated manner. Additional programs need
to be written to parse the generated output files, and this might become an
additional source of errors.  

Nupack-serve is a working proof-of-concept to address these limitations, and
bring NUPACK functionalities to reproducible bioinformatics pipelines.  
This package modifies only three of the many operations provided by NUPACK:
*mfe*, *complexes*, and *concentrations*; which str used within the context of
the [Triplexer](https://github.com/sbi-rostock/triplexer) pipeline to compute
the minimum free energy of secondary structures, the equilibrium base-pairing
properties for each complex in a test tube, and the equilibrium concentration
for each complex.
<p align="right"><a href="#top">&#x25B2; back to top</a></p>



## Practical benefits

Following, we illustrate the benefits of this solution using NUPACL's mfe as
an example.  

Mfe computes the minimum free energy (MFE) secondary structure(s), of
sequence(s) *S* over the ensemble of the complex *C*.  
For multiple strands, mfe's command line invocation should specify the
``-multi`` flag, and provide an input file of the form ``FILE.in`` with the
following information:
1. The number of distinct strand species
2. The sequence for each distinct strand species (each on a separate line)
3. As many integers as the range of 1 to # of *S* sequences, representing the
strand ordering in the complex *C*.  

So, to compute the MFE of the complex resulting from the binding of human
target gene [E2F1](https://www.ncbi.nlm.nih.gov/nuccore/NM_005225)'s binding
site with the two putatively cooperating miRNAs
[hsa-miR-205](http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MIMAT0000266)
and [hsa-miR-342-3p](http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MIMAT0000753),
we have to provide the following ``test.in`` file:
```
3
ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu
uccuucauuccaccggagucug
ucucacacagaaaucgcacccgu
1 2 3
```

and call mfe with the following command line invocation:
```
$ mfe -multi test
```

NUPACK stores the result of mfe's call in a file called ``test.mfe``, which
contains the following lines:
```
% NUPACK 3.2.2
% Program: mfe
% Start time: Thu Apr 23 16:32:20 2020
% Command: mfe -multi test
% Sequence:  ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu+uccuucauuccaccggagucug+ucucacacagaaaucgcacccgu
% v(pi): 1
% Parameters: RNA, 1995
% Dangles setting: 1
% Temperature (C): 37.0
% Sodium concentration: 1.0000 M
% Magnesium concentration: 0.0000 M

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
99
-50.563
.((((.((((.(.(((((((((((((......)))))..((((((((((((((.+.)))))))))).))))......+.))))))))..).)))).)))).
2   98
3   97
4   96
5   95
7   93
8   92
9   91
10  90
12  88
14  85
15  84
16  83
17  82
18  81
19  80
20  79
21  78
22  37
23  36
24  35
25  34
26  33
40  70
41  69
42  68
43  67
44  65
45  64
46  63
47  62
48  61
49  60
50  59
51  58
52  57
53  56
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
```

This flow of execution presents limitations that we address using nupack-serve:
<p align="right"><a href="#top">&#x25B2; back to top</a></p>



### CPU-bound execution

If we want to compute one million MFE calculations, we end up having to provide
one million ``.in`` files, and obtain another million ``.mfe`` files as result.
This behavior heavily impacts on performances.  
For this reason, we force mfe to avoid I/O operations on disk.  

```
$ mfe -multi
Requesting input manually.
Enter number of strands: 3 <--- automated input
Enter sequence for strand type 1:
ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu <--- automated input
Enter sequence for strand type 2:
uccuucauuccaccggagucug <--- automated input
Enter sequence for strand type 3:
ucucacacagaaaucgcacccgu <--- automated input
Enter strand permutation (e.g. 1 2 4 3 2):
1 2 3 <--- input
```
<p align="right"><a href="#top">&#x25B2; back to top</a></p>



### Human and machine readability

As humans, we are able to immediately tell which parts of the resulting
``.mfe`` file are reserved to report the MFE calculation, and which to the
metadata of the whole operation. However to retrieve specific details of each
run, *e.g.* NUPACK's version, the resulting complex's sequence, or the
simulation's temperature, we would need to write dedicated parsers to run against
each ``.mfe`` file. This represents an additional source for errors.  
For this reason, we instruct mfe to return results and provenance metadata in the
JSON format.

```
{
  "version": "3.2.2",
  "command": "mfe -multi",
  "cutoff": 0.001,
  "sequence": "ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu+uccuucauuccaccggagucug+ucucacacagaaaucgcacccgu",
  "dangles": "1",
  "temperature (C)": 37.0,
  "concentration Na (M)": 1.0000,
  "concentration Mg (M)": 0.0000,
  "pseudoknots": false,
  "sequence length (nt)": 99,
  "minimum free energy (Kcal/mol)": -50.563,
  "dot-bracket": ".((((.((((.(.(((((((((((((......)))))..((((((((((((((..)))))))))).)))).......))))))))..).)))).))))",
  "pairs": [
    (2,98),(3,97),(4,96),(5,95),(7,93),(8,92),(9,91),(10,90),(12,88),(14,85),(15,84),(16,83),(17,82),(18,81),(19,80),(20,79),(21,78),(22,37),(23,36),(24,35),(25,34),(26,33),(40,70),(41,69),(42,68),(43,67),(44,65),(45,64),(46,63),(47,62),(48,61),(49,60),(50,59),(51,58),(52,57),(53,56)
   ]
 }
```
<p align="right"><a href="#top">&#x25B2; back to top</a></p>



### Reproducibility

NUPACK has a list of software dependencies and requires users to compile its
source code. While this doesn't represent a problem for the computer savvy
users who have full control over their software installations, many do not fall
within this category.  
For this reason, we provide NUPACK's functions within a docker container.
<p align="right"><a href="#top">&#x25B2; back to top</a></p>



## Installation requirements

The only requirement is [Docker](https://www.docker.com/), which can be
installed in different ways depending on the underlying operative system:
- Unix users should follow the [Docker installation for Linux](https://docs.docker.com/compose/install/#install-compose-on-linux-systems#install-compose-on-linux-systems),
and install both Docker and Docker compose
- MacOS 10.13+ users should follow the [Docker installation for Mac](https://docs.docker.com/docker-for-mac/install/)
- Windows 10+ users, should follow the [Docker installation for Windows](https://docs.docker.com/docker-for-windows/install/)
- For legacy systems, users can rely on the [Docker Toolbox](https://docs.docker.com/toolbox/overview/).
<p align="right"><a href="#top">&#x25B2; back to top</a></p>



## Run nupack-serve

To run nupack-serve, you need to run its docker container. Type:
```
docker run -d -p 9000:8000 quay.io/bagnacan/nupack-serve
```

The nupack-serve app is now listening on your host's port 9000.  
You can test its docker container by trying one of the provided examples.  
In your browser, type:
```
localhost:9000
```

You should obtain something like this:
```
{
  "usage":"Nupack-serve is a wrapper for NUPACK: a software suite for the analysis and design of nucleic acid structures, devices, and systems (Wolfe et al. 2017).",
  "homepage":"https://github.com/sbi-rostock/nupack-serve",
  "license":"NUPACK Software License Agreement" <-- here cut for brevity
}
```

You can now:
- Incoporate NUPACK within computational pipelines
- Enhance the machine-readability of results (JSON)
- Avoid read/writes to disk
<p align="right"><a href="#top">&#x25B2; back to top</a></p>



## Usage example

Nupack-serve wraps NUPACK's *mfe*, *complexes* and *concentrations* in a python
application which returns:
1. NUPACK's license agreement
2. The exit status of the subprocess which executed the called NUPACK function
3. Results and provenance metadata in the JSON format

With reference to the example [above](#practical-benefits), you can compute
the MFE of human target gene E2F1's binding site sequence with the two
putatively cooperating miRNAs hsa-miR-205 and hsa-miR-342-3p.  

In your browser, type:
```
localhost:9000/example/mfe
```

You should obtain something like this:
```
{
  "license":"NUPACK Software License Agreement", <-- here cut for brevity
  "status":0, <-- exit status of the subprocess that executed mfe in the container
  "result":{ <-- mfe result
    "version":"3.2.2",
    "command":"/usr/local/bin/mfe -multi",
    "cutoff":0.001,
    "sequence":"ccgggggugaaugugugugagcaugugugugugcauguaccggggaaugaaggu+uccuucauuccaccggagucug+ucucacacagaaaucgcacccgu",
    "dangles":"1",
    "temperature (C)":37.0,
    "concentration Na (M)":1.0,
    "concentration Mg (M)":0.0,
    "pseudoknots":false,
    "sequence length (nt)":99,
    "minimum free energy (Kcal/mol)":-50.563,
    "dot-bracket":".((((.((((.(.(((((((((((((......)))))..((((((((((((((..)))))))))).)))).......))))))))..).)))).))))",
    "pairs":[
      [2,98],[3,97],[4,96],[5,95],[7,93],[8,92],[9,91],[10,90],[12,88],[14,85],[15,84],[16,83],[17,82],[18,81],[19,80],[20,79],[21,78],[22,37],[23,36],[24,35],[25,34],[26,33],[40,70],[41,69],[42,68],[43,67],[44,65],[45,64],[46,63],[47,62],[48,61],[49,60],[50,59],[51,58],[52,57],[53,56]
    ]
  }
}
```
<p align="right"><a href="#top">&#x25B2; back to top</a></p>
