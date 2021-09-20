# alice-rivet

Each of our proton-proton jet analyses should contribute a [Rivet](https://rivet.hepforge.org) routine, 
which allows MC generators to be easily compared to our experimental results and used for tuning.
Rivet is built around [HEPData](https://www.hepdata.net) for experimental data, and the HEPMC event format
for MC generators. Rivet uses an internal data format called [yoda](https://yoda.hepforge.org) 
and (annoyingly) uses its own histogramming rather than ROOT or python. 
You can find helpful documentation on the [Rivet Gitlab](https://gitlab.com/hepcedar/rivet)
or [ALICE tutorial](https://alice-doc.github.io/alice-analysis-tutorial/rivet/rivet-tutorial.html).

There are two steps to contribute each new analysis into Rivet:
1. Write a Rivet routine for our analysis.
2. Validate our analysis against an MC sample, and get it approved by ALICE.

# Workflow to write and run a Rivet routine on hiccup

### Create a template
Find your HEPData identifier. For [jet angularities](https://inspirehep.net/literature/1891385) it is `1891385`.
Create a directory in this repo, following the existing examples, then `cd` into it and run:
```
rivet-mkanalysis ALICE_2021_I1891385
```
This command will create several files, including yoda automatically generated from the HEPData entry.
Note that if your analysis is not officially in HEPData, you can upload your submission to the sandbox and download it in yoda format.

### Write your routine
Edit the `.cc` file to write your `init()`, `analyze()`, `finalize()` functions.

### Build your analysis
Enter docker:
```
docker run -it -v /homes/james/alice-rivet:/home/alice-rivet/ hepstore/rivet:3.1.4
```
Then `cd` to your folder and build:
```
rivet-build RivetALICE_2021_I1891385.so ALICE_2021_I1891385.cc
```

### Run over MC events and plot
To run over a single file:
```
rivet --pwd -a ALICE_2021_I1891385 -n 1000 -o Rivet.yoda /home/alice-rivet/generators/herwig/LHC_5020_MPI-S985111.hepmc
rivet-mkhtml --pwd Rivet.yoda
```
To run over multiple files, assume you have a script `validate.sh` that runs over N single files in parallel and produces 
N output files `Rivet1.yoda, Rivet2.yoda, ..., RivetN.yoda`:
```
./validate.sh
yodamerge -o final.yoda Rivet*.yoda
rivet-mkhtml --pwd final.yoda
```
