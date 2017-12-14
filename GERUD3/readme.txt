Compiling GERUD3.0 on Linux:

Ensure that all of the header files and the source code file are in your current directory.
The required files are:

Gerud3.cpp
AlleleFrequencyData.h
GS3functionsclasses.h
MTwisterFunctions.h
G3LoadProgArray.h
 
Type:

g++ Gerud3.cpp -o Gerud3

You may need to alter the permission to make it executable:

chmod u+x Gerud3


To run Gerud3, you need an input file (progeny array). It's advisable to also have an allele
frequency file, which contains population-level estimates of allele frequencies. 
Ensure the following:

(1) The locus names are identical, including case, between files.
(2) Every allele present in the progeny array is also present in the allele frequency file.
(3) The allele frequencies in the allele frequency file sum to 1.
(4) The progeny array file has no missing data. Offspring with missing genotypes should be dropped.

See the sample files for the format. The format for the allele frequency file is:

A first line indicating the number of loci: 4 loci

A second line with the name of the first locus, preceded by an exclamation point: !locus1

A third line with the designation of the first allele, a tab character, 
and the frequency: 165[tab]0.00171.

Additional lines for all of the allele frequencies for this first locus.

Data for the second locus are in the same file, and they start with the locus name, 
preceded by an exclamation point: !locus2.

Repeat this format for every locus.


PROGENY ARRAY FILE:

The first line indicates whether the mother is known: Known Mother or Unknown Mother
The second line indicates the number of embryos: 1000 embryos
The third line indicates the number of loci: 2 loci
The fourth line indicates the maternal genotype if it is known: Mom[tab]100/102[tab]222/226
The remaining lines indicate embryo genotypes: Emb_1[tab]100/110[tab]226/232

See the sample file for an example.


GERUD3 USAGE:

Put the GERUD3 compiled executable in the folder with your allele frequency file and progeny array file (or add it to your path).

./Gerud3 -i progenyarrayfile.txt -o progenyarrayoutput.txt -a allelefrequencies.txt -v n

Arguments:

-i: name of your unput file (default: Gerud3progenyarray.txt)
-o: name of your output file (default: Gerud3output.txt)
-a: name of your allele frequency file (default: Gerud3allelefrequencies.txt)
-v: verbose (Y or N) (default: yes)

The Gerud3 output will be saved in your output file. It will overwrite any existing file of the same name, so be careful!


