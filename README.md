To build mipfinder2

You need to have g++ compiler that supports C++17 installed.

Navigate to a directory where you want to have mipfinder2

Clone the mipfinder2 code by running the following command:
git clone https://github.com/ku-mip/mipfinder2.git

This will create a new folder in the current directory called "mipfinder2. Navigate to that folder and create a new folder called "build". Navigate to this folder using these commands:

cd mipfinder2
mkdir build
cd build

Create the build instructions for CMake using this command:

cmake ..

Now compile mipfinder2 using the following command:

make

After the compilation is done, you should see a new folder appear in mipfinder2 called "bin". In there you will find "mipfinder".

To run mipfinder2:

You need to make sure you have hmmsearch, hmmbuild, phmmer and clustalo installed and accessible through the PATH variable. Easiest way is to install the HMMER package and the clustalo package.

You need to provide the following files:
 - The proteome to search, in FASTA format. Set the "organism_proteins_fasta" in the configuration to this file.
 - A file of known microproteins from that organism in FASTA format. If there are none known, you can create an empty file, but it has to exist on the disk. Set the "known_mips_fasta" in the configuration to this file.
 - A list of InterPro domains and their types. You can download it from here: https://www.ebi.ac.uk/interpro/download/ , choose the "InterPro entry list" tsv file. Set the "interpro_database" in the configuration to this file.
 - Gene Ontology annotation file, you can download it from here: http://purl.obolibrary.org/obo/go.obo. Set the "go_database" in the configuration to this file.
 - A file that maps UniProt accessions to InterPro accessions. This has to be a tab-separated file where the first column is UniProt accession, the second column is the UniProt accession sequence version (i.e. which isoform it is) and the third column is a comma (;) separated list of InterPro identifiers.
You can create this by creating a custom query on UniProt for your organism of choice and downloading the results, e.g. for Arabidopsis thaliana see the following link: https://www.uniprot.org/uniprotkb?facets=reviewed%3Atrue&query=%28taxonomy_id%3A3702%29&fields=accession%2Csequence_version%2Cgo_id&view=table
 Set the "uniprot_to_interpro" in the configuration to this file.
 - A file that maps UniProt accessions to GO (Gene Ontology) accessions. This is in the same format as the above InterPro file. Set the "uniprot_to_go" in the configuration to this file. You can do this the same way as the InterPro data above, i.e. https://www.uniprot.org/uniprotkb?facets=reviewed%3Atrue&query=%28taxonomy_id%3A3702%29&fields=accession%2Csequence_version%2Cgo_id&view=table

Set the "organism_identifier" to whatever name you want to provide to the results.

Run the compiled mipfinder2 executable and provide the configuration.ini as the first command-line argument, e.g. if the executable and the file are in the same folder, just run:

./mipfinder configuration.ini

The results will be in the results_folder/final_results.txt