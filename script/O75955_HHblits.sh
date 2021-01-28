#run alignement on the database to found templats
hhblits -cpu 4 -i O75955.fasta. -d databases/uniclust30 -oa3m  O75955.a3m

#to add secondary structure information to the MSA
addss.pl  O75955.a3m

#generate a hidden Markov model (HMM) from this MSA
hhmake -i  O75955.a3m

#to convert on fasta file
reformat.pl a3m fas  O75955.a3m  O75955_blits.fas
