wd=test3

if [ -d $wd ]; then
    rm -r $wd
fi

# The expected runtime of this test is ~2:15 minutes. The fast runtime of
# this test is mostly caused by generating a low number of training genes.
# Note that this approach does not scale well with increasing genome size
# and the number of proteins in a protein database. The runtime on a full
# genome will be much slower than with the command used in test1.sh.
# Also, --skipOptimize is here included to speed up the test. Never enable
# this flag for annotating a real genome!


( time galba.pl --genome=../genome.fa --prot_seq=../proteins.fa --skipAllTraining --softmasking --workingdir=$wd --threads 8 --species=arabidopsis) &> test3.log
