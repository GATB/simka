

import os, sys

if sys.argv[1] == "reset":
    os.system("rm -rf /local/output/simka2_database/")
    os.system("python /Users/gbenoit/workspace/gits/gatb-simka/scripts/simka2/core/simka2_init.py -database-dir /local/output/simka2_database -kmer-size 31 -abundance-min 0")
else:
    pass

#os.system("python /Users/gbenoit/workspace/gits/gatb-simka/scripts/simka2/simka2-count.py -in /Users/gbenoit/workspace/gits/gatb-simka/example/simka_input_N-5.txt -out-tmp /local/output/simka2_test_tmp/ -database /local/output/simka2_database -simka-bin /Users/gbenoit/workspace/gits/gatb-simka/build/bin/")
os.system("python /Users/gbenoit/workspace/gits/gatb-simka/scripts/simka2/core/simka2-count.py -in ../example/simka_input.txt -out-tmp /local/output/simka2_test_tmp/ -database /local/output/simka2_database -simka-bin /Users/gbenoit/workspace/gits/gatb-simka/build/bin/")
#/local/data/metagenomic/input_metasoil_simka
#../example/simka_input.txt

command = "python ../scripts/simka2/core/simka2-distance.py -database-dir /local/output/simka2_database/ -simka-bin ./bin/"
os.system(command)

#./bin/simka2-distance-final -nb-partitions 200 -database-dir /local/output/simka2_database/ -kmer-size 31