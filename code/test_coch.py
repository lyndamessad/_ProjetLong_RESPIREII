import os
import shutil
import subprocess
import os.path as op
from module import *

def prep_folder(seq, ident, root_dir):
    """Take in a sequence string and prepares the folder for the I-TASSER run."""
    itasser_dir = op.join(root_dir, ident)

    if not op.exists(itasser_dir):
        os.makedirs(itasser_dir)

    tmp = {ident: seq}
    with open("cotch.txt", "w") as file: 
        file.write("indict=tmp,outname='seq',outext='.fasta',outdir=itasser_dir")
    return itasser_dir
prep_folder("KKJEFHEZFIH","O75955.fasta", '/home/lynda/Bureau/Projet_long/code/')
thedir = '../Proteines'
print(thedir)
folders = [name for name in os.listdir(thedir)]
folders = sorted(folders, reverse=True)
length = int(len(folders)) / int(2)
print(length)
for_ssb3 = folders[:int(length)]
#
for fo in for_ssb3:
    coach = open('%s_coach.sh' % fo, 'w')
    coach.write('#!/bin/bash\n')
    coach.write('#PBS -N %s\n' % fo)
    coach.write('-pkgdir /home/lynda/I-TASSER5.0/I-TASSERmod/runI-TASSER.pl \
-libdir /home/lynda/ITLIB \
-seqname %s \
-model model1.pdb \
-datadir /home/lynda/I-TASSER5.0/%s \
-GO true \
-runstyle parallel \
-ligth \n\n '% (fo, fo))
    coach.close()
    subprocess.call('qsub %s_coach.sh;' % (fo), shell=True)
    print('qsub %s_coach.sh;' % (fo)),
