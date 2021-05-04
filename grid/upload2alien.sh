#!/bin/bash

#alien-token-init

aHOME="/alice/cern.ch/user/h/hscheid"                   # alien home directory
DIR="test/Delphes"                        # directory where to store the FILEs

LOCALDIR=`pwd`                                            # local workdir
#LOCALDIR="$ALICE_ROOT/$DIR/$FILE"                      # local workdir
# FILEs to upload

# colors
Cm="\033[35m"
Cy="\033[33m"
Cc="\033[36m"
Cb="\033[34m"
Cg="\033[32m"
Cr="\033[31m"
Cw="\033[37m"
Cz="\033[m"
Br="\033[41m"
By="\033[43m"

FILELIST=$*

echo -e "${Cm}files to upload: ${Cb}$FILELIST${Cz}"
echo -e "${Cm}destination: $aHOME/$DIR ${Cz}"
####################

I=0
for FILE in $FILELIST ; do

  if [ -f "$LOCALDIR/$FILE" ] ; then
    echo -e "${Cg}uploading file $I: ${Cb}$FILE${Cz}"
    alien_rm "$aHOME/$DIR/$FILE"
	  alien_cp "$LOCALDIR/$FILE" alien://"$aHOME/$DIR/$FILE"
    let I=I+1
  fi

done

echo -e "${Cg}files now present on alien:${Cz}"
alien_ls -l "$aHOME/$DIR"
