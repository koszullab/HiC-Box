#!/bin/sh

# recuperer le genome de reference souhaite
# changer le nom du genome

GENOME_NAME=$3

BOWTIE_FOLDER=$1
echo ${BOWTIE_FOLDER}

BOWTIE_BUILD_EXE=${BOWTIE_FOLDER}/bowtie2-build
if [ ! -x "$BOWTIE_BUILD_EXE" ] ; then
    if ! which bowtie2-build ; then
        echo "Could not find bowtie2-build in current directory or in PATH"
        exit 1
    else
        BOWTIE_BUILD_EXE=`which bowtie2-build`
    fi
fi

OLDIFS=$IFS

#fournir ici le chemin d'acces au genome de reference en question

INPUTS=$2




IFS=$OLDIFS

CMD="${BOWTIE_BUILD_EXE} ${INPUTS} ${GENOME_NAME}"
echo Running $CMD
if $CMD ; then
    echo "${GENOME_NAME} index built"
else
    echo "Index building failed; see error message"
fi

#echo "Move the bowtie index to the bowtie folder"
mkdir -p ./index
mv *.bt2 ./index
