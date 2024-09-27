#!/bin/bash

#get command line arguments
while getopts ":f:t:c:h" flag; do
    case $flag in
        f) seqFile=$OPTARG
        echo "Using '$seqFile' as list of sequencing files to process"
        ;;
        t) targetFile=$OPTARG
        echo "Using '$targetFile' as list of spacer targets for which to search for hits"
        ;;
        c) colorFile=$OPTARG
        echo "Using '$colorFile' for the colors for the Chord Diagram that will be drawn"
        ;;
        h) echo ""
        echo "This script will run the Docker container an analyze MuRiCS data to identify"
        echo "unique spacer combinations using a program called 'auto_murcis.py'."
        echo ""
        echo "This requires three input files from the user:"
        echo "    '-f' to inidicate the sequencing files"
        echo "    '-t' to inidicate the spacer target file"
        echo "    '-c' to indicate the Chord Diagram color file"
        echo "    '-h' will call this detail descirption"
        echo ""
        echo "This requires Docker and the auto_murcis Docker image"
        echo ""
        echo "Please contact Kevin Myers (kmyers2@wisc.edu) for help."
        echo ""
        exit 1
        ;;
        \?) 
        echo ""
        echo "Invalid option: -$OPTARG."
        echo "Please use the following flags:"
        echo "    '-f' to inidicate the sequencing files"
        echo "    '-t' to inidicate the spacer target file"
        echo "    '-c' to indicate the Chord Diagram color file"
        echo "Quiting."
        echo ""
        exit 1
        ;;
        :)
        echo ""
        echo "Option -$OPTARG requires an argument. Quiting."
        echo "" 
        exit 1
    esac
done

cwd=$(pwd)

copySeqFile="${cwd}/${seqFile}"
copyTargetFile="${cwd}/${targetFile}"
copyColorFile="${cwd}/${colorFile}"

docker run -dt --name auto_murcis kevinmyers/auto_murcis
docker cp $copySeqFile auto_murcis:/var/auto_murcis/
docker cp $copyTargetFile auto_murcis:/var/auto_murcis/
docker cp $copyColorFile auto_murcis:/var/auto_murcis/
while read line; do docker cp ${cwd}/$line auto_murcis:/var/auto_murcis/; done < $seqFile
#docker run conda activate auto_murcis_env
docker exec auto_murcis /bin/bash -c "python /var/auto_murcis/murcs_script.py -f $seqFile -t $targetFile -c $colorFile"
#set date and output directory to use
dateVar=$(date +%d-%m-%Y)
outputDir=$"auto_murcis_output"-$dateVar
#copy output directory from Docker container to local drive
docker cp auto_murcis:/var/auto_murcis/$outputDir ./
docker container stop auto_murcis
docker rm auto_murcis