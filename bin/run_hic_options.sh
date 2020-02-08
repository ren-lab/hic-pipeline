#!/bin/bash
#######################################################################
### Copyleft (c) 2017 Bing Ren Lab
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#######################################################################


function usage(){
echo -e "Usage: $0 -c CONFIG_FILE -s server" 
echo -e "\t-c [config_file]: configuration file "
echo -e "\t-e [email]: email address."
echo -e "\t-s [server]: silencer or TSCC"
echo -e "\t-o [options]: valid_pairs or vanilla or all"
exit 1
}



### check options given to this command
while getopts "c:s:e:o:h:" OPT
do
    case $OPT in
    c) CONFIG_FILE=$OPTARG;;
    e) EMAIL=$OPTARG;;
    s) SERVER=$OPTARG;;
    o) OPTION=$OPTARG;;
    h) help ;;
    \?)
         echo "Invalid option: -$OPTARG" >&2
         usage
         exit 1
         ;;
     :)
         echo "Option -$OPTARG requires an argument." >&2
         usage
         exit 1
         ;;
    esac
done

if [ $# -eq 0 ]; then usage ;exit 1; fi

if ! [ -e $CONFIG_FILE ]; then echo File $CONFIG_FILE not exist; exit 1; fi 

if [ -z ${SERVER+x} ];
  then echo -e "Please tell us the server, eg. silencer, TSCC"; usage;exit; fi

NTHREADS=30
DIR=$(dirname $0)
LOG=run-$(date +%Y-%m-%d-%H-%M-%S).log
## validate the programs are installed.
. ${DIR}/validate_programs.sh

## load snakemake environment for Renlab
if [ $SERVER == "silencer" ]; then
  source /projects/ps-renlab/share/Pipelines/environments/python3env/bin/activate
  ### unlock the directory
  touch Snakefile
  snakemake --unlock
  rm Snakefile
  ## start analysis
  echo "$(date) # Analysis Began" > $LOG
  nice -n 19 snakemake $OPTION -p -k --ri --snakefile ${DIR}/../scripts/Snakefile \
  --configfile $CONFIG_FILE --cores $NTHREADS \
  --config BWA_INDEX_PATH=/projects/ps-renlab/share/bwa_indices/ \
  2> >(tee -a $LOG >&2)
  echo "$status"
  echo "$(date) # Analysis finished" >> $LOG
  [[ $EMAIL =~ @ ]] && (
  echo "See attachment for the running log.
  Your results are saved in:
  $(pwd)"  | mail -s "Hi-C analysis done" -a $LOG  $EMAIL )

elif [ $SERVER == "TSCC" ]; then
  ### load modules.
  module load python
  unset $PYTHONPATH
  #source /projects/ps-renlab/share/Pipelines/environments/python3env_TSCC/bin/activate
  ### unlock the directory
  touch Snakefile
  snakemake --unlock
  rm Snakefile
  ## started analysis
  if [ ! -d pbslog ]; then mkdir pbslog; fi
    echo "$(date) # Analysis Began" > $LOG
  snakemake $OPTION --snakefile ${DIR}/../scripts/Snakefile -p  -k -j 1000 --ri \
  --configfile $CONFIG_FILE \
  --config BWA_INDEX_PATH=/projects/ps-renlab/share/bwa_indices/ \
  --cluster-config ${DIR}/../cluster.json \
  --cluster "qsub -l nodes=1:ppn={threads},walltime={cluster.time} -N {rule} -q hotel -o pbslog/{params.pbsName}.{rule}.pbs.out -e pbslog/{params.pbsName}.{rule}.pbs.err" \
  --jobscript ${DIR}/../scripts/jobscript.pbs --jobname "{rulename}.{jobid}.pbs" \
  2> >(tee -a $LOG >&2)
  echo "$(date) # Analysis finished" >> $LOG
  [[ $EMAIL =~ @ ]] && (
  echo "See attachment for the running log.
  Your results are saved in:
  $(pwd)"  | mail -s "Hi-C analysis done" -a $LOG  $EMAIL
  )
else 
  echo -e "Invalid server option: $server"; exit 1;
fi 

