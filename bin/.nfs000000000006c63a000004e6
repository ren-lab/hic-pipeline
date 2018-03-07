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
echo -e "Usage: $0 -c CONFIG_FILE" 
}


while getopts "c:h:" OPT
do
    case $OPT in
    c) CONFIG_FILE=$OPTARG;;
    s) JOB_SCRIPT=$OPTARG;;
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

DIR=$(dirname $0)

#make -f ${DIR}/Makefile CONFIG_FILE=$CONFIG_FILE CONFIG_SYS=${DIR}/../system.configure
command -v snakemake >/dev/null 2>&1 || { echo >&2 "snakemake is not installed";exit 1; }

if [ -z ${JOB_SCRIPT+x} ]; then 
  source /mnt/silencer2/share/Piplines/environments/python3env/bin/activate
  snakemake vanilla  -p --snakefile ${DIR}/../scripts/Snakefile2 --configfile $CONFIG_FILE --cores 30
else 
  snakemake --snakefile ${DIR}/../scripts/Snakefile --configfile $CONFIG_FILE -p  -k -j 500 --cluster "qsub -l nodes=1:ppn={threads} -N {rule} -q hotel -o pbslog/{rule}.pbs.out -e pbslog/{rule}.pbs.err" --jobscript $JOB_SCRIPT --jobname "{rulename}.{jobid}.pbs"
fi 


