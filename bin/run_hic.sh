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

make -f ${DIR}/Makefile CONFIG_FILE=$CONFIG_FILE CONFIG_SYS=${DIR}/../install.configure
