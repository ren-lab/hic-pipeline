trap "echo aha broken pipe" SIGPIPE
#set -e
set -o pipefail
if [ -e pip ]; then rm pip; fi
mkfifo pip
#bzcat 20M.txt.bz2 > pip &
#tail pip &
bzcat 20M.txt.bz2 > pip | tail pip
echo $?
sleep 3
echo done

