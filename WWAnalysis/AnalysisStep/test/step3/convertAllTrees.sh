#!/usr/bin/env bash

usage()
{
if [ $# -gt 0 ] ; then a="Error: $*"; fi
cat << EOF
$a
usage: $0 options file1 [file 2 file 3 ...]

mwl Wrapper for CRAB

OPTIONS:
   -h      Show this message
   -v      Verbose
   -d      Output directory [default = .]
   -v      Verbose
   -t      Do not execute anything ... testing only

EOF
}

scriptPath=$(dirname $(readlink -f $0))

set -- `getopt "hvd:vt" "$@"`

VERBOSE=
TEST=
OUTPUTDIR="."
while [ "$1" != "--" ] ; do
  case "$1" in
    -h) usage; exit 1;;
    -v) VERBOSE=1;;
    -d) OUTPUTDIR=$2; shift;;
    -t) TEST="1";;
    -v) VERBOSE=1;;
  esac
  shift
done
shift


if [ $# -lt 1 ] ; then
    usage "I would appreciate something to do ... ask tar what they think."
    exit 3
fi

pids=
[ -n "$VERBOSE" ] && echo "Output directory: $OUTPUTDIR"
[ -n "$VERBOSE" ] && echo "Input files: $*"

for x in $*; do 
    [ -n "$VERBOSE" ] && echo "Converting `basename $x` from directory `dirname $x`"
    [ -z "$TEST"    ] && echo "Processing: $x" 
    [ -z "$TEST"    ] && root -b -l -q $scriptPath/convertToLatinos.C\(\"`dirname $x`\",\"`basename $x`\"\) > /dev/null 2>&1 &
    [ -z "$TEST"    ] && pids="$pids $!"
    [ -z "$TEST"    ] && { while [[ `ps | grep "root\.exe" | grep -v grep | wc -l` -ge 8 ]] ; do sleep 10s; done }
done
echo

[ -n "$VERBOSE" ] && echo "Waiting for pids: $pids"
[ -z "$TEST"    ] && { printf "All Jobs have started ... waiting for them to finish ... "; wait $pids; printf "and DONE!\n"; }
echo

pids=
for x in $*; do 
    tempFiles=$(for y in elel mumu elmu muel; do echo `dirname $x`/${y}_`basename $x`; done)
    tempFiles=`echo $tempFiles`
    [ -z "$TEST"    ] && echo "Beginning hadd of $OUTPUTDIR/latino_${x/tree_/}"
    [ -n "$VERBOSE" ] && echo " with files: " `for x in $tempFiles; do echo -e "    $x" ; done`
    [ -z "$TEST"    ] && hadd -f $OUTPUTDIR/latino_${x/tree_/} $tempFiles > /dev/null 2>&1 &
    [ -z "$TEST"    ] && pids="$pids $!"
    [ -z "$TEST"    ] && { while [[ `ps | grep hadd | grep -v grep | wc -l` -ge 8 ]] ; do sleep 10s; done }
done
echo

[ -n "$VERBOSE" ] && echo "Waiting for pids: $pids"
[ -z "$TEST"    ] && { printf "All Jobs have started ... waiting for them to finish ... "; wait $pids; printf "and DONE!\n"; }
echo

[ -z "$TEST" ] && printf "Removing temporary files ... "
for x in $*; do 
    tempFiles=$(for y in elel mumu elmu muel; do echo `dirname $x`/${y}_`basename $x`; done)
    tempFiles=`echo $tempFiles`
    [ -n "$VERBOSE" ] && echo "rm $tempFiles"
    [ -z "$TEST"    ] && rm $tempFiles
done
printf "and DONE!\n"

