#!/bin/bash
# Produce the scheduler/filesystem evidence for the stage 31 job claims.
#
# Committed BESIDE its output so the transcript's provenance claims can be
# checked rather than taken on trust.
#
#   bash capture_scheduler_evidence.sh > 2026-09-22_scheduler_check.txt
#   echo $?        # must be 0; a non-zero exit means the capture is incomplete
#
# What this script does NOT do, stated because earlier versions claimed
# otherwise:
#   * it does not filter the stream. An earlier version piped the ssh output
#     through `grep -v` to drop SSH banner lines and still called the result
#     "raw stdout+stderr". Whatever the remote end emits is what lands here,
#     banners included.
#   * `run()` guarantees that WITHIN ONE EXECUTION the printed command and the
#     evaluated command are the same string. It cannot prevent the transcript
#     being edited afterwards; git records any later byte change.
# Failure propagation, which is what makes a completed transcript meaningful:
#
#   local  : set -e makes a non-zero ssh abort this script.
#   remote : run() checks each command's status and exits 70 on the first
#            failure, printing WHICH command failed into the transcript. That is
#            deliberately not `set -e` on the remote side -- set -e would abort
#            before the handler could name the command, leaving a transcript that
#            simply stops.
#
# An earlier comment here credited `pipefail` with catching a failed ssh whose
# pipeline ended in a successful command. That pipeline was the `grep -v` filter,
# which no longer exists; the comment described a version of the script that was
# already gone. pipefail is kept because it costs nothing, but it is not what is
# load-bearing now.
set -euo pipefail
CLUSTER=xinruqiu@cluster.hpcc.ucr.edu
REPO_ROOT=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
RESULTS=/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2
WINDOW=2026-09-22T00:00:00

cat <<HDR
================================================================================
Captured query output supporting the stage 31 job claims.

Produced by run_evidence/capture_scheduler_evidence.sh, committed alongside this
file. Within one execution, each command is printed by the same string that is
then evaluated -- see the run() echoed below. That is a property of the display
mechanism; it does not certify that this file was never edited afterwards, for
which git history is the record.

The stream is unfiltered: this is the combined remote stdout+stderr exactly as
it arrived, SSH banner lines included. No path is abbreviated, nothing is piped
through head/tail, and sacct uses --parsable2 so no field is cut to a column
width.
================================================================================
cluster                      : ${CLUSTER#*@} (UCR HPCC, SLURM)
account                      : ${CLUSTER%@*}
captured from                : $(hostname) over ssh, non-interactive login shell
captured at                  : $(date '+%Y-%m-%d %H:%M:%S %Z %z')  (capture host clock)
capture-script local repo HEAD: $(git rev-parse HEAD 2>/dev/null || echo 'not a git checkout')
window                       : sacct --starttime=$WINDOW

NOTE: the line above is the LOCAL checkout the script ran from. The HPCC
checkout being inspected is reported separately, by the first remote command.
HDR

ssh -x "$CLUSTER" 'bash -l -s' <<REMOTE 2>&1
set -uo pipefail
run () {
  printf '\n\$ %s\n' "\$1"
  eval "\$1" 2>&1
  rc=\$?
  if [ \$rc -ne 0 ]; then
    printf '\n*** COMMAND FAILED (exit %s): %s\n' "\$rc" "\$1"
    printf '*** capture ABORTED -- this transcript is incomplete\n'
    exit 70
  fi
}
echo ""
echo "the run() used for every command below, as the shell itself reports it:"
declare -f run | sed 's/^/  /'
# `declare -f` prints the actual definition. An earlier version hand-wrote a
# one-line paraphrase under the word "verbatim" -- it dropped the failed
# command's text and the second abort line, and would have drifted again on the
# next edit.
echo ""
echo "cluster clock at query time: \$(date '+%Y-%m-%d %H:%M:%S %Z %z')"
echo "slurm version: \$(sinfo --version 2>&1)"

run "git -C $REPO_ROOT rev-parse HEAD"
run "git -C $REPO_ROOT status --short"
run "sacct --user=${CLUSTER%@*} --starttime=$WINDOW --allocations --parsable2 --format=JobID,JobName,Partition,Submit,Start,End,Elapsed,State,ExitCode,WorkDir"
run "squeue --user=${CLUSTER%@*} --format='%i|%j|%P|%T|%M|%V'"
run "sacct --jobs=29010326 --parsable2 --format=JobID,JobName,Submit,Start,End,Elapsed,State,ExitCode,MaxRSS,WorkDir,SubmitLine"
run "ls -la $RESULTS/"
run "ls -la $RESULTS/contamination_audit/"
run "ls -la $RESULTS/contamination_audit_v2_fixed_semantics_2026-09-22_v1/"
run "find $REPO_ROOT/scripts/31_donor_level_revision/logs -maxdepth 1 -type f -name 'doublet_detect_*' -print"
run "ls -1 $REPO_ROOT/scripts/31_donor_level_revision/logs/"
REMOTE
