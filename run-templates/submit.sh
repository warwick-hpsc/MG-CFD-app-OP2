echo "Processing job <JOB_NUM>/<NUM_JOBS>"
if [ -f <RUN_DIRPATH>/job-is-complete.txt ]; then
  echo " - already complete"
else
  basedir=`pwd`
  cd <RUN_DIRPATH>

  if [ "$submit_cmd" = "" ]; then
    ./<BATCH_FILENAME> 2>&1 | tee run.log
  else
    if [ -f "job-in-queue.txt" ]; then
      echo " - already in queue"
    elif [ -f "job-is-running.txt" ]; then
      echo " - already running"
    else
      if [[ `hostname` == *"login"* ]]; then
        if [ ! -f <BIN_FILEPATH> ]; then
          ## Run batch script on login node to perform compilation, which will exit before app execution.
          ./<BATCH_FILENAME>
        fi
        if ! eval "$submit_cmd" ./<BATCH_FILENAME> ; then
          echo "WARNING: Submission failed for: <RUN_DIRPATH>/<BATCH_FILENAME>"
        else
          touch <RUN_DIRPATH>/"job-in-queue.txt"
        fi
      else
        if ! eval "$submit_cmd" ./<BATCH_FILENAME> ; then
          echo "WARNING: Submission failed for: <RUN_DIRPATH>/<BATCH_FILENAME>"
        fi
      fi
    fi
  fi
  cd "$basedir"
fi
