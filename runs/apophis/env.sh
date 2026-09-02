# Source before building/running:  source env.sh
# Keeps builds reproducible across macOS and WSL/Linux.

export SYSTEM=gfortran

# macOS only: Anaconda's ld (ld64-530) cannot read the modern SDK's tbd-v4 stubs
# and breaks the link step. `conda deactivate` does NOT remove it from PATH.
case "$(uname -s)" in
  Darwin) export PATH="$(echo "$PATH" | tr ':' '\n' | grep -v anaconda | paste -sd: -)" ;;
esac

# OpenMP: phantom needs a large per-thread stack or it segfaults in the force loop.
export OMP_STACKSIZE=512M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$( (nproc 2>/dev/null || sysctl -n hw.physicalcpu) )}
ulimit -s unlimited 2>/dev/null || true

echo "SYSTEM=$SYSTEM  OMP_NUM_THREADS=$OMP_NUM_THREADS  OMP_STACKSIZE=$OMP_STACKSIZE"
