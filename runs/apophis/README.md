# Apophis DEM run templates

Input files for the Apophis SSDEM runs. These are **templates** — copy them into a
fresh scratch run directory, never run phantom in this folder. (phantom rewrites
`dumpfile=` in the `.in` after every dump, so a live run dir produces constant
spurious git diffs, and re-running would resume the old dump instead of starting
the new setup.)

## Files

| File              | What it is                                                    |
|-------------------|---------------------------------------------------------------|
| `apophis.setup`   | Roche-limit sweep template (Earth present, no spin, no cohesion) |
| `roche.setup`     | Roche variant                                                  |
| `apophis.shape`   | Points phantom at the mesh + scale: `mesh apophis_1.obj 0.170` |
| `itokawa.shape`   | Analytic ellipsoid alternative                                 |
| `apophis_1.obj`   | Apophis shape mesh — **required** by `apophis.shape`            |
| `Makefile`        | Rebuild helper; expects phantom checked out as a sibling dir    |

## Starting a run (same on macOS and WSL/Linux)

    ../../scripts_local/new_run.sh myrun          # or do it by hand:

    mkdir -p ~/runs/roche_rho2.2_N10000
    cd       ~/runs/roche_rho2.2_N10000
    cp ~/phantom/runs/apophis/{apophis.setup,apophis.shape,apophis_1.obj,Makefile} .
    ln -s ~/phantom ../phantom          # Makefile expects PHANTOMDIR=../phantom/
    make setup                          # builds + copies phantomsetup here
    ./phantomsetup apophis              # reads apophis.setup, queries JPL Horizons
    make                                # builds + copies phantom here
    ./phantom apophis.in

Name each run directory uniquely (density, N, machine) so runs from the laptop
and the office PC can never collide when they end up in the same Drive folder.

## Gotchas (learned the hard way)

- **Rebuild before trusting a run.** A stale `phantomsetup` silently ignores
  unknown `.setup` options — no error, no warning. Verify with:
  `strings phantomsetup | grep scale_rho`
- **`phantomsetup` needs internet** — it queries the JPL Horizons API for the
  ephemeris at the given epoch. Setup will fail on an offline/firewalled machine.
- **Re-run `phantomsetup` before each new run**, otherwise `./phantom X.in`
  resumes the previous dump rather than starting your new setup.
- `hfact = 1.000` is correct here (M_6 quintic kernel, hfact_default = 1.0).
- Use the energy-based unbound-mass analysis, not FoF clustering — FoF flips
  between 0 and exactly 50% on intact bodies.
