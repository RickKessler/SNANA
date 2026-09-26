#!/usr/bin/env bash
#
# Regression test for the D+U*U^T (Woodbury) factorized-covariance patch to
# wfit.c / create_covariance.py (new_branch2 / SNANA_woodbury). Two
# independent, automated checks -- run this before ever proposing the patch
# for merge, and after any future change to either file.
#
#   Check 1 (numerical equivalence): wfit with -mucov_factorized_file must
#     reproduce wfit with -mucovtot_inv_file (dense) on the SAME COVOPT to
#     tight tolerance -- confirms the Woodbury identity is evaluated
#     correctly, not just that it runs without crashing.
#
#   Check 2 (INFO.YML filename correctness -- regression test for the Sep
#     2026 bug fix): create_covariance.py --write_factorized together with
#     --write_format_cov text (i.e. NOT the npz default) must write a
#     covfactorized_NNN.npz file whose name EXACTLY matches what INFO.YML
#     records for that COVOPT. Before the fix, get_cov_filename() named the
#     factorized file using the dense-cov suffix (e.g. "covfactorized_000.
#     txt.gz") while write_covariance_factorized() always calls np.savez()
#     regardless, actually writing "covfactorized_000.txt.npz" -- any
#     consumer (submit_prog_cosmofit.py, or a user script) trusting
#     INFO.YML would get a silent FileNotFoundError at fit time.
#
# Usage: ./test_factorized_woodbury.sh   (exits 0 on pass, 1 on any failure)
#
# Ayan Mitra, Sep 2026

set -u
FAIL=0
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

TESTDIR=$(cd "$(dirname "$0")" && pwd)
if [ -x "$TESTDIR/bin/wfit.exe" ]; then WFIT=$TESTDIR/bin/wfit.exe
elif [ -x "$HOME/SNANA/bin/wfit.exe" ]; then WFIT=$HOME/SNANA/bin/wfit.exe
else echo "Error: cannot find patched wfit.exe" >&2; exit 1; fi

if [ -f "$TESTDIR/util/create_covariance.py" ]; then CC_PY=$TESTDIR/util/create_covariance.py
elif [ -f "$HOME/SNANA_woodbury/util/create_covariance.py" ]; then CC_PY=$HOME/SNANA_woodbury/util/create_covariance.py
else echo "Error: cannot find patched create_covariance.py" >&2; exit 1; fi

WOPTS="-cmb_sim -sigma_Rcmb 0.006 -wsteps 201 -omsteps 101 -ranseed_Rcmb 1"
TOL=1e-4

UNBINNED_RUNDIR=/pscratch/sd/d/desctd/PIPPIN_OUTPUT/LSST_ANALYSIS-2/7_CREATE_COV/LSST_UNBINNED_COV
UNBINNED_BBCDIR=/pscratch/sd/d/desctd/PIPPIN_OUTPUT/LSST_ANALYSIS-2/6_BIASCOR/BBC_SIMDATA_PHOTOZ/output
UNBINNED_OFFICIAL_OUTDIR=$UNBINNED_RUNDIR/output/LSST_UNBINNED_COV_BBC_SIMDATA_PHOTOZ_OUTPUT_BBCFIT-0001

echo "======================================================================"
echo " Check 1: factorized wfit == dense wfit (numerical equivalence)"
echo "======================================================================"
OUT1=$TESTDIR/output_BBCFIT-0001
if [ ! -f "$OUT1/covfactorized_000.npz" ] || [ ! -f "$UNBINNED_OFFICIAL_OUTDIR/covtot_inv_000.npz" ]; then
  echo "  SKIP: prerequisite covfactorized_000.npz / covtot_inv_000.npz not found"
  echo "        (run run_3way_wfit_test.sh or run_create_cov_to_wfit.sh first)"
else
  $WFIT $UNBINNED_OFFICIAL_OUTDIR/hubble_diagram.txt $WOPTS \
        -mucovtot_inv_file $UNBINNED_OFFICIAL_OUTDIR/covtot_inv_000.npz \
        -cospar_yaml $TMPDIR/check1_dense.yaml > $TMPDIR/check1_dense.log 2>&1
  $WFIT $OUT1/hubble_diagram.txt $WOPTS \
        -mucov_factorized_file $OUT1/covfactorized_000.npz \
        -cospar_yaml $TMPDIR/check1_fac.yaml > $TMPDIR/check1_fac.log 2>&1

  python3 - "$TMPDIR/check1_dense.yaml" "$TMPDIR/check1_fac.yaml" "$TOL" <<'PYEOF'
import sys, yaml
dense_file, fac_file, tol = sys.argv[1], sys.argv[2], float(sys.argv[3])
d = yaml.safe_load(open(dense_file))
f = yaml.safe_load(open(fac_file))
ok = True
for key in ("w", "om", "chi2"):
    dv, fv = float(d[key]), float(f[key])
    rel = abs(dv - fv) / max(abs(dv), 1e-12)
    status = "OK" if rel < tol else "FAIL"
    if status == "FAIL":
        ok = False
    print(f"  {key}: dense={dv}  factorized={fv}  rel_diff={rel:.2e}  [{status}]")
sys.exit(0 if ok else 1)
PYEOF
  [ $? -ne 0 ] && FAIL=1
fi

echo
echo "======================================================================"
echo " Check 2: INFO.YML covfactorized filename matches file on disk"
echo "          (regression test for the write_format_cov=text/csv bug fix)"
echo "======================================================================"
OUT2=$TMPDIR/check2
mkdir -p $OUT2
python3 $CC_PY $UNBINNED_RUNDIR/input_file.txt \
     --input_dir $UNBINNED_BBCDIR --version OUTPUT_BBCFIT-0001 \
     --outdir $OUT2 --unbinned --write_factorized \
     --write_mask_cov 0 --write_format_cov text \
     > $OUT2/log.txt 2>&1

RECORDED=$(grep "^  0:" $OUT2/INFO.YML | awk '{print $NF}')
if [ -n "$RECORDED" ] && [ -f "$OUT2/$RECORDED" ]; then
  echo "  OK: INFO.YML records '$RECORDED', and it exists on disk"
else
  echo "  FAIL: INFO.YML records '$RECORDED', but that file does NOT exist on disk"
  echo "        (actual covfactorized_000.* present: $(ls $OUT2 2>/dev/null | grep covfactorized_000))"
  FAIL=1
fi

echo
if [ $FAIL -eq 0 ]; then
  echo "ALL CHECKS PASSED"
else
  echo "SOME CHECKS FAILED"
fi
exit $FAIL
