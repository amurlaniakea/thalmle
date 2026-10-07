#!/bin/sh
# Assembles ~/work-thalmle/thalmle_repo from thalmle_real6 + this documentation bundle.
# Does not modify thalmle_real6. No git commands are run.
set -eu
cd "$HOME/work-thalmle"
SRC=thalmle_real6
DST=thalmle_repo
BUNDLE="$HOME/work-thalmle/repo_docs"
[ -d "$SRC/R" ] || { echo "FALTA $SRC/R"; exit 1; }
[ -f "$BUNDLE/README.md" ] || { echo "FALTA $BUNDLE/README.md"; exit 1; }
rm -rf "$DST"
cp -r "$SRC" "$DST"
mkdir -p "$DST/experiments" "$DST/results"
for f in check_*.R diag*.R paired*.R inspect_g*.R time_*.R smoke_api.R; do
  if [ -f "$f" ]; then cp "$f" "$DST/experiments/"; fi
done
for f in *.log *.rds; do
  if [ -f "$f" ]; then cp "$f" "$DST/results/"; fi
done
cp "$BUNDLE/quickstart.R" "$DST/experiments/quickstart.R"
for f in README.md LICENSE DESCRIPTION NAMESPACE NEWS.md .Rbuildignore .gitignore; do
  cp "$BUNDLE/$f" "$DST/$f"
done
cp "$BUNDLE/test-core.R" "$DST/tests/testthat/test-core.R"
Rscript -e 'sessionInfo()' > "$DST/results/sessionInfo.txt"
Rscript -e 'for (p in c("glmnet","hal9001","haldensify")) cat(p, as.character(packageVersion(p)), "\n")' > "$DST/results/package_versions.txt"
echo "REPO_ENSAMBLADO en $DST"
ls -a "$DST"
