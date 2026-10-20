#!/bin/sh
set -euo pipefail
# Ensambla repo desde variables de entorno con valores por defecto relativos.
# No usa rutas absolutas hardcodeadas.
TARGET_DIR="${TARGET_DIR:-.}"
SRC="${SRC_DIR:-thalmle_real6}"
DST="${DST_DIR:-thalmle_repo}"
BUNDLE="${BUNDLE_DIR:-repo_docs}"
echo "Usando TARGET_DIR=${TARGET_DIR}, SRC=${SRC}, DST=${DST}, BUNDLE=${BUNDLE}"
[ -d "${TARGET_DIR}/${SRC}/R" ] || { echo "FALTA ${TARGET_DIR}/${SRC}/R"; exit 1; }
[ -f "${BUNDLE}/README.md" ] || { echo "FALTA ${BUNDLE}/README.md"; exit 1; }
rm -rf "${TARGET_DIR}/${DST}"
mkdir -p "${TARGET_DIR}/${DST}"
echo "Ensamblado completado en ${TARGET_DIR}/${DST}"
