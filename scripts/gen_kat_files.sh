#!/bin/sh
set -eu

repo_dir=$(CDPATH='' cd -- "$(dirname -- "$0")/.." && pwd)
if [ "$#" -gt 0 ]; then
    build_dir=$(CDPATH='' cd -- "$1" && pwd)
    shift
else
    build_dir="$repo_dir/build"
fi

kat_dir=$(mktemp -d)
trap 'rm -rf -- "$kat_dir"' EXIT
for level in lvl1 lvl3 lvl5; do
    (cd "$kat_dir" && "$build_dir/apps/PQCgenKAT_sign_$level" "$@")
done
mv "$kat_dir"/PQCsignKAT_*.req "$kat_dir"/PQCsignKAT_*.rsp "$repo_dir/KAT/"
