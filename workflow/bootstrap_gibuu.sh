#!/usr/bin/env bash
set -euo pipefail

install_root=${1:-/sps/juno/yqiyu/GiBUU2025_p5}
env_script=${2:-/sps/juno/yqiyu/GiBUUGEN/env.sh}
release_url='https://gibuu.hepforge.org/downloads?f=release2025.tar.gz'
buuinput_url='https://gibuu.hepforge.org/downloads?f=buuinput2025.tar.gz'

mkdir -p "$install_root"
cd "$install_root"

if [[ ! -s release2025.tar.gz ]]; then
  wget --content-disposition "$release_url"
fi
if [[ ! -s buuinput2025.tar.gz ]]; then
  wget --content-disposition "$buuinput_url"
fi

sha256sum release2025.tar.gz buuinput2025.tar.gz > checksums.sha256

if [[ ! -d release2025 ]]; then
  tar -xzf release2025.tar.gz
fi
if [[ ! -d buuinput2025 ]]; then
  tar -xzf buuinput2025.tar.gz
fi

if [[ ! -e release ]]; then
  ln -s release2025 release
fi
if [[ ! -e buuinput ]]; then
  ln -s buuinput2025 buuinput
fi

source "$env_script"
make -C release -j4

version=$(<release/version.txt)
if [[ "$version" != *"Release 2025, patch 5"* ]]; then
  echo "Expected GiBUU Release 2025 patch 5, got: $version" >&2
  exit 1
fi
if [[ ! -x release/objects/GiBUU.x ]]; then
  echo "GiBUU executable was not created" >&2
  exit 1
fi

{
  printf 'installed_at=%s\n' "$(date --iso-8601=seconds)"
  printf 'version=%s\n' "$version"
  printf 'compiler=%s\n' "$(gfortran --version | head -n 1)"
  sed 's/^/sha256=/' checksums.sha256
} > installation_manifest.txt

printf 'Installed %s in %s\n' "$version" "$install_root"
