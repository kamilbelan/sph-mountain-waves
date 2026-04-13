#!/usr/bin/env bash

set -euo pipefail

REMOTE_HOST="snehurka"
REMOTE_PATH="/usr/users/belank/work/projects/sph-mountain-waves"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MOUNT_POINT="$SCRIPT_DIR/remote"

  # colors
  GREEN='\033[0;32m'
  RED='\033[0;31m'
  YELLOW='\033[1;33m'
  NC='\033[0m'

  is_mounted() {
	  mountpoint -q "$MOUNT_POINT"
  }

do_mount() {
	if is_mounted; then
		echo -e "${YELLOW}Already mounted.${NC}"
		return 0
	fi
	echo "Mounting $REMOTE_HOST:$REMOTE_PATH -> remote/ ..."
	sshfs "$REMOTE_HOST:$REMOTE_PATH" "$MOUNT_POINT" \
		-o reconnect,ServerAliveInterval=15,ServerAliveCountMax=3,follow_symlinks
	echo -e "${GREEN}Mounted.${NC}"
}

do_umount() {
	if ! is_mounted; then
		echo -e "${YELLOW}Not mounted.${NC}"
		return 0
	fi
	echo "Unmounting remote/ ..."
	fusermount -u "$MOUNT_POINT" 2>/dev/null || fusermount -uz "$MOUNT_POINT"
	echo -e "${GREEN}Unmounted.${NC}"
}

do_status() {
	if is_mounted; then
		echo -e "${GREEN}Mounted${NC}: $REMOTE_HOST:$REMOTE_PATH -> remote/"
		df -h "$MOUNT_POINT" | tail -1 | awk '{printf "  Remote disk: %s used of %s (%s free)\n", $3,
		$2, $4}'
	else
		echo -e "${RED}Not mounted.${NC}"
	fi
}

do_reconnect() {
	echo "Reconnecting..."
	do_umount
	do_mount
}

do_sync() {
      local dir="${1:-}"
      local dst="${2:-}"
      if [[ -z "$dir" ]]; then
          echo "Usage: $0 sync <path-under-data/sims/> [custom-destination]"
          exit 1
      fi

      local src="$REMOTE_HOST:$REMOTE_PATH/data/sims/$dir/"
      if [[ -z "$dst" ]]; then
          dst="$SCRIPT_DIR/data/sims/$dir/"
      else
          dst="${dst/#\~/$HOME}"
      fi

      mkdir -p "$dst"
      echo "Syncing $dir -> $dst ..."
      rsync -avz --progress "$src" "$dst"
      echo -e "${GREEN}Done${NC}"
  }

case "${1:-}" in
	mount)     do_mount ;;
	umount)    do_umount ;;
	status)    do_status ;;
	reconnect) do_reconnect ;;
	sync)      do_sync "${2:-}" "${3:-}" ;;
	*)
		echo "Usage: $0 {mount|umount|status|reconnect|sync}"
		exit 1
		;;
esac
