{
  description = "SPH Mountain Waves - A Julia framework for smoothed particle hydrodynamics and its applications in meteorology";

  # The exact nixpkgs revision is pinned in flake.lock (committed). Update it
  # deliberately: `nix flake update`, then `Pkg.update()`, rerun a reference
  # simulation, and commit flake.lock + Manifest.toml together.
  inputs.nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
  let
    system = "x86_64-linux";
    pkgs = import nixpkgs { inherit system; };
    # Julia packages are pinned by Manifest.toml, not by Nix. Nix only supplies
    # the Julia binary, whose minor version must match the Manifest's.
    juliaMinor = pkgs.lib.versions.majorMinor pkgs.julia-bin.version;
  in {
    devShells.${system}.default = pkgs.mkShell {
      buildInputs = with pkgs; [
        julia-bin
        git
      ];

      shellHook = ''
        export JULIA_NUM_THREADS=auto
        export JULIA_PROJECT=@.
        export JULIA_DEPOT_PATH="$PWD/.julia-depot:$HOME/.julia"

        manifestMinor=$(sed -n 's/^julia_version = "\([0-9]*\.[0-9]*\).*/\1/p' Manifest.toml)
        if [ "$manifestMinor" != "${juliaMinor}" ]; then
          echo "WARNING: julia ${pkgs.julia-bin.version} from flake.lock, but Manifest.toml was resolved with $manifestMinor"
        fi

        echo "SPH environment active — Julia ${pkgs.julia-bin.version}, threads: $JULIA_NUM_THREADS"
      '';
    };
  };
}
