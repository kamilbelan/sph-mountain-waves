{
  description = "SPH Mountain Waves - A Julia framework for smoothed particle hydrodynamics and its applications in meteorology";

  inputs.nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
  let
    system = "x86_64-linux";
    pkgs = import nixpkgs { inherit system; };
    gccLibPath = "${pkgs.gcc-unwrapped.lib}/lib";
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
        export LD_LIBRARY_PATH=${gccLibPath}:$LD_LIBRARY_PATH
        export LD_PRELOAD=${gccLibPath}/libquadmath.so.0
        echo "SPH environment active — Julia threads: $JULIA_NUM_THREADS"
      '';
    };
  };
}
