{
  description = "Utility classes to check the (robust) equilibrium of a system in contact with the environment.";

  inputs = {
    gepetto.url = "github:gepetto/nix";
    flake-parts.follows = "gepetto/flake-parts";
    nixpkgs.follows = "gepetto/nixpkgs";
    nix-ros-overlay.follows = "gepetto/nix-ros-overlay";
    systems.follows = "gepetto/systems";
    treefmt-nix.follows = "gepetto/treefmt-nix";
  };

  outputs =
    inputs:
    inputs.flake-parts.lib.mkFlake { inherit inputs; } (
      { lib, self, ... }:
      {
        systems = import inputs.systems;
        imports = [
          inputs.gepetto.flakeModule
          { gepetto-pkgs.overlays = [ self.overlays.default ]; }
        ];
        flake.overlays.default = _final: prev: {
          hpp-centroidal-dynamics = prev.hpp-centroidal-dynamics.overrideAttrs {
            src = lib.fileset.toSource {
              root = ./.;
              fileset = lib.fileset.unions [
                ./CMakeLists.txt
                ./include
                ./package.xml
                ./python
                ./src
                ./test
                ./test_data
              ];
            };
          };
        };
        perSystem =
          { pkgs, self', ... }:
          {
            packages = {
              default = self'.packages.hpp-centroidal-dynamics;
              hpp-centroidal-dynamics = pkgs.python3Packages.hpp-centroidal-dynamics;
            };
          };
      }
    );
}
