{
  description = "Utility classes to check the (robust) equilibrium of a system in contact with the environment.";

  inputs.gepetto.url = "github:gepetto/nix";

  outputs =
    inputs:
    inputs.gepetto.lib.mkFlakoboros inputs (
      { lib, ... }:
      {
        overrideAttrs.hpp-centroidal-dynamics = {
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
      }
    );
}
