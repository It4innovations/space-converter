# Developer documentation

User documentation (installation, command line, BSPACE add-on, local and remote workflow) is in the [project wiki](https://github.com/It4innovations/space-converter/wiki). The files here are for developers.

| File | Content |
|---|---|
| [SpaceConverter_Dev_Notes_Release1.md](SpaceConverter_Dev_Notes_Release1.md) | Specification of the Release1 code base: repository layout, build, data model, CLI, conversion pipeline, MPI, TCP protocol, readers and the BSPACE add-on. |
| [SpaceConverter_Dev_Notes_Release2.md](SpaceConverter_Dev_Notes_Release2.md) | What changed from Release1 to the current development branch: voxel managers, particle cache, CUDA/HIP, new readers (iPIC3D, PLUTO, RAMSES, BHAC, FIL, FIL_GRACE), CLI changes, add-on changes. |
| [SpaceConverter_Code_Analysis_2026-08.md](SpaceConverter_Code_Analysis_2026-08.md) | Code review and bug hunt of August 2026, with the fix plan and its implementation status. |

The code analysis is a dated report: file and line references point to the code as it was then, the GPU test plan it refers to (`GPU_Test_Plan.md`) is no longer part of the repository, and the `build_spaceconverter_*.sh`, `run_spaceconverter_*.sh` and `gpu_bind_*.sh` scripts it mentions are the developers' cluster scripts, kept outside this repository. The per-cluster scripts shipped with the code are in [scripts/space_converter/](../scripts/space_converter/).
