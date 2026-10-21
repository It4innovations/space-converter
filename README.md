# Space Converter

SpaceConverter is an MPI application that turns astrophysical and cosmological simulation output into volumetric OpenVDB/NanoVDB grids for rendering in Blender and Cycles. It reads particle (SPH, N-body, PIC) and grid (uniform, AMR) data from OpenGADGET3, ChaNGa, HACC, iPIC3D, PLUTO, RAMSES, BHAC, the Einstein Toolkit (Carpet HDF5) and GRACE, deposits it onto a sparse or dense grid with selectable SPH kernels, and scales across CPU nodes or NVIDIA (CUDA) and AMD (HIP) GPUs.

It runs either as a batch converter that writes files, or as a TCP server driven interactively from Blender through the **BSPACE** add-on (`blender_addons/bspace`). In remote mode the volumes go straight to the [**CyclesPhi**](https://github.com/It4innovations/cyclesphi) renderer on the cluster, so the data never leaves it, and Blender shows the rendered images through the [BRaaS-HPC RenderEngine add-on](https://github.com/It4innovations/braas-hpc-renderengine-addon) (built on the [BRaaS-HPC RenderEngine](https://github.com/It4innovations/braas-hpc-renderengine) library).

## Wiki

[**Home**](https://github.com/It4innovations/space-converter/wiki)

[**Installation Guide**](https://github.com/It4innovations/space-converter/wiki/Installation-Guide)
- [Linux Installation](https://github.com/It4innovations/space-converter/wiki/Linux-Installation)
- [macOS Installation](https://github.com/It4innovations/space-converter/wiki/MacOS-Installation)
- [Windows Installation](https://github.com/It4innovations/space-converter/wiki/Windows-Installation)

[**Usage Guide**](https://github.com/It4innovations/space-converter/wiki/Usage-Guide)
- [Supported Data Formats](https://github.com/It4innovations/space-converter/wiki/Supported-Data-Formats)
- [How to Use SpaceConverter](https://github.com/It4innovations/space-converter/wiki/How-to-Use-SpaceConverter)
- [How to Use BSPACE](https://github.com/It4innovations/space-converter/wiki/How-to-Use-BSPACE-Add‐on)

[**Application Workflow Overview**](https://github.com/It4innovations/space-converter/wiki/Application-Workflow-Overview)
- [Local Mode](https://github.com/It4innovations/space-converter/wiki/Local-Mode)
- [Remote Mode](https://github.com/It4innovations/space-converter/wiki/Remote-Mode)

Developer documentation is in [docs/](docs/), and the smoke tests are described in [tests/README.md](tests/README.md).

## License

SpaceConverter and the BSPACE add-on are copyright (C) 2023-2026 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava, and are released under the GNU General Public License, version 3 or (at your option) any later version (GPL-3.0-or-later). The full text is in [LICENSE](LICENSE).

Third-party components keep their own licenses:

| Component | Location | License |
|---|---|---|
| OpenVDB / NanoVDB | `submodules/openvdb` | Apache-2.0 |
| c-blosc | `submodules/c-blosc` | BSD-3-Clause |
| nanoflann | `submodules/nanoflann` | BSD-2-Clause |
| cudaKDTree | `submodules/cudaKDTree` | Apache-2.0 |
| braas-hpc-renderengine | `submodules/braas-hpc-renderengine` | GPL-3.0-or-later |
| GenericIO, SZ (HACC reader) | `src/hacc/utility` | UChicago Argonne open source licenses, see `COPYING` and `SZ/COPYRIGHT.txt` there |

Other code bundled under `src/*/utility` (the OpenGADGET3 I/O subset and the ChaNGa/Tipsy structures) comes from the respective simulation codes and remains under the terms of those projects.

## Acknowledgement

This work was supported by the SPACE project under grant agreement No 101093441. The project is supported by the European High-Performance Computing Joint Undertaking and its members (including top-up funding by the Ministry of Education, Youth and Sports of the Czech Republic ID: MC2304).

This work was also supported by the Ministry of Education, Youth and Sports of the Czech Republic through the e-INFRA CZ (ID:90254).

## References

- Jaros, M., Strakos, P., Barbosa, J., Riha, L., and Rizzi, S. (2026). [Visualizing the Cosmos: Bridging Astrophysical Simulations with Volume Rendering](https://doi.org/10.2312/visgap.20261003). VisGap 2026.
- Jaroš, M., Strakoš, P., & Riha, L. (2026). SpaceConverter: Scientific Simulation Visualization Gallery [Video recording]. Zenodo. https://doi.org/10.5281/zenodo.23394314
- Taffoni, G., et al. (2024). [Redesign of astrophysical codes for exascale computing: the SPACE experience](https://doi.org/10.1117/12.3020361). Proc. SPIE 13101, Software and Cyberinfrastructure for Astronomy VIII, 131010X.

## Development

AI agent harnesses were used for selected development, testing, and optimization tasks.
