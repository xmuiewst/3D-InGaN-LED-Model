# Model-Constrained 3D CL Reconstruction Open-Source Package

This folder contains the reviewer-facing reproduction code and data for the 5 kV example workflow associated with:

**Model-Constrained Nanoscale 3D Cathodoluminescence Reconstruction of InGaN LEDs via Coupled Electron–Photon Transport**

## Requirements

Use MATLAB with the Image Processing Toolbox. This package has been smoke-tested with MATLAB R2025b.

Required:

- MATLAB
- Image Processing Toolbox

Optional:

- Parallel Computing Toolbox, used only when MATLAB detects a supported GPU and enables `gpuArray` acceleration.

The open-source 5 kV reproduction entry point does not require a GPU. For full-scale simulations matching the manuscript hardware setting, a workstation with two NVIDIA Tesla V100 GPUs is recommended.

Run all commands from this folder:

```matlab
cd Three_Dim_Model_OpenSource
```

## Self-Tests

Run the lightweight internal tests before running the full reconstruction:

```matlab
addpath(pwd)
ok1 = precies.pixelwise_inversion('selftest');
result2 = precies.structural_inversion('selftest');
assert(ok1)
assert(isstruct(result2))
```

From a system shell, the same check can be run non-interactively:

```bash
matlab -batch "addpath(pwd); ok1=precies.pixelwise_inversion('selftest'); result2=precies.structural_inversion('selftest'); assert(ok1); assert(isstruct(result2));"
```

## Reproduce the 5 kV Reconstruction

Run:

```matlab
run_reconstruction
```

The workflow automatically:

1. Finds `data/5KV-Example.tif`.
2. Uses the full TIFF pixel grid as the ROI.
3. Uses 5 kV electron energy.
4. Simulates every selected pixel with the GUI-default parameters: 0.3 s integration time, 1000 rays, 100 pA beam current, 100 nm depth step, 10 MQW wells, 8 nm barriers, and 5 nm wells.
5. Prints simulation progress as completed pixels over total pixels.
6. Saves simulated mapping TIFF files and the reconstruction metrics.
7. Exports `output/activity_structure_package.mat`.
8. Asks whether to run the forward ablation experiment.

The main outputs are written to `output/`:

| File | Content |
| --- | --- |
| `5KV_simulation_mapping.tif` | Spectral TIFF stack for the simulated mapping |
| `5KV_simulation_mapping_integrated.tif` | Band-integrated simulated mapping |
| `5KV_metrics.csv` | NCC, nRMSE, peak shift, FWHM shift, and intensity-ratio metrics |
| `Table1_validation_metrics.csv` | Table 1-style validation metrics for the 5 kV example |
| `activity_structure_package.mat` | Activity-Structure package equivalent to the ActStruct MAT export |
| `5KV_reconstruction_result.mat` | Full reproduction result bundle |
| `Table2_forward_ablation_metrics.csv` | Created only if forward ablation is run |

For a non-interactive run without the final ablation prompt:

```matlab
run_reconstruction(struct('promptForAblation', false))
```

For a non-interactive run that performs ablation:

```matlab
run_reconstruction(struct('promptForAblation', false, 'runAblation', true))
```

## 3D Modeling

Run:

```matlab
run_3d_modeling
```

Select either `output/activity_structure_package.mat` after running reconstruction or `results/three_dim_model_package.mat`.

This opens two separate figures:

| Figure | Method | Controls |
| --- | --- | --- |
| Legacy Grid Slice Model | Legacy Grid 2D slice reconstruction | Color map, slice position, and X/Y/Z slice direction |
| Activity-Structure 3D Model | Activity-Structure constrained 3D volume | Color map and alpha |

You can also pass the MAT file directly:

```matlab
run_3d_modeling(fullfile('output', 'activity_structure_package.mat'))
```

## Static Paper Results

The `results/` folder contains the static result data referenced in the manuscript.

For GitHub distribution, this folder is intended to be provided as a Release asset named `results.zip`, because several static result files are too large for normal GitHub tracking. Download `results.zip` from the repository release page, then extract it into the repository root so the following paths exist:

| Folder | Content |
| --- | --- |
| `results/Fig7` | Experimental/simulated TIFF data and Fig. 7 comparison outputs |
| `results/Fig8` | Depth-resolved mapping panels |
| `results/Fig9` | 3D visualization outputs |
| `results/Table_1&2` | Cleaned manuscript Table 1 and Table 2 data only |
| `results/three_dim_model_package.mat` | Activity-Structure data package used for manuscript 3D modeling |

`results/Table_1&2` now contains only:

| File | Content |
| --- | --- |
| `Table1_validation_metrics.csv` | Manuscript Table 1 validation metrics |
| `Table2_forward_ablation_metrics.csv` | Manuscript Table 2 M0-M3 ablation metrics |
| `Table1_Table2.xlsx` | Workbook containing the two manuscript tables |

## Citation

If you use this code or the accompanying data, please cite:

Model-Constrained Nanoscale 3D Cathodoluminescence Reconstruction of InGaN LEDs via Coupled Electron–Photon Transport.

The formal journal citation, DOI, and BibTeX entry should be added here after publication.

## License

Before public release, add a standard open-source license file named `LICENSE` at the repository root. For academic reproduction code, BSD-3-Clause or MIT are both common permissive choices. If the `results/` data are distributed separately through GitHub Releases, also state the data reuse terms in the release notes or a `DATA_LICENSE` file.
