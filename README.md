# BU-RSI

BU-RSI is a low-order rotor-stator broadband interaction noise prediction tool. It includes an annular acoustic frequency-selection utility, the RSI solver build workflow, example input and output data, and Python scripts for running multiple cases.

## Repository layout

| Path | Purpose |
| --- | --- |
| `findfreqs_src/` | Annular eigenfrequency calculation and MATLAB frequency selection |
| `rsi_src/` | RSI solver build files and per-case batch runner |
| `libsk3/nasalib/` | Included NASA numerical routines and static library |
| `libsk3/linpacklib/` | Included LINPACK routines and static library |
| `examples/rsi_input/ADP_2022/` | Example RSI namelist inputs |
| `examples/rsi_output/ADP_2022/` | Example processed spectra |
| `runrsi_batch.py` | Top-level multi-case workflow |
| `get_output.py` | Collects and sorts `fort.33` results from completed runs |

The original manual refers to the `FINDFREQS` and `CODE` directories. In this repository, their corresponding directories are `findfreqs_src/` and `rsi_src/`.

## Requirements

The exact requirements depend on whether you use the included binaries and libraries or rebuild everything.

- Linux or a compatible HPC environment
- A Fortran toolchain compatible with the supplied makefiles
  - `pgf77` for the included numerical libraries
  - `pgf95` for `comp2` and `rsi95`
- `make` and `ar`
- MATLAB for `findfreqs_src/freqs2.m`
- Python 3
- Python packages: `numpy` and `pandas`
- Sun Grid Engine-compatible `qsub` for the supplied parallel batch workflow

## Workflow overview

1. Generate annular eigenvalues with `comp2`.
2. Run `freqs2.m` to select nonsingular acoustic frequencies and write `fs.dat`.
3. Prepare the three RSI namelist files.
4. Place the inputs and `fs.dat` beside the `rsi95` executable.
5. Run `rsi95`, directly or through the supplied batch scripts.
6. Read the spectrum from `fort.33`.

## 1. Select acoustic frequencies

The frequency-check module is in `findfreqs_src/`.

### Build `comp2`

An executable is included. To rebuild it:

```bash
cd findfreqs_src
make -f Makecomp2.pgi
```

The build uses `CompXmn2.f95` and `../libsk3/nasalib/nasalib.a`.

### Run `comp2`

Edit `findfreqs_src/aero95const.input`, then run:

```bash
cd findfreqs_src
./comp2
```

An IEEE warning at the end of the calculation may be ignored, according to the original documentation.

The main output is `fort.102`, a three-column ASCII file containing:

1. circumferential mode number;
2. radial mode number;
3. annular eigenfrequency parameter.

### Run `freqs2.m`

Before running `findfreqs_src/freqs2.m`, set:

- `M_ann`: duct bulk Mach number;
- `Rtip`: duct tip radius in feet;
- `T`: duct temperature in Rankine;
- `nfreq0` through `nfreq8`: desired numbers of frequencies in the frequency bands defined by the script.

Keep `fort.102` in the same directory and run the MATLAB script:

```matlab
run('freqs2.m')
```

The script writes:

- `fs.sing`: singular frequencies for acoustic propagation in the annulus;
- `fs.dat`: selected nonsingular frequencies for the RSI calculation.

Inspect the plotted distribution and edit `fs.dat` when needed to obtain a suitable set of calculation frequencies.

## 2. Build and run RSI

The RSI executable is named `rsi95`. Its makefile is `rsi_src/Makersifull.pgi`:

```bash
cd rsi_src
make -f Makersifull.pgi
```

The build expects these source files:

- `rsi.f95`
- `get_Pmns.f95`
- `casc_k3_no_r.f`

It links against:

- `rsi_src/libsk3/nasalib/nasalib.a`
- `rsi_src/libsk3/linpacklib/linpack.a`

Before building, verify that `rsi.f95` is available. It is referenced by the makefile and batch workflow but is not included in the current repository tree.

For a direct run, place the following files in the working directory with `rsi95`:

- `aero95const.input`
- `aero95mat.input`
- `turb95mat.input`
- `fs.dat`

Then run:

```bash
./rsi95
```

The output spectrum of interest is written to `fort.33`.

## Input reference

The input files use Fortran namelist format. Example cases are available under `examples/rsi_input/ADP_2022/`.

### `aero95const.input`

This file defines scalar geometry, flow, and numerical settings.

| Parameter | Description |
| --- | --- |
| `NV` | Number of stator vanes |
| `NB` | Number of rotor blades |
| `NRdat` | Number of spanwise input positions, including the endwalls |
| `NRcalc` | Number of equally spaced spanwise strips used to calculate vane lift response |
| `NC` | Number of chordwise points used to calculate vane lift response |
| `NF` | Number of frequencies read from `fs.dat` and computed |
| `Ng` | Number of terms in the `k3` discretization |
| `T` | Spanwise mass-averaged static temperature at the vane leading edge, in Rankine |
| `rho` | Spanwise mass-averaged duct density, in `lbm/ft^3` |
| `M_ann` | Duct bulk Mach number, assuming plug flow |
| `iSTDYL` | Incidence-angle correction switch; the documented default is `0` |
| `tipc_to_tipr` | Vane tip chord divided by duct radius at the vane leading edge |
| `hub_to_tip` | Duct hub-to-tip radius ratio at the vane leading edge |
| `sweepdeg` | Vane sweep angle in degrees, measured from radial; aft sweep is positive |
| `RDUCT` | Duct radius at the stator leading edge, in inches |
| `BW` | Spectral bandwidth in Hz |
| `conv_tol_xmn` | Root-finding convergence tolerance for the nondimensional eigenvalue equation |

Additional guidance:

- `NRcalc` is independent of the `NRdat` locations. Input values are interpolated from the `NRdat` positions to the calculation strips.
- Choose `NRcalc` with the radial turbulent length scale in mind; the strip size should be roughly comparable to that scale.
- Use `NC >= 80`. The original guidance recommends approximately `80` for `k_z <= 60` and more than `80` for `k_z > 60`.
- `k3` spans from `0` to `24` in increments of `24/Ng`. Change `dk3` in the solver source to alter this range.
- `T` is documented as being used only for `dB_ref`; the legacy code otherwise sets the speed of sound to `1115 ft/s`.
- `rho` is used only in the calculation of `dB_ref`.
- `iSTDYL = 1` enables the incidence-angle correction in the Mach number normal to the vane leading edge. The original author always used `iSTDYL = 0` and did not validate the alternate setting.
- Compute `hub_to_tip` from the leading-edge radii as `r(1) / r(NRdat)`.
- The bandwidth correction is `BWc = 10 log10(BW)`. For a prediction rather than a match to test-data bandwidth, the original guidance recommends `BW = 0` and applying any correction afterward.
- The example convergence tolerance is `conv_tol_xmn = 0.0001`.

### `aero95mat.input`

The arrays in this file each contain `NRdat` values.

| Parameter | Description |
| --- | --- |
| `radius` | Radial positions normalized by the tip radius |
| `chord` | Vane chord at each radius normalized by the tip chord |
| `stagger` | Vane stagger angle at each radius, in degrees |
| `Mx` | Mean absolute axial Mach number at each radius at the vane leading edge |
| `Mt` | Mean absolute tangential Mach number at each radius at the vane leading edge |


### `turb95mat.input`

The turbulence arrays each contain `NRdat` values and are evaluated at or upstream of the vane leading edge. They may come from fan-stage CFD, hot-wire measurements, or laser velocimetry. If no case-specific data are available, the original documentation recommends using data from a similar fan-stage configuration.

| Parameter | Description |
| --- | --- |
| `backurms` | Average-passage RMS turbulence intensity normalized by the local mean absolute velocity |
| `maxurms` | Legacy Gaussian-fit parameter; not used and should be set to zero |
| `urmsgausswidth` | Legacy Gaussian-fit parameter; not used and should be set to zero |
| `Lxdat` | Longitudinal integral length scale normalized by duct tip radius |
| `Ltdat` | Upwash integral length scale normalized by duct tip radius |
| `Lrdat` | Radial integral length scale normalized by duct tip radius |

The current model uses only the average-passage turbulence intensity in `backurms`. Set `maxurms` and `urmsgausswidth` to zero.

Length-scale directions are defined as follows:

- `Lxdat`: along the vane chord, not the axial direction;
- `Ltdat`: normal to the vane chord, not the tangential direction;
- `Lrdat`: radial direction.

For isotropic turbulence inferred from CFD, the documented assumptions are:

```text
Ltdat = Lxdat / 2
Lrdat = Lxdat / 2
```

## Output format

`fort.33` is a whitespace-delimited, three-column spectrum:

| Column | Quantity |
| --- | --- |
| 1 | Frequency in Hz |
| 2 | Inlet power in dB |
| 3 | Exhaust power in dB |

Example:

```text
2320.45000000000  95.8128864000951  102.416565511673
```

Only the exhaust-power calculation has been validated. The inlet-power column is produced by the code but should not be treated as validated.

The documented power-density conversion in `rsi.f95` is:

```text
sPWR = 10 log10(sPWR) + dB_ref + BWc
```

Here, `dB_ref` dimensionalizes the result and `BWc` accounts for the experimental data-collection bandwidth. The original workflow generally applied the bandwidth adjustment after the RSI run.

## Batch workflow

`runrsi_batch.py` prepares and submits multiple cases. Review and edit its configuration block before running it, especially:

- `parafreq`: number of frequency groups submitted in parallel; `0` selects one sequential solver job;
- `max_runtime`: requested runtime for the sequential job;
- `nfreq`: requested frequency counts across the script's nine frequency bands;
- `nf`: number of frequencies for the sequential path;
- `input_dir`, `findfreq_dir`, `rsi_main_dir`, and `rsi_sub_dir`;
- `maxf`, `ng`, `nc`, and `nrcalc` for parallel runs;
- `input_name`: optional explicit case names.

Create the run and output parent directories if they do not already exist:

```bash
mkdir -p examples/rsi_run/ADP_2022
mkdir -p examples/rsi_output/ADP_2022
```

Install the Python dependencies and start the workflow from the repository root:

```bash
python3 -m pip install numpy pandas
python3 runrsi_batch.py
```

For each case, the script:

1. copies `rsi_src/` into a case-specific run directory;
2. renames the case inputs to the standard RSI filenames;
3. updates `freqs2.m` and generates `fs.dat` with MATLAB;
4. updates numerical settings in the copied input files;
5. builds and submits one or more RSI executables through `qsub`.

After all jobs complete, review the paths and missing-file behavior in `get_output.py`, then collect and sort the results:

```bash
python3 get_output.py
```

The collector writes case-level `.33` files under `examples/rsi_output/ADP_2022/`.
