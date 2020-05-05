# TriWild: Robust Triangulation With Curve Constraints

![](figures/teaser_row.jpg)
Yixin Hu, Teseo Schneider, Xifeng Gao, Qingnan Zhou, Alec Jacobson, Denis Zorin, Daniele Panozzo.
ACM Transactions on Graphics (SIGGRAPH 2019).

[![Build Status](https://travis-ci.org/wildmeshing/TriWild.svg?branch=master)](https://travis-ci.org/wildmeshing/TriWild)
[![Build status](https://ci.appveyor.com/api/projects/status/3k8lru312sw46hs5/branch/master?svg=true)](https://ci.appveyor.com/project/YixinHu69838/triwild/branch/master)
[![Build Status](https://img.shields.io/docker/cloud/build/yixinhu/triwild.svg)](https://hub.docker.com/r/yixinhu/triwild)

## Important Tips

💡💡💡 We also have 3D version of "TriWild" - **TetWild**! It's the parent of TriWild. TetWild can generate linear tetrahedral meshes robustly and automatically. Check it out 👉 **[TetWild](https://github.com/Yixin-Hu/TetWild)**.

💡💡💡 If you are interested in the algorithm details, please refer to our **[paper](https://dl.acm.org/doi/pdf/10.1145/3306346.3323011)** first. We provide plenty of examples and statistics in the paper.

```
@article{Hu:2019:TRT:3306346.3323011,
 author = {Hu, Yixin and Schneider, Teseo and Gao, Xifeng and Zhou, Qingnan and Jacobson, Alec and Zorin, Denis and Panozzo, Daniele},
 title = {TriWild: Robust Triangulation with Curve Constraints},
 journal = {ACM Trans. Graph.},
 issue_date = {July 2019},
 volume = {38},
 number = {4},
 month = jul,
 year = {2019},
 issn = {0730-0301},
 pages = {52:1--52:15},
 articleno = {52},
 numpages = {15},
 url = {http://doi.acm.org/10.1145/3306346.3323011},
 doi = {10.1145/3306346.3323011},
 acmid = {3323011},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {curved triangulation, mesh generation, robust geometry processing},
} 
```

💡💡💡 Check our **[license](https://github.com/wildmeshing/TriWild#license)** first.

## Dataset

💡💡💡 **Please kindly cite our paper when using our pre-generated data.**

### Examples in the Paper

Download [zip](https://drive.google.com/file/d/13xZqYpBz1cV1JaakgkcSO6hSbV9or5V4/view?usp=sharing).

💡💡💡Quickly try TriWild on some small exmaples here!!

### 20k Openclip Dataset

Input: [19686 meshes (.obj) each with a curved feature file (.json)](https://drive.google.com/file/d/1yhrJtfCVMahwgPc0pmX0D8GAJgZ9M7v9/view?usp=sharing)

(For your reference, [here](https://drive.google.com/open?id=1RWzbLKqXeWIQeNYaAfzE8BEcqJZ0R6b6) is original 20k SVG images. Those with animation are not converted to obj/json.)

Output with curved constrains: [19685 meshes (.msh)](https://drive.google.com/open?id=189OP5v5EJNP9QMqpWw_XuGRK_MjMThuJ)

Output with linear constrains(todo James): [19686 meshes (.msh)]()

## Installation

You can use TriWild either by pulling a Docker image or compiling the source code with CMake.

### via Docker

Install Docker and run Docker. Pull TetWild Docker image and run the binary:

```bash
docker pull yixinhu/triwild
docker run --rm -v "$(pwd)":/data yixinhu/triwild /app/TriWild/build/TriWild [TriWild arguments]
```

### via CMake
Our code was originally developed on MacOS and has been tested on Linux and Windows. We provide the commands for installing TriWild in Unix OS: 

- Clone the repository into your local machine:

```bash
git clone https://github.com/wildmeshing/TriWild
```
- Compile the code using cmake (default in Release mode):

```bash
cd TriWild
mkdir build
cd build
cmake ..
make -j
```

- Check the installation:

```bash
./TriWild --help
```
This command should show a list of TriWild parameters.

## Usage

**Input**：

- Linear constraints (required): segment soup in `.obj` format.

- Curved constraints: Bezier curves in `.json` format.

**Output**: Linear/high-order triangle mesh in `.msh` format.

Please check dataset above for examples.

### Quick Try

You can try TriWIld quickly with default parameters by running

```
./TriWild --input input.obj
```
for linear constrains, or

```
./TriWild --input input.obj --feature-input input.json
```
for curved constrains.

### Command Line Switches

```
Usage: ./TriWild [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --input TEXT (REQUIRED)     Input segments in .obj format.
  --output TEXT               Output path.
  --postfix TEXT              Add postfix into outputs' file name.
  --feature-input TEXT        Input feature json file.
  --stop-quality FLOAT        Specify max AMIPS energy for stopping mesh optimization.
  --max-its INT               Max number of mesh optimization iterations.
  --stage INT                 Specify envelope stage
  --envelope-r FLOAT          relative envelope epsilon_r. Absolute epsilonn = epsilon_r * diagonal_of_bbox
  --feature-envelope-r FLOAT  Relative feature envelope mu_r. Absolute mu = mu_r * diagonal_of_bbox
  --target-edge-length FLOAT  Absolute target edge length l.
  --target-edge-length-r FLOAT
                              Relative target edge length l_r. Absolute l = l_r * diagonal_of_bbox
  --log-file TEXT             Output a log file.
  --min-angle FLOAT           Desired minimal angle.
  --mute-log                  Mute prints.
  --cut-outside               Remove "outside part".
  --skip-eps                  Skip saving eps.
  --cut-holes TEXT            Input a .xyz file for specifying points inside holes you want to remove.
  --output-linear-mesh        Output linear mesh for curved pipeline.
```

More details about some important parameters:

* **`--feature-input`**

We provide a [python script](https://github.com/teseoch/svg2obj) for converting a svg to curves in `.json` format.

* **`--envelope`**

Relative surface envelope <a href="https://www.codecogs.com/eqnedit.php?latex=$\epsilon_r$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\epsilon_r$" title="$\epsilon_r$" /></a> (1e-3 in default). Absolute surface envelope <a href="https://www.codecogs.com/eqnedit.php?latex=$\epsilon&space;=&space;\epsilon_r&space;d$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\epsilon&space;=&space;\epsilon_r&space;d$" title="$\epsilon = \epsilon_r d$" /></a>, where <a href="https://www.codecogs.com/eqnedit.php?latex=$d$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$d$" title="$d$" /></a> is the length of the diagonal of the bounding box of input.

* **`--feature-envelope`**

Relative feature envelope <a href="https://www.codecogs.com/eqnedit.php?latex=$\mu_r$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\mu_r$" title="$\mu_r$" /></a> (1e-3 in default with linear constraints and 2e-3 for curved constraints). Absolute feature envelope <a href="https://www.codecogs.com/eqnedit.php?latex=$\mu&space;=&space;\mu_r&space;d$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\mu&space;=&space;\mu_r&space;d$" title="$\mu = \mu_r d$" /></a>.

* **`--target-edge-length-r`**

Relative targeted edge length <a href="https://www.codecogs.com/eqnedit.php?latex=$\ell_r$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\ell_r$" title="$\ell_r$" /></a> (0.05 in default). Absolute targeted edge length <a href="https://www.codecogs.com/eqnedit.php?latex=$\ell&space;=&space;\ell_r&space;d$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\ell&space;=&space;\ell_r&space;d$" title="$\ell = \ell_r d$" /></a>.

## License

TriWild is MPL2 licensed and free for both commercial and non-commercial usage. However, you have to cite our work in your paper or put a reference of TriWild in your software. Whenever you fix bugs or make some improvement of TriWild, you should contribute back.

## Gallery
![](figures/mosaic_new.png)

<!--## Acknowledgements

(todo)-->
