# ABCDGraphGenerator.jl
Artificial Benchmark for Community Detection (ABCD) - A Fast Random Graph Model with Community Structure

Bogumił Kamiński, Paweł Prałat, François Théberge

### Julia API

The package does not export functions and types. The public API is the following:
* `ABCDGraphGenerator.ABCDParams`: type holding information about sampled degrees,
  sampled cluster sizes and required mode of ABCD graph generation
* `ABCDGraphGenerator.gen_graph`: ABCD graph generator that uses `ABCDParams`
  specification
* `ABCDGraphGenerator.sample_degrees`: sample degrees of vertices following power law
* `ABCDGraphGenerator.sample_communities`: sample community sizes following power law
* `ABCDGraphGenerator.get_ev`: get expected value of truncated discrete power law distribution
* `ABCDGraphGenerator.find_v_min`: find the lower truncation given expected value
  and upper truncation of truncated discrete power law distribution

The reason for such split of the functionality is that generation of vertex degrees
and community sizes is fast, while the generation of the final graph is the most expensive step.

### Using ABCDGraphGenerator.jl from R and Python

The functions provided in the package can be directly called from R and Python.

Instructions how to seamlessly integrate Julia into R session are given in the [JuliaCall](https://cran.r-project.org/web/packages/JuliaCall/index.html) package documentation.

An interface to call Julia directly from Python is provided by the [PyJulia](https://github.com/JuliaPy/pyjulia) package.

### Command Line Interface

The `utils/` folder contains command line interface utilities that are aimed at users that want
to the use package without using the API directly. A requirement for these
utilities to be run is to have the Julia language in version at least 1.0 installed on a computer.
It contains the following files:
* `install.jl`: installs all required packages
* `abcd_sampler.jl`: generates an ABCD graph following a configuration file
* `deg_sampler.jl`: samples degrees of vertices in the graph
* `com_sampler.jl`: samples communitiy sizes in the graph
* `graph_sampler.jl`: samples edges and community assignments in the graph

The main file intended to be used is `abcd_sampler.jl`.
Here is an example configuration file, named `example_config.toml`, in this guide:
```
seed = "42"                   # RNG seed, use "" for no seeding
n = "10000"                   # number of vertices in graph
t1 = "3"                      # power-law exponent for degree distribution
d_min = "5"                   # minimum degree
d_max = "50"                  # maximum degree
d_max_iter = "1000"           # maximum number of iterations for sampling degrees
t2 = "2"                      # power-law exponent for cluster size distribution
c_min = "50"                  # minimum cluster size
c_max = "1000"                # maximum cluster size
c_max_iter = "1000"           # maximum number of iterations for sampling cluster sizes
# Exactly one of xi and mu must be passed as Float64. Also if xi is provided islocal must be set to false or omitted.
xi = "0.2"                    # fraction of edges to fall in background graph
#mu = "0.2"                   # mixing parameter
islocal = "false"             # if "true" mixing parameter is restricted to local cluster, otherwise it is global
isCL = "false"                # if "false" use configuration model, if "true" use Chung-Lu
degreefile = "deg.dat"        # name of file do generate that contains vertex degrees
communitysizesfile = "cs.dat" # name of file do generate that contains community sizes
communityfile = "com.dat"     # name of file do generate that contains assignments of vertices to communities
networkfile = "edge.dat"      # name of file do generate that contains edges of the generated graph
nout = "100"                  # number of vertices in graph that are outliers; optional parameter
                              # if nout is passed and is not zero then we require islocal = "false",
                              # isCL = "false", and xi (not mu) must be passed
                              # if nout > 0 then it is recommended that xi > 0
```
In this file all parameters required to generate an ABCD graph and store to on disk are passed.

Note that if `nout` is passed and is greater than 0 then the first community consists of outliers
and its size is equal to `nout`.

Here is an output from an example session using CLI in the ABCD-generation mode using the above file:
```
$ julia abcd_sampler.jl example_config.toml
[ Info: Usage: julia abcd_sampler.jl config_filename
[ Info: For the syntax of config_filename see example_config.toml file
[ Info: Expected value of degree: 8.327743727955891
[ Info: Expected value of community size: 156.5613820733916
$ shasum -a 256 edge.dat #sha256sum edge.dat on Linux
a2102912153726dc20da30ee7dbfc479ee6b8fa785a11a04e5d044320f50c0fa  edge.dat
```
After the program terminates four files, `deg.dat`, `cs.dat`, `com.dat` and `edge.dat`
are created in the working directory.

`deg_sampler.jl`, `com_sampler.jl` and `graph_sampler.jl` files are provided
mainly to facilitate comparisons with LFR algorithm.
Here is an output from an example session using CLI in the LFR-comparison mode.
SHA values are computed under Julia 1.7.

Note that in this mode if `nout` is passed must be passed to community
generation and graph generation functions and must be the same.

```
$ julia install.jl
  Updating registry at `~\.julia\registries\General`
  Updating git-repo `https://github.com/JuliaRegistries/General.git`
  Updating git-repo `https://github.com/bkamins/ABCDGraphGenerator.jl`
  Updating git-repo `https://github.com/bkamins/ABCDGraphGenerator.jl`
 Resolving package versions...
  Updating `~\.julia\environments\v1.3\Project.toml`
  [4c9194b5] ~ ABCDGraphGenerator v0.1.0 #master (https://github.com/bkamins/ABCDGraphGenerator.jl)
  Updating `~\.julia\environments\v1.3\Manifest.toml`
  [4c9194b5] ~ ABCDGraphGenerator v0.1.0 #master (https://github.com/bkamins/ABCDGraphGenerator.jl)

$ julia deg_sampler.jl degrees.dat 3 5 50 10000 1000 42 100
[ Info: Usage: julia deg_sampler.jl filename τ₁ d_min d_max n max_iter [seed]
[ Info: Example: julia deg_sampler.jl degrees.dat 3 5 50 10000 1000 42
[ Info: Expected value of degree: 8.327743727955891
$ shasum -a 256 degrees.dat #sha256sum degrees.dat on Linux
85abec430e2d03076af9408548dbc09343cd115aa585b57cecc8558cccef1343  degrees.dat

$ julia com_sampler.jl community_sizes.dat 2 50 1000 10000 1000 42 100
[ Info: Usage: julia com_sampler.jl filename τ₂ c_min c_max n max_iter [seed] [nout]
[ Info: Example: julia com_sampler.jl community_sizes.dat 2 50 1000 10000 1000 42 100
[ Info: Expected value of community size: 156.5613820733916
$ shasum -a 256 community_sizes.dat #sha256sum community_sizes.dat on Linux
af57a5a060c0c48884dbec0138c6dbb121311701cd10cbae66d3577323dfb21e  community_sizes.dat

$ julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat xi 0.2 false false 42 100
[ Info: Usage: julia graph_sampler.jl networkfile communityfile degreefile communitysizesfile mu|xi fraction isCL islocal [seed] [nout]
[ Info: Example: julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat xi 0.2 true true 42 100
$ shasum -a 256 network.dat
4fd2e419f5c0bec89c77edaaeb876242b5d2f098e790d3c6febb4a35f3b43d34  network.dat
```
After running these commands you will have the following files in your working directory (all data is 1-based)):
* `degrees.dat` a sequence of vertex degrees (in descending order)
* `community_sizes.dat` a sequence of cluster sizes (in descending order)
* `community.dat` a sequence of vertex number-community number pairs
* `network.dat` a sequence of generated edges sorted lexicographically as pairs of vertices (in increasing order)

If `nout` was passed then the first community consists of outliers.

You can check the generated graph using the `graph_check.jl` file. Here is an example (truncated output):
```
$ julia graph_check.jl degrees.dat community_sizes.dat community.dat network.dat false
[ Info: Example sage: julia graph_check.jl degrees.dat community_sizes.dat community.dat network.dat [isCL]
[ Info: Number of nodes: 10000
[ Info: Number of communities: 69
[ Info: mean required degree: 8.3672
[ Info: min required degree: 5
[ Info: max required degree: 50
[ Info: mean generated degree: 8.3672
[ Info: min generated degree: 5
[ Info: max generated degree: 50
[ Info: mean graph level internal fraction: 0.7955945930519239
[ Info: Internal fractions per community:
[ Info: Community 1 has size 100 and internal fraction 0.052109181141439205
[ Info: Community 2 has size 605 and internal fraction 0.8067756549143195
[ Info: Community 3 has size 564 and internal fraction 0.8110303548525011
[ Info: Community 4 has size 501 and internal fraction 0.80732635585157
[ Info: Community 5 has size 398 and internal fraction 0.8056426332288401
...
[ Info: Community 65 has size 58 and internal fraction 0.7992202729044834
[ Info: Community 66 has size 55 and internal fraction 0.8016528925619835
[ Info: Community 67 has size 54 and internal fraction 0.7933491686460807
[ Info: Community 68 has size 53 and internal fraction 0.7888040712468194
[ Info: Community 69 has size 53 and internal fraction 0.8066037735849056
```
As you can see first community, that had outliers has much lower fraction of
internal edges. All other communities have it around 80% percent as expected.

Note that in extreme cases of very low `xi` the outliers will form a community
(as there are not enough edges from proper communities that allow forming
background graph). In such a case you will get a warning, but the generator
was designed to produce the result. Here is an example:
```
$ julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat xi 0.0 false false 42 100
[ Info: Usage: julia graph_sampler.jl networkfile communityfile degreefile communitysizesfile mu|xi fraction isCL islocal [seed] [nout]
[ Info: Example: julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat xi 0.2 true true 42 100
┌ Warning: Because of low value of ξ the outlier nodes form a community. It is recommended to increase ξ.
└ @ ABCDGraphGenerator graph_sampler.jl:346

$ julia graph_check.jl degrees.dat community_sizes.dat community.dat network.dat false
[ Info: Example sage: julia graph_check.jl degrees.dat community_sizes.dat community.dat network.dat [isCL]
[ Info: Number of nodes: 10000
[ Info: Number of communities: 69
[ Info: mean required degree: 8.3672
[ Info: min required degree: 5
[ Info: max required degree: 50
[ Info: mean generated degree: 8.3714
[ Info: min generated degree: 5
[ Info: max generated degree: 51
┌ Warning: Nodes with not matching degrees are NamedTuple{(:node, :expected, :actual), Tuple{Int64, Int64, Int64}}[(node = 1, expected = 50, actual = 51), (node = 2, expected = 50, actual = 51), (node = 3, expected = 50, actual = 51), ..., (node = 265, expected = 24, actual = 25), (node = 277, expected = 24, actual = 25)]
└ @ Main graph_check.jl:46
[ Info: mean graph level internal fraction: 0.9998998335095137
[ Info: Internal fractions per community:
[ Info: Community 1 has size 100 and internal fraction 0.9950372208436724
[ Info: Community 2 has size 605 and internal fraction 0.9996061441512406
[ Info: Community 3 has size 564 and internal fraction 1.0
[ Info: Community 4 has size 501 and internal fraction 1.0
[ Info: Community 5 has size 398 and internal fraction 1.0
...
[ Info: Community 65 has size 58 and internal fraction 1.0
[ Info: Community 66 has size 55 and internal fraction 1.0
[ Info: Community 67 has size 54 and internal fraction 1.0
[ Info: Community 68 has size 53 and internal fraction 1.0
[ Info: Community 69 has size 53 and internal fraction 1.0
```
As you can see in the truncated output, since we set `xi=0` (asking no
background graph) all communities (including outliers) have almost 100% of
internal edges. The small deviations are due to the fact that in Configuration
Model we need to make sure that each community graph has even sum of degrees.
This also means that we get a warning that we needed to change by at most 1
degrees of some nodes in comparison to their generated degrees to meet this
condition.
