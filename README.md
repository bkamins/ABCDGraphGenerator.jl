# ABCDGraphGenerator.jl
Artificial Benchmark for Community Detection (ABCD) - A Fast Random Graph Model with Community Structure

Bogumił Kamiński, Paweł Prałat, François Théberge

---

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

The resason for such split of the functionality is that generation of vertex degrees
and community sizes is fast, while the generation of the final graph is the most expensive step.

The `utils/` folder contains CLI utilities that are aimed at users that want
to the use package without using the Julia language. A requirement for these
utilities to be run is to have the Julia language in version at least 1.0 installed on a computer.
It contains the following files:
* `install.jl`: installs all required packages
* `deg_sampler.jl`: samples degrees of vertices in the graph
* `com_sampler.jl`: samples communitiy sizes in the graph
* `graph_sampler.jl`: samples edges and community assignments in the graph

Here is an output from an example session using CLI:
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

$ julia deg_sampler.jl degrees.dat 3 5 50 10000 1000
[ Info: Usage: julia deg_sampler.jl filename τ₁ d_min d_max n max_iter
[ Info: Example: julia deg_sampler.jl degrees.dat 3 5 50 10000 1000
[ Info: Expected value of degree: 8.327743727955891

$ julia com_sampler.jl community_sizes.dat 2 50 1000 10000 1000
[ Info: Usage: julia com_sampler.jl filename τ₂ c_min c_max n max_iter
[ Info: Example: julia com_sampler.jl community_sizes.dat 2 50 1000 10000 1000
[ Info: Expected value of community size: 156.5613820733916

$ julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat 0.2 true true
[ Info: Usage: julia graph_sampler.jl networkfile communityfile degreefile communitysizesfile μ isCL islocal
[ Info: Example: julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat 0.2 true true
```
After running these commands you will have the following files in your working directory (all data is 1-based)):
* `degrees.dat` a sequence of vertex degrees (in descending order)
* `community_sizes.dat` a sequence of cluster sizes (in descending order)
* `community.dat` a sequence of vertex number-community number pairs
* `network.dat` a sequence of generated edges sorted lexicographically as pairs of vertices (in increasing order)
