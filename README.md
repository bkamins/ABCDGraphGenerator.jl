# ABCDGraphGenerator.jl
Artificial Benchmark for Community Detection (ABCD) - A Fast Random Graph Model with Community Structure

Bogumił Kamiński, Paweł Prałat, François Théberge

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
