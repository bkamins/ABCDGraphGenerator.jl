# Description of procedure allowing to reproduce Figures 2 and 4 from the companion paper

## Instructions to reproduce Figure 2

### For LFR plot
1. Download and compile code from https://github.com/eXascaleInfolab/LFR-Benchmark_UndirWeightOvp
2. run `benchm -N 250000 -k 25 -maxk 1500 -t1 2.5 -t2 1.5 -minc 50 -maxc 2500 -muw 0.2`
3. run `properties-mixings.ipynb` Jupyter Notebook

### For ABCD plots
1. Install ABCD generator following instructions at https://github.com/bkamins/ABCDGraphGenerator.jl
2. run `julia deg_sampler.jl degrees.dat 2.5 10 2500 250000 1000`
3. run `julia com_sampler.jl community_sizes.dat 1.5 50 2500 250000 1000`
4. run `julia graph_sampler.jl network.dat community.dat degrees.dat community_sizes.dat 0.2 true true`
   (where the last two parameters should be adjusted to match the required version of ABCD)
5. run `properties-mixings.ipynb` in Jupyter Notebook

## Instructions to reproduce Figure 4
1. Follow the installation instructions as for Figure 2
2. Generate LFR or chosen version of ABCD graph with parameters (save the network file to networks.dat):
    * N = 100000
    * t1 = 3
    * mink = 13
    * maxk = 500
    * t2 = 2
    * minc = 50
    * maxc = 2000
    * mu = 0.2
3. Run command `properties_vector.py -a ALG_NAME >> 100k.stats` where `ALG_NAME`
stands for a name under which the algorithm should be presented.
4. Repeat steps 1 and 2 for selected algorithm multiple times (30 was used in the paper)
5. run `properties-example.ipynb` in Jupyter Notebook
