# Running the benchmarks

The benchmarks are run by calling one of the benchmark functions: 
* `benchmark_n_p`: Benchmark on random graphs with `n` vertices with probability `p` of being connected, with `n` ranging from `n_min` to `n_max`.
* `benchmark_n_fixed_m`: Benchmark on random graphs with `n` vertices and `m` edges,  with `n` ranging from `n_min` to `n_max`.
* `benchmark_m_fixed_n`: Benchmark on random graphs with `n` vertices and `m` edges,  with `m` ranging from `m_min` to `m_max`.
* `benchmark_geo`: Benchmark on random **geometric** graphs with `n` vertices with probability `p` of being connected, with `n` ranging from `n_min` to `n_max`.
* `benchmark_dimacs`: Benchmark on the graphs from the DIMACS dataset.

Examples of running all of these functions is given in the file `main.cpp`.

The DIMACS benchmark assumes the following files are present in the `dimacs/` directory:

* USA-road-d.NY.gr
* USA-road-d.BAY.gr
* USA-road-d.COL.gr
* USA-road-d.FLA.gr
* USA-road-d.NW.gr
* USA-road-d.NE.gr
* USA-road-d.NY.co
* USA-road-d.BAY.co
* USA-road-d.COL.co
* USA-road-d.FLA.co
* USA-road-d.NW.co
* USA-road-d.NE.co

These can be downloaded from http://www.diag.uniroma1.it//challenge9/download.shtml or https://web.archive.org/web/20220208020227/http://www.diag.uniroma1.it//challenge9/download.shtml.