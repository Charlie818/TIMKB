## Build && Usage

__Input__

- edges file

The topic-aware edges file, ending with `.tedge` extension, starts with two ints: number of edges, number of topics at the first line. The following lines each contains two ints each denotes one user id, and `# of topic` floats each represents one topic dimension.

Example .tedge file:

```
390160 2
1 2 0.2 0.8
...
132423 49939 0.2 0.8
...
```

_Note_, I use space, not tab.

You can intepret the above file as:
This file has 390160+1 lines, starting from the second line, each line represents one edge in the graph, in total 390160 edges. This graph has 2 topics. The first edge (the second line) is from user 1 to user 2 (directional), and the topic vector is `<0.2, 0.8>` which means this relationship emphasis more on topic 2.

- query file

The query file, ending with `.query.items` extension, contains queries you can feed to the program at once. It starts with two ints at the first line, number of queries and number of topics respectively. The following lines each contain `# of topics` floats each represents one topic dimension.

Example .query.items file

```
20 2
0.2 0.8
0.3 0.7
```

You can intepret the above file as:
This file has 20 queries, each query has `2-` dimension topic. The second line represents query NO.1, which has the query topic vector `<0.2, 0.8>`. `k` is specified in command line, not in this file.


__Commands__

Given the input file, run these commands sequentially:

```
# preprocessing
$ make clean && make offline
$ ./offline -l ./data/diggs.tedge 1000 10 ./data/diggs.list

# 1000 and 10 are the two theta thresholds which controls the influence spread. 1000 is converted to 1/1000 in the code (0.001, only influence strength less than 0.001 is ignored), we use 1000 only to make it easier to type in shell.

# querying (queries are read from a file)
$ make clean && make experiments
# best effort
$ ./experiments -b ./data/diggs.tedge ./data/diggs.list 1000 10 ./data/diggs.query.items 50
```

__Output__

The printed output contains other information that is only useful in debugging/experiments. The useful result in each line is the `6th and 7th` floats, which corresponds to time cost and influence spread respectively.

```
956, 153, 0.056244, 0.094, 51, 0.152592, 214.741, 214.191, 1000, 10, 50
962, 158, 0.051039, 0.088306, 51, 0.141398, 299.757, 299.248, 1000, 10, 50
1030, 158, 0.053409, 0.090449, 51, 0.145898, 375.602, 374.731, 1000, 10, 50
999, 177, 0.055424, 0.073973, 51, 0.131366, 154.649, 154.384, 1000, 10, 50
...
```

You can intepret this output as:
Query No.1 takes `0.152592` seconds, with influence spread `214.741 `, and `50` seeds are selected.

**Note**

As you run this program, you may not get the result as quicly as the value printed out, for example, it seems to take more than `0.152592` seconds for the example query NO.1. That's because it takes time to calcuate the influence spread using monte carlo simulations. Please be aware of this.

## Code architecture

### brief introduction to the files.

* `PMIA/` folder contains the `PMIA` baseline source code.
* `src/` folder contains source code for this paper.
	* `experiments.cpp` is the `main` function. Take a look at `experiments.*.sh` to see how to use it, `experiments.sh` defines only which dataset to use.
	* `global.h` defines some macros that are set mainly for debugging.
	* `offline.cpp` and `offline.sh` runs the offline index building process.
	* `pmia_no_topic_spread.cpp` and `pmia_no_topic_spread.cpp.sh` builds graph for the `PMIA` baseline. (Invoking the binary built in `PMIA/` folder)
	* `tests.cpp` is only used during dubugging.
	* `tic.cpp` is where graph definition resides.
	* `tim.cpp` contains all the algorithms proposed in the paper.
	* `kmeans.py` is the topic sample selection method in the paper.
* `data/` contains datasets used in the experiment.
* `figs/` will contain the experiment result files after you run experiments. I know it should really be named `exps/` or `res/`, I was just thinking of experiment figures, so, that name.

### A walk through the code

To give you a feeling how things are grouped together and to generate experiment results shown in paper, I will walk you through this process right from the beginning (use shell scripts to invoke binaries).

First of all, look at those `experiments.*.sh` shell scripts. Those are written to streamline the experiment processes. You can tell from the name which part of experiment it is running. For example, `approx` means Approximation (Topic-Sample Approach in the paper, section 5). `experiments.bound.sh` compares different bounds, in figure 4. `itex` is actually `MIS` baseline, `greedy` and `pmia` both use code in `PMIA/` code, not our code. `inflex` uses our implementation in `src/`. `greedy` results were not shown in the paper (maybe because it's too slow if I remember it correctly). All `experiment.*.sh` scripts follow a similar pattern, define which dataset to use, which queries to run, set parameters to test, define the output filename (in `figs/`) which contains the meaning of the experiment, for example `figs/diggs.approx.sample` means it is the result of Topic sample approach on diggs dataset with different sampling ratios (figure 5 in the paper). All experiments will run queries in a batch, ie, run all queries defined in a query file. All experiments will invoke the binary `experiments` (compiled and linked from `experiments.cpp`), except that `pmia` and `greedy` only use `experiments` to convert a topic-unaware graph into a topic-aware graph under the query, and uses `PMIA/` code to actually run the algorithms. `MIS` will invoke `TopicAwareCode/`.

Next we talk about dataset preparation. Above section `Build & Usage` gives the input data format, but there are some scripts in `data/` need clarification. `dblp_generate_itemset.py` generates "action logs" because dblp don't have such things. `dblp_generate_queries` generates queries. `itemset_2_queries.py` is more general, it converts itemset into queries file. (itemset is another name of topic-sample).

Now we discuss offline preparation process. There are three offline index building processes. `./offline -l` means influence historgram (the predecessor of local graph method, but not mentioned in the paper), it generates a `$(DATA).list` file (for example, `diggs.list`), each row is a histogram of a vertex's influence with different, descending theta value (lower theta, less cut-off, more influence). `-uset` generates index for in-edge (neighborhood-based in paper) method, named `$(DATA).uset`, each row is the maximum influence of a vertex under each topic. `-us` and `-ls` generetes upper sample and lower samples (hints) in Topic sample approach, it needs `*.up.sample` and `*.low.sample` files to tell it which topic distributions to find seeds on. `*.up/low.sample` files are the same to `itemset` files, each row is a topic distribution vector.

Now we can move on to `experiments.cpp` file. It supports several command line arguments. `-p` is used for PMIA preprocessing, to assign topic-aware propagation probability on edges under a query (actual PMIA algorithm is in `PMIA/`). `-l` is used for preprocessing process of `MIS`. __But it is the only part I can't remember correctly right now, some files are missing, it is commented out in `experiments.itex.sh`__. `-n` are actually not used. `-b`, `-f`, `-u` and `-r` are all best-effort algorithms with different bounds. `-f` is originally designed to represent a best-effort algorithm without bound estimation (direct compute), but its results was not shown in paper. `-b` stands for best-effort algorithm with local graph bound estimation method, `-u` stands for neighborhood based estimation method, and `-r` stands for precomputation based bound estimation. `-a` is the Topic sample approach (approximation in the code). `-i` is Inflex we implemented.

We now move on to some important functions/variables in the code.

in `global.h`, macro `VERBOSE`/`EXAMPLE`/`EXP` are used to control the output for debugging (switch on/off in different scenario). This dumb, actually. macro `NBARS` defines the influence historgram granularity, `NBARS`=100 means there will be 100 bars.

`pmia_no_topic_spread` is originally written so see how well PMIA alone (no topic-aware) can do in our topic-aware setting. It turned out poorly, and the results were not shown in the paper.

`tests.cpp` can be ignored.

`tic.cpp` is where topic-aware graph is built. The most tricky part is that graph is represented as list of edges, and sorted by `u` (edge source, since edge is directed). Also, since the lower bound estimation of Topic-sample approach needs a quick computation of margainal ab(activation probability) of each influencees, there are outedges list and inedges list. It contains a `Dijkstra` method to get the initial influence.

`tim.cpp` should be easy to understand given the paper. Main confusion comes from misnomers. `Local Graph` in paper is `Bound/Bounded` in the code. `Neighborhood-base` is `InEdge` in code. `LB` always stands for lower bound, same goes for `UB` upper bound. One important implementation mistake is the `InEdge` function, which after several different versions under many discussions is not the same as in the paper. Also, I remember `Bound()` method uses influence histogram, as I recall that used to be the first-to-go bound estimation method. Another important detail is that, in this code, ap (activation probability) calculation `MarginalAPOf(const int v, const int p, const double w)` calculates the marginal influence of `p` over `v`, with the propagation probability `w`. `p` is always `v`'s in-neighbor (one hop). There is a problem here, `p`'s own ap (`ap_[p]` in the code) is its global ap, not just inside `v`'s MIIA (in the word of PMIA paper). So I have to say this is not right. But that's what got into our paper. This one is a mistake, I remembered I argued about this at that time.

Other things I think might be helpful to a new-comer (who hasn't seen Chen Wei's PMIA code before) is the meaning of those class variables.

* `ap_` stands for current activation probability of a vertex, * `tap_` stores a temporary ap for calculating the marginal influence. 
* `H_` is the heap in best-effort framework. `HNode` is a node in the heap. I forget why there is a `LNode`, as far as I can remember it is not used.
* `Q` is the query data structure.
* `theta_` is the cut-off threshold, and `theta_2` is used in local graph method (defines the local area).
* `round_` is current top-`i` iteration of all `k` iterations.
* `dist_`, `seen_`, `seen_idx_`, `children_`, `pred_` are for dijkstra. `pred_` stands for predecessor. `dist_` is distance.
* `is_border_` means the border of a local graph.
* `hist_` historgram.
* `L_` is the offline topic-unaware initial influence of a vertex. In best-effort framework, vertex `u` is added into the heap (ordered by) according to their initial influence (`L_[u]`). `cursorL_` keeps the last added vertex.
* `max_topic_` is used in neighborhood-based estimation.
* `used_` means selected as seed.
* `lap_` is similar to `tap_`, but used in lower bound estiamtion for topic sample approach.
* `ubsamples_` and `lbsamples_` are for topic sample approach also.

`Inflex` is implemented by a `WCopeland` rank function. The other advanced algorithm described in that paper is not included, only one that claimed to be the best is used.

_Details can be found in the comments_


If there is still anyting I didn't cover, please let me know.