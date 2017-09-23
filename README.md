# DSSA-Importance-Sampling-SKIS

Information:
--------------------------------------------------------
Implementation of D-SSA Algorithms on SKIS Importance Sketch for Influence Maximization under Independent Cascade(IC) or Linear Threshold(LT) model. For more details about D-SSA, please read our paper:

[H. T. Nguyen, M. T. Thai, and T. N. Dinh, Stop-and-Stare: Optimal Sampling Algorithms for Viral Marketing in Billion-scale Networks , in Proc. of the annual ACM SIGMOD/PODS Conference, 2016] (https://arxiv.org/abs/1605.07990)

and of SKIS, please refer to other paper:

[H. T. Nguyen, T. P. Nguyen, NhatHai Phan, T. N. Dinh, Importance Sketching of Influence Dynamics in Billion-scale Networks, IEEE International Conference on Data Mining series (ICDM), 2017](https://arxiv.org/abs/1709.03565)

Contact Authors: Hung T Nguyen (hungnt@vcu.edu)
		 Thang N Dinh (tndinh@vcu.edu)


Terms of use:
--------------------------------------------------------
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.


Requirements:
--------------------------------------------------------
In order to compile all the tools, it requires GCC 4.7.2 and later.


Compile:
--------------------------------------------------------

Use `make' command to compile everything


How to use:
--------------------------------------------------------
This package offers a set of functions to use in order to find a seed set of given size k:
0. (Optional) Computing edge weights (probabilities) as described in the experiments:

	./format <input file> <output file> 1

	<input file>: the path to text file in edge list format with no weights: the first line contains the number of nodes n and number of edges m, each of the next m lines describes an edge following the format: <src> <dest>. Node index starts from 1.
	<output file>: the path to text output file with edge probabilities
	The last parameter (1) means the input graph is considered as directed.

1. Conversion from a text format to binary file

	./el2bin <input file> <output file>

    	<input file>: the path to text file in weighted edge list format: the first line contains the number of nodes n and number of edges m, each of the next m lines describes an edge following the format: <src> <dest> <weight>. Node index starts from 1.
    	<output file>: the path to binary output file

2. Run SSA to find the seed sets

        ./dssa_skis [Options]

    Options:

        -i <binary graph file>
            specify the path to the binary graph file (default: network.bin)

        -o <seed output file>
            specify the path to the output file containing selected seeds (default: network.seeds)

        -k <number of seeds>
            number of selected seed nodes in the seed set(default: 1)

        -epsilon <epsilon value used>
            epsilon value in (epsilon,delta)-approximation (see our paper for more details, default: 0.1)

        -delta <delta value used>
            delta value in (epsilon,delta)-approximation (default: 1/n)

        -m <model>
            diffusion model (LT or IC, default: LT)

     Output format:

        The outputs are printed on standard output stream in the following order

		Total samplings: <total samples generated>
                Seed Nodes: <list of selected seed nodes>

                Influence: <Spread of Influence of the select seed set>
                Time: <running time in seconds>
                Memory: <peak memory used>

3. (Optional) Verify influence spread of a seed set - returns a (epsilon, 1/n)-estimate of the influence:

        ./verifyInf <binary graph file> <seed file> <epsilon> <number of threads> <model: LT or IC>

********************************************************************************************************

Examples on DBLP network to find seed set of 100 seed nodes:

	1. Convert to binary file:

                ./el2bin dblp_format.txt dblp.bin

        2. Run DSSA+SKIS with k=100, epsilon=0.1,delta=0.01:

                ./dssa_skis -i dblp.bin -k 100 -epsilon 0.1 -delta 0.01 -m IC

        The output:

           	Total samplings: 40960
		Seed Nodes: 1300 4246 144 560 2692 404 315 200 7277 1256 1362 934 2290 26 3428 445 435 1844 2130 261 4251 141 11049 8223 4 4681 3340 241 2588 7217 7001 270 3498 7347 606 2452 322 959 54 921 925 13105 5356 4248 3069 3866 935 7975 6241 3895 211 6207 9191 1602 3489 5092 1968 16 391 7871 11380 7985 14360 2021 3293 1563 97 4279 1162 6330 2100 910 2647 4281 6761 5445 3460 11077 4562 2749 12206 21186 8946 6431 6990 4242 5207 3057 1671 7487 2500 2963 9465 3017 3540 5665 11421 6982 8741 7124 

		Influence: 73468.25
		Time: 0.27s
		Memory: 119.684
