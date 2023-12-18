---
layout: post
title: 📢 Genetic-demultiplexing with Souporcell -  GPU implementation
date: 2023-11-15
description: speed up demultiplexing with GPU Tensorflow 
tags: math genetic-demultiplexing
categories: Bioinformatics
related_posts: false
thumbnail: /assets/img/tf.jpg
---
I rewrote souporcell so it can be run on GPU with half of the time required. This can be particularly useful if you need to tune parameters several times for an experiment

Then the whole [souporcell.py](https://github.com/wheaton5/souporcell/blob/master/souporcell_pipeline.py) can be broken down do three components:
read_mtx(), cluster_step(). The cluster step can benefit from running on GPU. Basically we can replace ***tf.*** with ***tf.compat.v1***, and run the Gradient Descent step with the with     ***with tf.device("/GPU:0")*** :
The first step is to read in the matrices (ref.mtx and alt.mtx), filter out loci that are expressed in at least $$min_{alt}$$ and $$min_{ref}$$ cells, and finally create a list of possible loci and their corresponding indices to be used. We call this the "pool of indices", from where we will go on to the cluster step and do more filtering:

```python
def read_mtx(min_alt = 5,min_ref = 5,K = 7,max_loci =1024,alt_matrix="alt_simulated.mtx",ref_matrix="ref_simulated.mtx"):
    cell_index = {}
    total_lost = 0
    loci_counts = {}
    cell_counts = {}
    with open(alt_matrix) as alt:
        alt.readline()
        alt.readline()
        tokens = alt.readline().strip().split()
        cells = int(tokens[1])
        cell_counts = {x: {} for x in range(1, cells + 1)}
        total_loci = int(tokens[0])
        for line in alt:
            tokens = line.strip().split()
            locus = int(tokens[0])
            cell = int(tokens[1])
            cell_counts.setdefault(cell, {})
            count = int(tokens[2])
            cell_counts[cell][locus] = [0, count]
            loci_counts.setdefault(locus, [0, 0])
            if count > 0:
                loci_counts[locus][1] += 1
    with open(ref_matrix) as alt:
        alt.readline()
        alt.readline()
        alt.readline()
        for line in alt:
            tokens = line.strip().split()
            locus = int(tokens[0])
            cell = int(tokens[1])
            count = int(tokens[2])
            cell_counts[cell][locus][0] = count
            loci_counts.setdefault(locus, [0, 0])
            if count > 0:
                loci_counts[locus][0] += 1

    used_loci_set = set()
    used_loci = []
    for (locus, counts) in loci_counts.items():
        if counts[0] >= min_ref and counts[1] >= min_alt:
            used_loci.append(locus - 1)
            used_loci_set.add(locus - 1)
    used_loci = sorted(used_loci)
    used_loci_indices = {locus:i for (i, locus) in enumerate(used_loci)}
    loci = len(used_loci)
    return used_loci_indices,used_loci_set,used_loci,loci,loci_counts,cell_counts
```

The next step is called the cluster step, where we build the fractional matrices for all cell. However, each cell will have their own set of loci and might not be overlapping with other sets from other cells. This is quite interesting! It is also important that one can use up to $$max_{loci}$$ loci (default =1024). Regardless, the set of loci for each cell must satisfy:
 - they must be pooled from the pool of cells identified earlier (used_loci_indices)
 - the loci are considered if they are expressed in the given cell ( ref_count > 0 and alt_count > 0)

I am adding a parameter for known_cells, where the indices of known cells and their corresponding cluster number must be provided (a dataframe)
This works by assigning a maximal weight to the known cluster:
```python
weightshape_np=np.array(np.transpose(np.broadcast_to(weights.T,(K,max_loci,weights.shape[0]))))

   
        
    if known_cells  is not None:
        assert len(known_cells.keys()) <= K
        assert known_cells.columns==["index","cluster"]
        #did a quick check for unique assignment to a cluster:
        freq=known_cell["cluster"].value_counts()
        assert np.max(freq)==1
        
        for i,j in enumerate(known_cells.keys()):
            for k in range(K):
                if k!=i:
                    weightshape_np[known_cells[j]-1,:,k]=weightshape_np[known_cells[j]-1,:,k]*10000    
```
and so the likelihood of seeing a known cell's genotype vector is:
\begin{equation}
 P(x=v_c) = e^{-\frac{1}{2}\|v_c-\theta_k\|}
\end{equation}
given that cell $$c$$ is drawn from  a Gaussian centered at $$\theta_k$$
Instead of assigning the weights with 10000, perhaps there is a better way to represent "maximal distance" of $$v_c$$ from other $$\theta_1$$, $$\theta_2$$,... but not $$\theta_k$$...
```python



def cluster_step(max_loci,K,training_epochs,repeats,cell_counts,loci_counts,used_loci_indices,known_cells=False,min_ref=5,min_alt=5,lr=.1):
    print("loci being us based on min_alt, min_ref, and max_loci "+str(loci))
    cells = len(cell_counts)
    total_lost = 0

    cell_data = np.zeros((cells, max_loci))
    cell_loci = np.zeros((cells, max_loci))
    
    weights = np.zeros((cells, max_loci))
    
    for cell in cell_counts.keys():
        index = 0
        single_cell_counts = cell_counts[cell]
        
        #prioritize locus that is highly expressed across cells:
        #for this given cell, get the set of 
        
        this_cell_locus=list(cell_counts[cell].keys())
       
             
        for locus in this_cell_locus:
            locus_counts = single_cell_counts[locus]
            if loci_counts[locus][0] >= min_ref and loci_counts[locus][1] >= min_alt:
                if index < max_loci:
                    ref_c = locus_counts[0]
                    alt_c = locus_counts[1]
                    if ref_c + alt_c > 0:
                        cell_data[cell - 1][index] = float(ref_c)/float(ref_c + alt_c) if ref_c + alt_c > 0 else 0.0
                        cell_loci[cell - 1][index] = used_loci_indices[locus - 1]
                        weights[cell - 1][index] = 1.0
                        index += 1
                    total_lost += 1
    ### set 0 weights for cells from other clusters:
    random_size=100


    data = cell_data
    data_loci = cell_loci
    weightshape_np=np.array(np.transpose(np.broadcast_to(weights.T,(K,max_loci,weights.shape[0]))))

   
        
    if known_cells  is not None:
        assert len(known_cells.keys()) <= K
        assert known_cells.columns==["index","cluster"]
        #did a quick check for unique assignment to a cluster:
        freq=known_cell["cluster"].value_counts()
        assert np.max(freq)==1

        for i,j in enumerate(known_cells.keys()):
            for k in range(K):
                if k!=i:
                    weightshape_np[known_cells[j]-1,:,k]=weightshape_np[known_cells[j]-1,:,k]*10000    
    
    cells = data.shape[0]

    #save this for investigation:
    print("save cell_data and cell_loci:")
     
    rng = np.random
    
    
    import tensorflow as tf
    
    session_conf = tf.compat.v1.ConfigProto(
          allow_soft_placement=True)
    session_conf.gpu_options.allow_growth = True
    session_conf.gpu_options.per_process_gpu_memory_fraction =1
    tf.compat.v1.disable_eager_execution()
    
    tf.compat.v1.reset_default_graph()
    
    with tf.device("/GPU:0"):
    

        #init = tf.compat.v1.constant(sample_genotypes.T)
        #phi = tf.compat.v1.get_variable(name="phi", initializer = init, dtype = tf.float64)
        
        phi = tf.compat.v1.get_variable(name = "phi",shape=(loci, K), initializer = tf.initializers.random_uniform(minval = 0, maxval = 1), dtype = tf.float64)

        input_data = tf.compat.v1.placeholder("float64", (cells, max_loci)) #tf.constant("input",np.asmatrix(data))
        input_loci = tf.compat.v1.placeholder("int32", (cells, max_loci))
        loci_per_cell = tf.compat.v1.placeholder("float64", (cells))
        trans = tf.compat.v1.transpose(input_data)
        broad_trans = tf.compat.v1.broadcast_to(trans,[K,max_loci,cells])
        untrans = tf.compat.v1.transpose(broad_trans)
        xtest = untrans-tf.compat.v1.gather(phi,input_loci)
        weight_data = tf.compat.v1.placeholder("float64", (cells,max_loci,K)) #tf.constant("weights",np.asmatrix(weights))
    
        weighted = weight_data*xtest
        powtest = -tf.compat.v1.pow(weighted,2)
        post = tf.compat.v1.reduce_sum(powtest,axis=1)
        logsum = tf.compat.v1.reduce_logsumexp(post,axis=1)
        cost = -tf.compat.v1.reduce_sum(logsum)
    
        optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate=lr).minimize(cost)
            
        posteriors = []
        min_cost = None
        weightsshape=[]
        logsum_list=[]
        cluster_list=[]
    for repeat in range(repeats):
        init = tf.compat.v1.global_variables_initializer()
        print("repeat "+str(repeat))
        training_epochs = 1000
        last_cost = None
        with tf.compat.v1.Session(config = session_conf) as sess:
            sess.run(init)
            for epoch in range(training_epochs):
                sess.run(optimizer, feed_dict={input_data:data, weight_data:weightshape_np, input_loci:data_loci})
    
                if epoch % 10 == 0:
                    c = sess.run(cost, feed_dict={input_data:data, weight_data:weightshape_np, input_loci:data_loci})
                    print("epoch "+str(epoch)+" "+str(c))
                    #if last_cost and ((last_cost-c)/c) < 0.0001:
                    if min_cost and last_cost and c > min_cost and (last_cost - c)/(c - min_cost) < 0.005:
                        print("bailing out, too little progress toward minimum so far")
                        break
                    if last_cost and last_cost - c < 1:
                        last_cost = None
                        break
                    last_cost = c
            if min_cost:
                min_cost = min(min_cost, c)
            else:
                min_cost = c
    
            posterior = sess.run(post, feed_dict={input_data:data, weight_data:weightshape_np, input_loci:data_loci})
            posteriors.append((c,posterior))
                
                
    
    sess.close()
    
    posterior = sorted(posteriors, key=lambda x: x[0])
    posterior = posterior[0][1]    
    clusters = np.argmax(posterior,axis=1)
    return(clusters)

```
And then from your GPU-equipped computer, run the following:

```python
used_loci_indices,used_loci_set,used_loci,loci,loci_counts,cell_counts=read_mtx(alt_matrix="alt.mtx",
                                                                                ref_matrix="ref.mtx")
clusters=cluster_step(max_loci=100,K=7,training_epochs=1000,repeats=30,cell_counts=cell_counts,
                      loci_counts=loci_counts,used_loci_indices=used_loci_indices,
                      cluster_tmp="cluster_simulated.tsv",known_cells=False,min_ref=5,min_alt=5,lr=.1)
```

I randomly select 100 cells from each known cluster from an earlier experiment, and here is the performance in the case of imbalance-dataset:

