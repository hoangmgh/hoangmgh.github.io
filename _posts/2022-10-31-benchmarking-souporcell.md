---
layout: post
title:  ⚖️ How class imbalance might affect genetic demultiplexing of scRNA-seq data
date: 2023-10-31 21:01:00
description: Benchmarking souporcell for scRNA-seq data in the case of class imbalance
tags: genetic-demultiplexing troubleshooting
categories: sample-posts
thumbnail: assets/img/soup_bad_example3_11500.png
---


 <font size="5">1. Troubleshooting previous known problems </font> 


From my experience with souporcell, I have seen that in some cases, one donor might be completely missing from the assignments. Instead of having 8 assigned donors, for example, you will have 7 donors with an equal number of cells, while one donor have around 20 cells - almost nothing. For example, here is one of our earlier experiment, where each cell has a hashtag identifying its origin (HTO), and we also used souporcell to assign a cell to a  donor:

<div class="row" style="text-align: center">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/soup_misassignment_example.png" title="Soup missasignment example" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Fig.1: An example illustrating missing donors from souporcell experiments 
</div>

In two of our batches, we have 8 samples: the columns are the identities given by the hashtag, and the rows are the assignments from souporcell.  As you can see, $$AO39PO$$ and $$AO38PO$$ are merged as $$Donor_2$$, while $$Donor_7$$ is trivial. The same phenomenon is also evident in the other Batch, DIA_3P_3C. What is similar in both cases is the low number of cells from the messed up Donors, i.e. 200-500 cells. By looking at the formula of the loss function and the fact that gaussian mixture might be very prone to class imbalance, I suspect that this is really the case.

**2. Baseline scenario**

And so, I attempted to simulate this situation by using our in-house data. From our Thyroid Autoimmune disease scRNA-seq dataset, I extracted the top 2000 cells from 7 Donors. All of these are CD45+ cells extracted from the Thyroid tissue. These cells are the one with the highest UMI and $$n_{genes}$$ expressed, so at least we can be ascertain that coverage is not the problem.  As a baseline, I run souporcell demultiplexing with 14,000 cells. In one run, here is the result that I have:

<div class="row" style="text-align: center">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/14k_cells_soup.png" title="Soup missasignment example" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Fig.2: An example showing the stochastic nature of the algorithm
</div>
This is quite interesting! I expect that most of the cells will fall into its desired bin, i.e. the diagonal entries will consist each of approximately 2000 cells. Instead, what I see is the merging of 2 donors into one assignment, and the splitting of one authentic donor into two clusters! ❗❗**This is very alarming!** given that most people will take the souporcell result at face's value. However, the truth is that  you are seeing  just 6 out of 7 donors. 

I then went on to repeat the same run again, but this time with a different parameter. Instead of using $$1024$$ loci as the maximum number of features, I used $$200$$ loci. Hopefully this can fix the issue! Indeed, we observe in the next run that most of the cells got  accurately assigned:

<div class="row" style="text-align: center">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/soup_good_example.png" title="Soup missasignment example" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Fig.3: A highly accurate result from one run with 14k cells, 2000 from each of 7 donors
</div>

I think this is quite alarming, because most labs will likely run souporcell one time  for one experiment and accept the result AS-IS. I also suspect that if you run this multiple times, perhaps you will arrive at a better solution than the one given in Fig.2. However, it might be true that one single run might give a suboptimal solution, local minima. When one tune the parameter drastically, the result changed, even though it might look like we have somewhat an equal number of cells. 

**Imbalance scenario 1:**

Having finished this experiment, I went on to simulate an imbalance situation, where I subset 500 cells randomly from one donor, and keep 2000 cells from each of 7 donors. I then ran souporcell clustering procedure with the following parameters: $$n_{loci} = 200$$, $$n_{loci} = 1024$$, each with several repeats ($$n_{repeats} = 50$$):
<div class="row" style="text-align: center">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/soup_bad_example_12500.png" title="Soup missasignment example" class="img-fluid rounded z-depth-1" %}
    </div>
</div>

<div class="caption">
    An example of imbalance class in which 500 cells are from one Donor and 12000 from the remaining 6 donors. 
</div>

In another simulation of subsetting 500 cells from one donor, here is the result:
<div class="row" style="text-align: center">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/soup_bad_example2_12500.png" title="Soup missasignment example" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Another example of imbalance class in which 500 cells are from one Donor and 12000 from the remaining 6 donors. 
</div>
While we have 500 cells fall nicely into one category, they were merged with cells from another donor :(.
**Imbalance scenario 2:**
<div class="row" style="text-align: center">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/soup_bad_example3_11500.png" title="Soup missasignment example" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Another example of imbalance class in which 500 cells are from one Donor, 1000 from another and 10000 from the remaining 5 donors. 
</div>

**Using known genotypes to improve clustering**

Given that this problem is an optimization problem, where we are trying to find the center of the clusters, which are essentially the genotypes of the donors, we can try to feed  the known genotype data to the algorithm, if such an information is available. Souporcell implements this by using the known genotype as the initial value for the sparse-mixture gaussian model. This makes a lot of sense! Moreover, hypothetically speaking, this should also solve the imbalance class situation as a good initial value has been shown to eleviate this. A good initial value might help us get closer to the true $$\theta_i$$. 

I attempted to create the genotype data for these 7 samples by using the paired-tissue data. While we can generate the genotype from the BAM files that we extract the cells from, this might cause some biases. With these genotypes information, here is the result in the 12500-cell class imbalance case:

<div class="row" style="text-align: center">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="/assets/img/soup_good_example_genotype2.png" title="A good example" class="img-fluid rounded z-depth-1" %}
    </div>
</div>

**What if some cells' identity  are known ?** 

In our lab, we used to generate the corresponding hashtag data of the cells from GEX experiments. The hashtag are not quite reliable to infer the origin of the cells. Given this, I want to incorporate this information to help demultiplex the cells. Obviously, the hashtag identities can't be used as ground-truth - it is in fact quite noisy, but it can get us closer to the ground truth. To simulate this scenario, I randomly select 100 cells from each sample, and assign a highest weight to its cluster of origin while 0 for the remaining clusters. Here is how souporcell's code work. For each mixture centered at $$\theta_i$$, its Tensorflow representation of the probability function is the following:
\begin{equation}
N(x=v_c | \theta_i)  = e^{-\frac{\sqrt{(v_c-\theta_i)^TW(v_c-\theta_i)}}{2} }
\end{equation}

If an entry at loci $$j$$ of $$v_c$$ is not observed, then basically the $$j$$ entry in the diagonal of $$W$$ will be $$0$$, so that unobserved loci will be "muted out". If we consider cell $$v_c$$ to be coming from the cluster with $$\theta_i$$, then we want to set other $$N(x=v_c | \theta_j)=0$$ where $$j \neq i$$. We can do this by setting the diagonal of $$W$$ to be all 1000 or some higher value. That way
$$\sqrt{(v_c-\theta_i)^TW(v_c-\theta_i)}$$ will be large enough. Regardless, I think there might be a better way to convey this in tensorflow. In the future, we might want to encode, for example, some weights to be $$\infty$$, and so the result will be shrink down to $$0$$. Similarly, we can use 2 different weight matrices. One matrix would be the one to convey the missingness of the data, and the other one to convey the true weights of the mixture gaussian.

With this set-up, I tested the imbalance scenarior (500 x1 and 2000 x 6):

<div class="row" style="text-align: center">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/soup_good_example_genotype_100.png" title="Soup good assignment example" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Fig.3: A highly accurate result from one run with 100 cells known from each cluster
</div>

For future experiments, we might want to assign weights from "high-quality" cells obtained from hashtag to assist demultiplexing. 

<font size="5">2. Conclusion </font> 

I think there are a few things I learned from these experiments with souporcell:
- the result given after a single run might be very inaccurate, and is driven by stochasticity. With that being said, rerunning an experiment with a different set of parameter might help
- If one can afford to genotype the samples, it is always better to include the known genotype as part of the run, instead of, for example, matching the found genotype taken from souporcell and match that with a known genotype later on
- I think we should also ask if imbalance class is a weakness of most unsupervised genetic demultiplexing tools out there. This is what I am trying to benchmark on.
