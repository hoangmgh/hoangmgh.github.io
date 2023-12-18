---
layout: post
title: 📢 📏 Attention-based Transformer, a mathematical primer
date: 2021-11-06
description: Learning the mathematical assumption behind the SoA Deep Learning model  
tags: math genetic-demultiplexing
categories: Bioinformatics
related_posts: false
featured: false
thumbnail: /assets/img/attention.png
---

There are a lot of tutorial on the net that teach the how and why of the dot-product attention Transformer model, such as that from [jalammar.github.io](https://jalammar.github.io/illustrated-transformer/). This is useful if you want to understand transformer by deriving the input, output, intermediate objects through matrix multiplication. While this is absolutely useful, most tutorials out there can be very confusing, and doesn't really paint the big picture of attention mechanism. At the end of the day, I think one has to seperate the attention mechanism from transformer, in which attention is the **engine**. Most tutorials out there are also confusing matrix embedding with the attention mechanism itself!

At the very core, formally speaking, the attention mechanism can be summarized by Bahdanau's equation:

\begin{equation}
Attention(q,D) = \sum_{i=1}^{m} \alpha (q,k_i) v_i
\end{equation}
where we define $$D$$ to be a "database" or dictionary of $$m$$ tuples associating keys to values. Moreover, we denote $$q$$ a query. Here we also have $$\alpha(q,k_i) \in \mathbb{R}$$ to be a kind of mapping, which take in $$q$$ and $$k_i$$ and "spit out" attention weights. The whole process can be viewed as attention pooling. To make it simpler to understand, let's consider the simple case where exactly one of the weights is $$1$$ while the rest are $$0$$. This is very similar to traditional database query, where only one key $$k_i$$ are matched for a given query $$q$$. In this case Attention function will "spit out" the matched value, i.e. $$v_i$$. Thus, we basically map $$q$$ to $$v_i$$. If all weights are equal, this is basically averaging across all database (average pooling). In most case, some of the weights are diminished, while only a few are positive weights and contribute to the output. We can also normalize the weights to sum up to $$1$$ by dividing their sum:
\begin{equation}
\alpha(q,k_i) = \frac{\alpha(q,k_i)}{\sum_{i=1}^m \alpha (q,k_j)}
\end{equation}

In deep learning setting, what we are trying to learn form the data is a set of attention weights, associating a given set of query $$q$$s to all keys $$k_i$$. 

In the case of self-attention, the queries and the values are the same set. In a way, you are trying to infer some kind of correlation/cooccurence between queries and values from the weights, where the co-occurence weight matrix might not be necessarily symmetric. 

Most tutorials out there teaches the concept of self-attention from query, key, value matrices embedding. In the case of self-attention, you are trying to embed the same input matrices of $$n$$ observations into three matrices. If we are start with $$n$$ observations stored in matrix $$X$$ of the size $$\mathbb{R^{n \times f}}$$,  then by multiplying them with $$Q,K,V \in \mathbb{R^{f \times d}}$$,we arrive at 
$$Q = XW_Q$$, $$K=XW_K$$, and $$V=XW_V$$, respectively. What this is doing is abstracting $$n$$ initial observations as a set of $$n$$ queries (in $$Q$$), $$n$$ keys in $$K$$, and $$n$$ values in $$V$$. 

<div class="row">

    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/attention_2.png" title="example image" style="text-align: center" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Starting from one input matrix X, one arrives at 3 different representation Q,K, and V.
</div>

Thus, from the same set of observations, we are representing them differently. However, to tie this back to our Bahdanau's equation, we "bulk" calculate the attention matrix by doing
\begin{equation}
attention(Q,K,V)= softmax(\frac{QK^T}{\sqrt{d_k}})V
\end{equation}
Here you can think of the left component $$softmax(\frac{QK^T}{\sqrt{d_k}})$$ as the matrix of attention weights mentionedin equation (1), where $$V$$ is the matrix of values. The above matrix is actually a kind of transformation, and given a new query embedded as $$Q$$, that s how you map $$Q$$ to the corresponding $$V$$.  Similarly, if given one single query $$q$$, you can arrive at a new value by doing:
\begin{equation}
Attention_{K,V}(q) = softmax(\frac{qK^T}{\sqrt{d_k}})V 
\end{equation}

In conclusion, with the same set of $$n$$ observations, the attention mechanisms allows us to relate one observation to another (and also itself!) via the quasi-similarity matrix $$softmax(\frac{QK^T}{\sqrt{d_k}})$$. Most importantly, this matrix is not necessarily symmetric! Therefore, you can think of the attention as a great tool that helps us learn about the correlation/co-occurence pattern between one input observation and another. In biology, co-occurence of input features is a 