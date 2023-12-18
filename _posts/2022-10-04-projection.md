---
layout: post
title:  📏 Projection in linear subspace
date: 2022-10-04
description: Abstract algebra - projection in terms of inner-products
tags: math notes
categories: Bioinformatics
related_posts: false
featured: true
thumbnail: /assets/img/pca_thumb.png
---

Often times, we want to find an orthogonal projection of a vector onto a subspace. Instead of writing this in matrix multiplication, here is what I like to write, using abstract algebra terms such as inner-product. The inner-product is a great tool to abstractly think about concepts such as angle and distance.

First, a subspace S of $$rank(S)=k$$ can be defined solely based on a chosen set of orthogonal basis $$V=(\vec{w_1},\vec{w_2},...\vec{w_k})$$.

The projection of an arbitrary vector $$x$$ onto $$S$$ is defined as a vector in $$S$$ that is the closest to $$x$$, and it can be written as 

\begin{equation}
proj_S(x)= proj_{w_1} (x) + proj_{w_2} (x) + ... + proj_{w_k} = \sum_{i=1}^k \alpha_i w_i (x)
\end{equation}

Indeed, the projection must be within the span of $$V$$, and  $$\alpha_i w_i = proj_{w_1} (x)$$. Similarly, we write

\begin{equation}
proj_S(x) =  \sum_{i=1}^k \frac{\langle x,v_i \rangle}{\|v_i\|} v_i = \sum_{i=1}^k \langle x,v_i \rangle v_i
\end{equation}
since our basis is chosen to be orthonormal, i.e $$\| v_i \| =1$$ for all $$i=1,2...k$$.
Now, there are different kind of inner-product. In the euclidean space, we are "endowed" with the dot-product inner-product, and so 
\begin{equation}
proj_S(x) = \sum_{i=1}^k x v_i^T v_i  = VV^T x  
\end{equation}
where we define $$V$$ as 
\begin{bmatrix}
V=[w_1 \\ w_2 \\ ... \\ w_n]
\end{bmatrix}

This will be quite useful for more advanced concepts such as SVD (PCA).
