---
layout: post
title:  📈 Krockner product
date: 2023-02-17 10:25:00
description: "Formality" for Tensor Decomposition
tags: math
categories: Bioinformatics
related_posts: false
featured: true
thumbnail: 
---

To understand the tensor decomposition, the first thing I tried to learn is the Kronecker product, which can be seen as a generalization of the vector scaling. Instead of scaling the vector with a value, you are scaling a matrix with multiple values:


For two matrices $$A$$ and $$B$$ living in $$\mathbb{R^{m \times n}}$$, we have the matrix product to be defined as 
$$A B= A [u_1 u_2 ... u_n ] = [ Au_1 Au_2 ... Au_n] $$
 
Sometimes you can think of matrix multiplication as a set of operations which are all  linear combinations ($$Au_i$$) (which I think is much easier to memorize)

For Kronecker product, instead of multiplying each vector with $$A$$, you are actually getting a multitude of $$B$$ using all  entries in $A$. This means you are making $$p x q$$ different copies of the same matrix.

