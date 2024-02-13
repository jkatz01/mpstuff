# A2: Parallel Array Merging

```
Partition A into r groups
each with k = logn elements
Group 1: A[1]	... A[k]
Group 2: A[k+1]	... A[2k]
Group i: A[(i-1)k+1]	... A[ik]
Group r: A[(r-1)k+1]	... A[rk]

Find r integers j(1) ... j(r)
j(1) is greatest index so A[k]  >= B[j(1)]
j(2) is greatest index so A[2k] >= B[j(2)]
j(i) is greatest index so A[ik] >= B[j(i)]
j(r) is greatest index so A[rk] >= B[j(r)]

This partitions B into r groups
b[1]	      ...b[j(1)], b[j(1)+1]  ...b[j(2)]
b[j(i-1)+1]...b[j(i)], b[j(r-1)+1]...b[j(r)]

Assign processor i to merge group i of A 
 		     & group i of B

This guarantees that the elments of B have
reached their final position in C(1:2n)
```
