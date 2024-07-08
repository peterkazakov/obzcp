# OBZCP
Algorithm for finding Odd Binary Z-Complimentary Pairs

# Description
For all questions regarding this algorithm, please refer to 
Peter Kazakov, Zilong Liu, "A Computer Search of Primitive OBZCPs of Lengths up to 49"
(link will be published soon)
or contact directly the authors.

# Usage
Prerequisites: Install (latest) Rust on your machine

Go to `src/main.rs` and modify values of `N` and `MAX_ACC`, in case search for non-optimal OBZCP is needed.

Run 
`cargo run --release` and wait for the output.

# Notes

Minimum value of `N` is 5.
Since the complexity of this algorithm is `O(2^N)`, executions above `N>40` are slow and above `N>50` might be impossible.

