# Minseps

Code to count minimal separating sets and their embeddings in surfaces.

Maximum genus to search is the input to main() in line 66 of minseps.jl

Genus <= 4 should run fine on personal computers, but genus 5 takes several weeks and needs over 600GB of memory in current version.

If you want to run it for genus 5 and have around the recommended memory but run out of memory, replacing the division by n with
division by (n-1) in line 249 of hypermap_search.jl should give a slight reduction in memory requirements.  A modular version 
more suitable to distributed work over a cluster is in progress.
