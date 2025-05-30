# Minseps

Code to count minimal separating sets and their embeddings in surfaces.

Maximum genus to search is the input to main() in line 66 of minseps.jl

Genus <= 4 should run fine on personal computers, but genus 5 takes just under 2 weeks and needs somewhere about 600GB of memory in current version (recent changes have reduced the exact amount needed somewhat and this will be 
updated once we've computed the new memory bound).
Partial results can be obtained more quickly for genus 5 by running minseps_modular, which only considers a single edge count per run in genus 5.
This edge count should be specified in line 153 of minseps_modular.jl
Memory requirements for minseps_modular.jl are much lower for E< 14, but the same as minseps.jl for larger edge counts.

If you want to run it for genus 5 and have around the recommended memory but run out of memory, replacing the division by n with
division by (n-1) in line 249 of hypermap_search.jl should give a slight reduction in memory requirements.  A version which never stores any
minimal separating sets, only isomorphism classes of underlying graphs, and which will be more suitable for distributed computing is in progress.
