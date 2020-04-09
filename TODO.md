# To do

* The base cases are most important: most of their distances will be short (in many cases). Hence, we should avoid
  selecting "representatives" (taxa names like "#17# in my implementation) in favor of the "local candidates".
  There are cases when almost all vertices are representatives, but that is not the norm, and should be avoided.
    - Implement this!
* Better distance estimation. Poisson model sucks
    - ML on aligned sequences done.
    - ML on unaligned sequences is ongoing (student).
* Model chosen by input type. Right now assuming *protein* alignments, but could be DNA alignments too.
* When two sequence do not share columns or they are completely different, I
  default to distance 2.5. This is not good in all models (like Poisson!).
    - Is this still the case??
