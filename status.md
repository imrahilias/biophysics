# todo

- why does the illumination order change?
- parallelise clustering
- compile clustering
- find more efficient 3d clustering

# problems

- variable illumination order has to be adjusted manually unless fixed
- rigid projection demands sets of the same size so some regions are skipped
- rigid projection does not cover sheering/scaling/etc
- more efficient 3d clustering
- affine projection is probably needed anyway for sheering/scaling

# done

- find efficient rigid 3d projection optimisation routine
- find efficient affine 3d projection optimisation routine
- simulate rigid transform > rekover rigidly
- simulate rigid transform > rekover affinely
- simulate affine transform > rekover affinely
- simulate affine transform > rekover rigidly (fails)
