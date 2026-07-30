#!/bin/zsh
# Per-model surface-type breakdown: the acceptance corpus for stage-2 guards.
# A quick 14-model look. tools/sweep.sh is the full harness (every model,
# --compare distances, baselines and diffs) and is what the measurement
# protocol requires; this stays because it is fast.
for m in \
  onshape/stepped_block_toroidal_fillets_fine onshape/stepped_block_toroidal_fillets_medium \
  onshape/rounded_cube_10_r2_fine onshape/part_rounded_cube_10_r2_coarse \
  onshape/pipe_elbow_10_medium onshape/nosecone_medium \
  onshape/pill_fine onshape/chamfered_cube_10_c1_medium \
  ccad/generated/cone_cylinder ccad/generated/filleted_cylinder ccad/generated/filleted_pipe \
  nurbs/bspline_wave nurbs/elliptical_arch nurbs/elliptical_cylinder ; do
  b=$(./target/release/dump_regions tests/$m.stl surfaces 2>/dev/null \
      | awk '{for(i=1;i<=NF;i++) if($i=="plane")p++; else if($i=="cyl")c++; else if($i=="cone")k++;
              else if($i=="sphere")s++; else if($i=="torus")t++; else if($i=="bspline")b++}
             END{printf "%d/%d/%d/%d/%d/%d", p+0,c+0,k+0,s+0,t+0,b+0}')
  printf "%-17s %s\n" "$(basename $m | cut -c1-17)" "$b"
done
