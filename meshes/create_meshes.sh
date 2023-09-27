

# ran on chicoma and copied files over
# source ~/repos/compass/master/load_dev_compass_*_gnu_mpich.sh

for N in 16 32 64 128 256
do
  let DC=64000*16/$N
  echo 'creating mesh with N='$N', DC='$DC
  planar_hex --nx ${N} --ny ${N} --dc $DC -o base_mesh_${N}x${N}.nc
  MpasMeshConverter.x base_mesh_${N}x${N}.nc mpas_mesh_${N}x${N}.nc
done
