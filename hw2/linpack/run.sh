# download and unzip
wget https://software.intel.com/content/dam/develop/external/us/en/documents/l_onemklbench_p_2021.2.0_109.tgz
tar -xvzf l_onemklbench_p_2021.2.0_109.tgz

# run benchmark
benchmarks_2021.2.0/linux/mkl/benchmarks/linpack/xlinpack_xeon64 -i config

# clean up
rm -rf benchmarks_2021.2.0 l_onemklbench_p_2021.2.0_109.tgz