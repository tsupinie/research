#!/bin/csh

# Start a new run
#python run_real_data_case.py --n-ens 40 --base-path /lustre/scratch/tsupinie/05June2009/ --job-name 1kmf-zs25-offtime-05XP --n-cores 960 --mpi-config 2 12 --ens-start 10800 --ens-end 14400 --ens-step 300 --assim-step 300 --initial-conditions /lustre/scratch/tsupinie/05June2009/3km-ensemble-tar/3kmf-control/ --subset-ic --covariance-inflation 0:mult=1.03 5400:mult=1.20,adapt=0.90 --piecewise --split-files --fcst-req 0:20 --init-fcst-req 0:50 --debug --restart

# Restart a run
#python run_real_data_case.py --n-ens 40 --base-path /lustre/scratch/tsupinie/05June2009/ --job-name 1kmf-ebr-no-mm-05XP-snd --n-cores 960 --mpi-config 2 12 --ens-start 10800 --ens-end 14400 --ens-step 300 --assim-step 300 --covariance-inflation 0:mult=1.03 5400:mult=1.20,adapt=0.90 --piecewise --split-files --restart --debug

# Do the ensemble forecast
#python run_real_data_case.py --n-ens 40 --base-path /lustre/scratch/tsupinie/05June2009/ --job-name 1kmf-zs25-offtime-05XP --n-cores 960 --mpi-config 2 12 --ens-start 14400 --ens-end 18000 --ens-step 300 --assim-step 300 --piecewise --chunk-size 1800 --split-files --subset-ic --boundary-conditions /lustre/scratch/tsupinie/05June2009/3km-ensemble-tar/3kmf-control/ --restart --free-forecast --debug

# Fix shit when kraken fucks up ...
python run_real_data_case.py --members 1 40 --base-path /lustre/scratch/tsupinie/05June2009/ --job-name 1kmf-zs25-offtime-05XP --n-cores 48 --mpi-config 2 12 --ens-start 14700 --ens-end 15000 --ens-step 300 --assim-step 300 --split-files --split-init auto --piecewise --chunk-size 300 --free-forecast --free-fcst-req 0:20 --restart --debug

#### 3km run, r0h=12km
#python run_real_data_case.py --n-ens 40 --base-path /lustre/scratch/tsupinie/05June2009/ --job-name 3kmf-mult=1.03 --n-cores 960 --mpi-config 2 12 --ens-start 0 --ens-end 10800 --ens-step 1800 --assim-step 1800 --covariance-inflation 0:mult=1.03 5400:mult=1.03,adapt=0.90 --piecewise --fcst-req 0:45 --init-fcst-req 1:30 --chunk-size 1800 --split-files --debug

#python run_real_data_case.py --n-ens 40 --base-path /lustre/scratch/tsupinie/05June2009/ --job-name 3kmf-mult=1.03 --n-cores 960 --mpi-config 2 12 --ens-start 10800 --ens-end 14400 --ens-step 300 --assim-step 300 --covariance-inflation 0:mult=1.03 5400:mult=1.03,adapt=0.90 --piecewise --fcst-req 0:20 --init-fcst-req 0:40 --chunk-size 300 --split-files --debug --restart

#### 3km forecasts
#python run_real_data_case.py --n-ens 40 --base-path /lustre/scratch/tsupinie/05June2009/ --job-name 3kmf-mult=1.03 --n-cores 960 --mpi-config 2 12 --ens-start 14400 --ens-end 18000 --ens-step 300 --assim-step 300 --piecewise --chunk-size 1800 --split-files --restart --free-forecast --debug

#### 3km run, n0r=8e5
#python run_real_data_case.py --n-ens 40 --base-path /lustre/scratch/tsupinie/05June2009_alt/ --job-name 3kmf-n0r=2e6 --n-cores 960 --mpi-config 2 12 --ens-start 0 --ens-end 10800 --ens-step 1800 --assim-step 1800 --covariance-inflation 0:mult=1.20,adapt=0.90 --piecewise --fcst-req 0:45 --init-fcst-req 1:30 --chunk-size 1800 --split-files --debug --restart

#python run_real_data_case.py --n-ens 40 --base-path /lustre/scratch/tsupinie/05June2009_alt/ --job-name 3kmf-n0r=2e6 --n-cores 960 --mpi-config 2 12 --ens-start 10800 --ens-end 14400 --ens-step 300 --assim-step 300 --covariance-inflation 0:mult=1.20,adapt=0.90 --piecewise --fcst-req 0:20 --init-fcst-req 0:40 --chunk-size 300 --split-files --debug --restart

#python run_real_data_case.py --n-ens 40 --base-path /lustre/scratch/tsupinie/05June2009_alt/ --job-name 3kmf-n0r=2e6 --n-cores 960 --mpi-config 2 12 --ens-start 14400 --ens-end 18000 --ens-step 300 --assim-step 300 --piecewise --chunk-size 1800 --split-files --restart --free-forecast --debug
