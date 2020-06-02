
# Code for "Financial Fragility with SAM?"
## Daniel Greenwald, Tim Landvoigt, and Stijn Van Nieuwerburgh<sup>1</sup>


Fast Steps to Reproduce Benchmark Model
=======================================

-   In Matlab, execute the script "`main_create_env`".

    -   If you have the Matlab Symbolic Math Toolbox and use version
        2018b or earlier, leave the flag "`useJacobian`" in line 43 on.
        Otherwise set to false to use numerical differentiation.

    -   Note: generating the analytic Jacobian for the benchmark model
        takes approximately 5 minutes with version 2018b, and can take
        longer for the other experiments.

    -   "`main_create_env`" will create a file "`env_bench_ini0`" that
        contains the experiment definition for the benchmark economy.

    -   The code archive contain the pre-computed file.

-   Execute script "`main_run_exper`".

    -   You can set the number of parallel workers to be started in the
        separate script "`open_parpool`"

    -   Set to zero if you want to run it with a single process.

    -   On a computer with sixteen cores (and 16 parallel workers) the
        benchmark model converges in about 3 hours.

    -   "`main_run_exper`" creates a results file named
        "`res_[current_date_time]`" that contains the converged policy
        functions.

    -   The code archive contain the pre-computed named
        "`res_20191112_bench`".

-   Simulate the model using "`sim_stationary`" and
    "`sim_trans_cluster`".

    -   "`sim_stationary`" simulates the model policies contained in
        "`res_20191112_bench`" for 10,000 periods and writes out the
        resulting time-series and several statistics. The main output is
        a file named "`sim_res_20191112_bench`".

    -   "`sim_trans_cluster`" reads both "`res_20191112_bench`" and
        "`sim_res_20191112_bench`", and simulates generalized IRFs.

    -   To plot IRFs, run "`plot_trans`".

For More Details See readme.pdf
=================



<sup>1</sup>: Greenwald: Massachussetts Institute of Technology Sloan School;
    email: dlg(at)mit.edu. Landvoigt: University of Pennsylvania Wharton
    School, NBER, and CEPR; email: timland(at)wharton.upenn.edu. Van
    Nieuwerburgh: Columbia University Graduate School of Business, NBER,
    and CEPR, 3022 Broadway, Uris Hall 809, New York, NY 10027; email:
    svnieuwe(at)gsb.columbia.edu.
