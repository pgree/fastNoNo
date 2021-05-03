        implicit real*8 (a-h,o-z)
        character(1000) arg

        call prini(6,13)

        ix = 1
        call get_command_argument(ix,arg, nlength)
ccc        call prina('args*', args, 100)

        call dense_eval_w_err(arg, nlength)

        stop
        end
c
c
c
c
c
        subroutine dense_eval_w_err(dir, n)
        implicit real *8 (a-h,o-z)
        real*8 a(10 000 000), x(100 000), y(1000 000), dexps(10 000),
     1     dsums2(10 000), dsums(100 000), dds(10 000), stds(10 000),
     1     stds2(10 000)
        character(1000) filename
        character(n) dir

c
c        parameters
c
        filename = dir // '/params.dat'
        call read_params(filename, n, k1, k2, a, y)
        k = k1+k2
        call prinf('n*', n, 1)
        call prinf('k1*', k1, 1)
        call prinf('k2*', k2, 1)

c
c        integrate
c
        nn = 60
        nn_theta = 10
        call dense_eval(nn_theta, nn, n, k1, k2, a, y, dsums, dsum, 
     1     stds)

c
c        double nodes in all directions
c
        nn2 = 2*nn
        nn_theta2 = 2*nn_theta
        call dense_eval(nn_theta2, nn2, n, k1, k2, a, y, dsums2, dsum2, 
     1     stds2)

c
c        check error
c
        call dd_abs_max(dsums2, dsums, k+3, dd_max1)
        call prin2('max posterior mean error*', dd_max1, 1)

        call dd_abs_max(stds2, stds, k+3, dd_max)
        call prin2('max posterior std error*', dd_max, 1)

c
c        write results to file
c
        filename = dir // '/exps.dat'
        call write_dd_exps_stds(filename, k, dd_max1, dsums, stds)

        
        return
        end

