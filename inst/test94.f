c
c
c
c
c
        subroutine test_integration()
        implicit real *8 (a-h,o-z)
        real*8 a(20 000 000), y(500 000)
        real *8, allocatable :: dsums(:), dexps(:),
     1     dsums2(:), dds(:), stds(:), stds2(:)
        character(100) csv_file, filename

c
c        parameters
c
        csv_file = 'test_dense_params.dat'
        csv_file = 'abort_params.dat'
        csv_file = 'params.dat'
        call read_params(csv_file, n, k1, k2, a, y)
ccc        call prin2('a = *',a,n*k)
ccc        call prin2('y = *',y,n)
        k = k1+k2
        call prinf('n*', n, 1)
        call prinf('k1*', k1, 1)
        call prinf('k2*', k2, 1)

        k3 = k+3
        allocate(dsums(k3))
        allocate(dsums2(k3))
        allocate(dexps(k3))
        allocate(dds(k3))
        allocate(stds(k3))
        allocate(stds2(k3))

c
c        integrate
c
        nn = 60
        nn_theta = 20
        nn_theta = 4
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
        call dd_abs_max(dsums2, dsums, k+3, dd_max)
        call prin2('max posterior mean error*', dd_max, 1)

        call dd_abs_max(stds2, stds, k+3, dd_max)
        call prin2('max posterior std error*', dd_max, 1)

        filename = 'exps.dat'
        call write_exps_stds(filename, k, dsums, stds)
ccc        call prin2('dsums*', dsums,k+3)

c
c        stan comparison
c
c        csv_file = 'means.dat'
c        csv_file = 'test_means_dense_stan2.dat'
c        csv_file = 'test_means_dense_stan.dat'
c        call read_means(csv_file, k, dexps)
cccc        call prin2('stan expectations*', dexps, k+3)
c
c        call dd_abs_max(dexps, dsums, k, dd_max)
c        call prin2('stan max error*', dd_max, 1)
        
        return
        end
c
c
c
c
c
        subroutine utest2()
        implicit real *8 (a-h,o-z)
        real*8 a(10 000 000), x(100 000), y(1000 000), dexps(10 000),
     1     dsums2(10 000), dsums(100 000), dds(10 000), stds(10 000),
     1     stds2(10 000)
        character(100) csv_file, filename

c
c        compare the results of the algorithm of this file with 
c        the hardcoded results from a high accuracy run from this
c        same algorithm
c

c
c        read in unit test parameters
c
        csv_file = 'test_dense_params2.dat'
        call read_params(csv_file, n, k1, k2, a, y)
        k = k1+k2
        call prinf('n*', n, 1)
        call prinf('k1*', k1, 1)
        call prinf('k2*', k2, 1)

c
c        integrate
c
        nn = 80
        nn_theta = 60
        call dense_eval(nn_theta, nn, n, k1, k2, a, y, dsums, dsum, 
     1     stds)

c
c        compare the result computed here with hardcoded results
c
        filename = 'test_means_dense_fast2.dat'
        call read_exps_stds(filename, k, dsums2, stds2)

        call dd_abs_max(dsums2, dsums, k+3, dd_max)
        call prin2('max posterior mean error*', dd_max, 1)

        call dd_abs_max(stds2, stds, k+3, dd_max)
        call prin2('max posterior std error*', dd_max, 1)

c
c        stan comparison
c
        filename = 'test_means_dense_stan2.dat'
        call read_means(filename, k, dexps)
ccc        call prin2('stan expectations*', dexps, k+3)

        call dd_abs_max(dexps, dsums, k, dd_max)
        call prin2('stan max error*', dd_max, 1)
        
        return
        end
c
c
c
c
c
        subroutine utest()
        implicit real *8 (a-h,o-z)
        real*8 a(10 000 000), x(100 000), y(1000 000), dexps(10 000),
     1     dsums2(10 000), dsums(100 000), dds(10 000), stds(10 000),
     1     stds2(10 000)
        character(100) csv_file, filename

c
c        compare the results of the algorithm of this file with 
c        the hardcoded results from a high accuracy run from this
c        same algorithm
c

c
c        read in parameters from unit test
c
        csv_file = 'test_dense_params.dat'
        call read_params(csv_file, n, k1, k2, a, y)
        k = k1+k2
        call prinf('n*', n, 1)
        call prinf('k1*', k1, 1)
        call prinf('k2*', k2, 1)

c
c        integrate
c
        nn = 80
        nn_theta = 60
        call dense_eval(nn_theta, nn, n, k1, k2, a, y, dsums, dsum, 
     1     stds)

c
c        compare to hardcoded results
c
        filename = 'test_means_dense_fast.dat'
        call read_exps_stds(filename, k, dsums2, stds2)

        call dd_abs_max(dsums2, dsums, k+3, dd_max)
        call prin2('max posterior mean error*', dd_max, 1)

        call dd_abs_max(stds2, stds, k+3, dd_max)
        call prin2('max posterior std error*', dd_max, 1)

c
c        and compare to stan
c
        filename = 'test_means_dense_stan.dat'
        call read_means(filename, k, dexps)
ccc        call prin2('stan expectations*', dexps, k+3)

        call dd_abs_max(dexps, dsums, k, dd_max)
        call prin2('stan max error*', dd_max, 1)
        
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging code and the beginning of the
c       least squares subroutines proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c         This file contains 8 user-callable subroutines: nrleastsq,
c         nrleastsq2, nrleamatlr, rleamatr, nrleamatrr, nrleamatll,
c         nrleamatrl, nrleas_proj. Following is a brief description 
c         of these subroutines.
c
c   dense_eval - evaluates two-group normal normal model with 
c         normal+(0, 1) hyperpriors on the scale parameters 
c 
c
c   eval_rho_int - compute the inner integral after a change of 
c         variables of the density of this problem 
c
c
c
c
c
c
c
c
        subroutine dense_eval(nn_theta, nn, n, k1, k2, a, y, dsums, 
     1     dsum, stds)
        implicit real *8 (a-h,o-z)
        real*8 a(*), y(*), dsums(*), stds(*),
     1     thetas(nn_theta),rhos(nn),phis(nn), s(k1+k2),s2(k1+k2), 
     1     xs(k1+k2),
     1     ysmall(k1+k2),ys2(k1+k2),vmoms(100 000), dsumsi(k1+k2+3+10),
     1     whts_thetas(nn_theta+10),whts_phis(nn+10),whts_rhos(nn+10),
     1     dsum_xsi(k1+k2+10), dsum_vars(k1+k2+10)
        real *8, allocatable :: u(:,:),vt(:,:),asca(:,:),ut(:,:),
     1     b(:,:),ynew(:),v(:,:),bb(:,:), fmaxs(:), fs(:,:),
     1     dsums_cov(:,:), dsums_covi(:,:), wwti(:,:), xxti(:,:),
     1     gs(:,:), c(:,:), ct(:,:), tmp_mat(:,:)

c
c        this subroutine computes posterior expectations of the two
c        group bayesian linear regression model
c
c          y ~ normal(A1*x1 + A2*x2, sig1)
c          x1 ~ normal(0, sig2)
c          x2 ~ normal(0, sig3)
c          sig{1,2,3} ~ normal+(0, 1)
c
c        this subroutine computes E[x_i], var(x_i), E[sig{1,2,3}]. 
c        after a change of variables, from sig1, sig2, sig3 to 
c        spherical coordinates, we integrate over a 3-dimensional 
c        grid using a tensor product of gaussian nodes. 
c        for each value of the outer integral
c        we do a diagonalization where the conditional density is 
c        a gaussian with a diagonal covariance.
c

        done = 1.0d0
        pi = 4*atan(done)

        k = k1+k2

        allocate (u(n,k))
        allocate (ut(k,n))
        allocate (vt(k,k))
        allocate (v(k,k))
        allocate (asca(n,k))
        allocate (b(k,k))
        allocate (bb(k,k))
        allocate (ynew(k))
        allocate (fmaxs(nn))
        allocate (fs(nn_theta,nn))
        allocate (gs(nn_theta,nn))
        allocate (dsums_cov(k,k))
        allocate (dsums_covi(k,k))
        allocate (wwti(k,k))
        allocate (xxti(k,k))
        allocate (c(k,k))
        allocate (ct(k,k))
        allocate (tmp_mat(k,k))

c
c        before computing integral find least squares solution to
c        ax = y as well as projection of y onto the
c        columns of a
c
        call dsvd(n,k, a, u, s, vt)
        call mat_trans(u, n, k, ut)
        call mat_mult(ut, k, n, a, k, b)
        call y_proj(y,u,s,n,k,ysmall,resid)
        call mat_vec_mult(ut, k, n, y, ynew)

c
c        construct grid for outer parameter of integration, theta
c
        theta0 = 0.0
        theta1 = pi/2.0d0
        call theta_lege_nodes_whts(nn_theta, theta0, theta1, thetas, 
     1     whts_thetas)

c
c        initialize
c
        ktmp = k+3
        call set_zero(ktmp, dsums)
        dsum = 0
        call set_zero(k*k, dsums_cov)
        fm = -1.0d250

c
c        theta integral
c       
        do 130 i=1,nn_theta
        t = thetas(i)
        wht_theta = whts_thetas(i)
c
c        compute the diagonal form of a^t*a for each theta which
c        will be used to convert the integrand to a diagonal 
c        gaussian for each (theta, phi, rho) 
c
        call rescale_a(t,b,asca,k,k,k1,k2)
        call dsvd(k,k, asca, u, s, vt)
        call y_proj(ynew,u,s,k,k,ysmall,res2)
        call entry_sqr(ysmall,ys2,k)
        call entry_sqr(s,s2,k)

c
c        initialize sums
c
        dsumi = 0
        ktmp = k+3
        call set_zero(ktmp, dsumsi)
        call set_zero(k**2, dsums_covi)
        call set_zero(k**2, wwti)
        call set_zero(k, dsum_vars)
        fmi = -1.0d250

c
c        for each theta, get upper integration bound for phi
c
        call get_phi1_fmax(nn, n, k, t, ysmall, ys2, s, s2, resid,
     1     phi1i, fmax_theta)
        fmaxs(i) = fmax_theta
ccc        call prin2('phi1i*', phi1i, 1)
ccc        call prin2('fmax*', fmax, 1)

c
c        ... and lay down nodes 
c
        phi0 = 0.0
        call lege_nodes_whts(nn, phi0, phi1i, phis, whts_phis)
ccc        call prin2('phis new*', phis, nn)

c
c          phi integral 
c
        do 120 j=1,nn
        phi = phis(j)
        call compute_beta_alpha_prefact(phi,ysmall,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)
        call compute_mjs(k,s2,ysmall,phi,vmoms)

        wht_phi = whts_phis(j)

c
c          rho integral
c
        call eval_rho_int(nn, n, k, exp_fact, 
     1     dsum_rho, dsum_rho1, dsum_rho2, fmax)
c
c        for fixed theta, compute integral over phi by a sum of 
c        the form \sum_i exp(fi)*gi so that
c        we have an expression exp(fmi)*dsum, the reason we do 
c        this is that exp(fi) is often smaller than 10^{-250}
c       
        fi = prefact + fmax - fmaxs(1)
        gi = dsum_rho
        gi1 = dsum_rho1
        gi2 = dsum_rho2
        wt = wht_phi
        fs(i, j) = fi
        gs(i, j) = gi
        if (fi .gt. fmi) then
          dsumi = dsumi * exp(fmi-fi) + gi*wt

          do ijk=1,k
          dsumsi(ijk) = dsumsi(ijk) * exp(fmi - fi) + gi*wt*vmoms(ijk) 
          enddo        

          dsumsi(k+1) = dsumsi(k+1)*exp(fmi-fi) + gi1*cos(phi)*wt
          dsumsi(k+2) = dsumsi(k+2)*exp(fmi-fi) + gi1*sin(phi)*cos(t)*wt
          dsumsi(k+3) = dsumsi(k+3)*exp(fmi-fi) + gi1*sin(phi)*sin(t)*wt

          do ijk=1,k
          coef = 1/(s2(ijk)/cos(phi)**2 + 1/sin(phi)**2)
          dsum_vars(ijk) = dsum_vars(ijk) * exp(fmi-fi) + coef*gi2*wt
          enddo        

          do i1=1,k
          do i2=1,k
          wwti(i1, i2)=wwti(i1,i2)*exp(fmi-fi)+gi*vmoms(i1)*vmoms(i2)*wt
          enddo
          enddo

         fmi = fi
        else
          dsumi = dsumi + exp(fi-fmi) * gi*wt

          do ijk=1,k
          dsumsi(ijk) = dsumsi(ijk) + exp(fi - fmi) * gi*wt*vmoms(ijk) 
          enddo        

          dsumsi(k+1) = dsumsi(k+1) + exp(fi-fmi)*gi1*cos(phi)*wt
          dsumsi(k+2) = dsumsi(k+2) + exp(fi-fmi)*gi1*sin(phi)*cos(t)*wt
          dsumsi(k+3) = dsumsi(k+3) + exp(fi-fmi)*gi1*sin(phi)*sin(t)*wt

          do ijk=1,k
          coef = 1/(s2(ijk)/cos(phi)**2 + 1/sin(phi)**2)
          dsum_vars(ijk) = dsum_vars(ijk) + exp(fi-fmi)*coef*gi2*wt
          enddo        

          do i1=1,k
          do i2=1,k
          wwti(i1, i2)=wwti(i1,i2)+exp(fi-fmi)*gi*vmoms(i1)*vmoms(i2)*wt
          enddo
          enddo

        endif

 120    continue

c        
c        for covariance we need to convert back to original coordinate
c        system. we have a different change of variables for each value
c        of theta, so do this for each theta
c
        call mat_trans(vt, k, k, v)
        call get_xs_to_ws_matrix(k,k1,k2,v,t,c,ct)
        call mat_mult_dense_diag(c, k, k, dsum_vars, tmp_mat)
        call mat_mult(tmp_mat, k, k, ct, k, dsums_covi)

ccc        call prin2('vars*', dsum_vars, k)
ccc        call prin2('vars*', dsums_covi, 1)

c
c        recover matrix E[x*x^t]
c
        call mat_mult(c, k, k, wwti, k, tmp_mat)
        call mat_mult(tmp_mat, k, k, ct, k, xxti)

        call get_xs_from_ws(dsumsi,k,k1,k2,vt,t, dsum_xsi)
ccc        call prin2('xs 2*', dsums_xsi, k)
ccc        call prin2('dsumi*', dsumi, 1)

c
c        due to underflow issues, integrate over theta by computing 
c        sum of the form \sum_i exp(fi)*gi such that at the end
c        we have an expression exp(fm)*dsum 
c
        fi = fmi
        wt = wht_theta
        if (fi .gt. fm) then
          dsum = dsum * exp(fm-fi) + dsumi*wt

          do ijk=1,k
          dsums(ijk) = dsums(ijk) * exp(fm-fi) + dsum_xsi(ijk)*wt
          enddo

          dsums(k+1) = dsums(k+1) * exp(fm-fi) + dsumsi(k+1)*wt
          dsums(k+2) = dsums(k+2) * exp(fm-fi) + dsumsi(k+2)*wt
          dsums(k+3) = dsums(k+3) * exp(fm-fi) + dsumsi(k+3)*wt

          do i1=1,k
          do i2=1,k
          tmp = dsums_covi(i1,i2)+xxti(i1,i2)
          dsums_cov(i1,i2)=dsums_cov(i1,i2)*exp(fm-fi) + tmp*wt
          enddo
          enddo
  
          fm = fi
        else
          dsum = dsum + exp(fi-fm) * dsumi*wt

          do ijk=1,k
          dsums(ijk) = dsums(ijk) + exp(fi-fm)*dsum_xsi(ijk)*wt
          enddo

          dsums(k+1) = dsums(k+1) + exp(fi-fm)*dsumsi(k+1)*wt
          dsums(k+2) = dsums(k+2) + exp(fi-fm)*dsumsi(k+2)*wt
          dsums(k+3) = dsums(k+3) + exp(fi-fm)*dsumsi(k+3)*wt

          do i1=1,k
          do i2=1,k
          tmp = dsums_covi(i1,i2)+xxti(i1,i2)
          dsums_cov(i1,i2)=dsums_cov(i1,i2) + exp(fi-fm)*tmp*wt
          enddo
          enddo
  
        endif

 130    continue



c
c        esimate error by looking at the coefficients of our
c        function in a 2d legendre expansion
ccc        call get_err(nn_theta, nn, fs, gs)

c
c        scale first and second moments by normalizing constant
c
        do i=1,k+3
        dsums(i)=dsums(i)/dsum
        enddo

        do i=1,k
        do j=1,k
        dsums_cov(i, j)=dsums_cov(i, j)/dsum
        enddo
        enddo

c
c        dsums_cov now contains E[xx^t], adjust to get covariance 
c
        do i=1,k
ccc        tmp_mat(i) = dsums_cov(i, i)
        do j=1,k
        dsums_cov(i, j)=dsums_cov(i, j) - dsums(i)*dsums(j)
        enddo
        enddo

ccc        call prin2('dsums_cov*', dsums_cov, k**2)
ccc        call prin2('dsums*', dsums, k)
ccc        call prin2('tmp_mat*', tmp_mat, k)


c
c        get stds
c
        do i=1,k
        stds(i) = sqrt(dsums_cov(i, i))
        enddo

c
c        plot
c        
ccc        call prin2('fs*', fs, nn_theta*nn)
ccc        call plot_heatmap(nn_theta, nn, fs)


        return
        end
c
c
c
c
c
        subroutine get_err(nnt, nn, fs, gs)
        implicit real *8 (a-h,o-z)
        real*8 fs(nnt, nn), gs(nnt, nn), hs(nnt, nn), u(100 000),
     1     v(100 000), x(10 000), u_theta(100 000), v_theta(100 000),
     1     coefs(nnt, nn), coefs_log(nnt, nn), whts(100 000),
     1     ut(100 000), tmp_mat(100 000)

c
c        this subroutine is not finished! it is in the middle of 
c        development
c

        call max_vec(nnt*nn, fs, fmax, ind)
        
        do i=1,nnt
        do j=1,nn
        hs(i, j) = fs(i, j) - fmax + log10(gs(i, j))
        hs(i, j) = exp(fs(i, j) - fmax)*gs(i, j)
        enddo
        enddo
ccc        call plot_heatmap(nnt, nn, hs)

c
c        legendre expansion
c
        itype=2
        call legeexps(itype,nnt,x,u_theta,v_theta,whts)
        call legeexps(itype,nn,x,u,v,whts)

        call mat_trans(u, nn, nn, ut)

        call mat_mult(u_theta, nnt, nnt, hs, nn, tmp_mat)
        call mat_mult(tmp_mat, nnt, nn, ut, nn, coefs)

ccc        call prin2('coefs*', coefs, nnt*nn)


        dsum = 0.0d0
        ndsum = 0
        tot_coef = 0.0d0
        do i=1,nnt
        do j=1,nn
        coefs_log(i, j) = log10(abs(coefs(i,j)))
        tot_coef = tot_coef + abs(coefs(i,j))
        
        if ((i .gt. nnt-3) .or. (j .gt. nn-5)) then
          dsum = dsum + abs(coefs(i,j))
          ndsum = ndsum + 1
        endif

        enddo

        if (i .eq. nnt) then
          call prin2('coefs(i, j)*', coefs(i, j), 1)
        endif

        enddo

ccc        call prin2('coefs_log*', coefs_log, nnt*nn)
        call prin2('coefs*', coefs, nnt*nn)
        call plot_heatmap(nnt, nn, coefs_log)


c
c        error estimate
c
ccc        call prinf('ndsum*', ndsum, 1)
        avg_tail_coef = dsum / (nn*3 + nnt*5 - 3*5)
        call prin2('avg_tail_coef*', avg_tail_coef, 1)
        err = avg_tail_coef / tot_coef 

        call prin2('err*', err, 1)

        fmax = -1.0d250
        do i=1,nn
        if (fmax .lt. coefs(nnt, i)) fmax = coefs(nnt, i)
        enddo

        do i=1,nnt
        if (fmax .lt. coefs(i, nn)) fmax = coefs(i, nn)
        enddo

ccc        call prin2('fmax*', fmax, 1)
ccc        call prin2('fmax2*', coefs(nnt, nn), 1)

        return
        end
c
c
c
c
c
        subroutine get_xs_to_ws_matrix(k,k1,k2,v,t,a,at)
        implicit real *8 (a-h,o-z)
        real*8 s1(k),v(k,k),a(*),at(*)

c
c        construct matrix that goes from diagonal coordinate
c        system (depending on t) to original coordinate system
c

        do i=1,k1
        s1(i) = cos(t)
        enddo
        
        do i=1,k2
        s1(i+k1) = sin(t)
        enddo

        call mat_mult_diag_dense(s1, k, v, a)
        call mat_trans(a, k, k, at)

        return
        end
c
c
c
c
c
        subroutine get_phi1_fmax(nn, n, k, t, ysmall, ys2, s, s2, resid,
     1     phi1i, fmax)
        implicit real *8 (a-h,o-z)
        real*8 rhos(nn+10),whts_rhos(nn+10),phis(nn+10),
     1     whts_phis(nn+10), ysmall(*),ys2(*),s(*),s2(*),
     1     fmaxs(nn+10)

c
c        find the upper integration bound of the integral with 
c        respect to phi. do this by evaluating the integral on  
c        a sparse grid and checking when the value decreases below
c        10^{-18} of its maximum. also return the maximum value of 
c        the integrand 
c

        done = 1.0d0
        pi = 4*atan(done)

c
c        lay down nodes
c
        phi0 = 0.0d0
        phi1 = pi/2.0d0
        call lege_nodes_whts(nn, phi0, phi1, phis, whts_phis)

c
c        tabulate function
c
        fmax = 0.0d0
        do 119 j=1,nn
        phi = phis(j)
        call compute_beta_alpha_prefact(phi,ysmall,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)

c
c        get maximum of integrand and use that maximum for the 
c        first value of phi for scaling
c
        dn = n - 2
        rho_max = sqrt(1.0d0/2.0d0*(sqrt(4.0d0*exp_fact + dn*dn) - dn))
        call eval_logdens_rho(rho_max,n,exp_fact,fmaxmax)
        fmaxs(j) = fmaxmax

c
c        compute integral in rho
c
        call eval_rho_int(nn, n, k, exp_fact, 
     1     dsum_rho, dsum_rho1, dsum_rho2, fmaxi)
        fi = log(dsum_rho) + prefact + fmaxi - fmaxs(1)

        if (fi .gt. fmax) fmax = fi

        thresh = log(1.0d-22)
        if (fi - fmax .lt. thresh) then
          phi1i = phi
          goto 121
        endif
 119    continue

c
c        must use full interval 
c
        phi1i = pi/2.0d0

 121    continue

c
c        we're going to use the log maximum of the integrand later,
c        so return that value
c
        fmax = fmax + fmaxs(1)

        return
        end
c
c
c
c
c
        subroutine eval_rho_int(nn, n, k, exp_fact, 
     1     dsum_rho, dsum_rho1, dsum_rho2, fmax)
        implicit real *8 (a-h,o-z)
        real*8 whts_rhos(nn),rhos(nn)

c
c        compute the integral from 0 to infinity of 
c         f(rho), rho*f(rho), and rho**2*f(rho)
c
c        where 
c
c         f(rho) = exp(-rho**2/2 - exp_fact/(2*rho**2))/rho**(n-2)
c

c
c        first find integration bounds
c       
        call get_rho0(exp_fact, n, rho0)
        call get_rho1(exp_fact, n, rho1)
        call lege_nodes_whts(nn, rho0, rho1, rhos, whts_rhos)

c
c        get maximum of integrand 
c
        dn = n - 2
        rho_max = sqrt(1.0d0/2.0d0*(sqrt(4.0d0*exp_fact + dn*dn) - dn))
        call eval_logdens_rho(rho_max,n,exp_fact,fmax)
        
c
c        compute integral
c
        dsum_rho = 0
        dsum_rho1 = 0
        dsum_rho2 = 0
        do 110 i=1,nn
        rho = rhos(i)
        call eval_logdens_rho(rho,n,exp_fact,f)
        wht_rho = whts_rhos(i)

        dsum_rho = dsum_rho + exp(f-fmax)*wht_rho
        dsum_rho1 = dsum_rho1 + rho*exp(f-fmax)*wht_rho
        dsum_rho2 = dsum_rho2 + rho**2*exp(f-fmax)*wht_rho

 110    continue


        return
        end
c
c
c
c
c
        subroutine compute_mjs(k,s2,ys,phi,vmoms)
        implicit real *8 (a-h,o-z)
        dimension vmoms(*),s2(*),ys(*)

c
c        compute conditional expecations of gaussian
c

        sp = sin(phi)
        sp2= sp*sp
        cp = cos(phi)**2

        do i=1,k
        vmoms(i) = ys(i)*sqrt(s2(i))*sp2/
     1          (s2(i)*sp2+cp)
        enddo

        return
        end
c
c
c
c
c
        subroutine eval_logdens_rho(rho,n,exp_fact,val)
        implicit real *8 (a-h,o-z)

c
c        evaluate the log of the integrand of the inner integral 
c

        val = -exp_fact/2.0d0/(rho*rho)
        val = val -rho*rho/2.0d0
        val = val - (n-2)*log(rho)

        return
        end
c
c
c
c
c
        subroutine compute_beta_alpha_prefact(phi,ys,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)
        implicit real *8 (a-h,o-z)
        dimension ys(*),ys2(*),s(*),s2(*)

c
c        for convenience, compute quantities that appear in 
c        several locations in the integrand. prefact is the
c        part of the integrand that doesn't depend on rho
c

        done = 1.0d0
        pi   = atan(done)*4

        sp = sin(phi)**2
        cp = cos(phi)**2

c
c        alpha
c
        alpha = 0

        do i=1,k
        alpha = alpha - log(cp+sp*s2(i))
        enddo

        alpha = alpha/2 + k/2.0d0*log(2*pi)

c
c        beta
c
        beta = 0

        do i=1,k
        beta = beta + ys2(i)/(s2(i)*sp+cp)
        enddo

c
c        prefact
c
        prefact = -(n-k)/2.0d0*log(cp)+alpha + log(sin(phi))
        exp_fact= resid/cp+beta

        return
        end
c
c
c
c
c
        subroutine get_xs_from_ws(ws,k,k1,k2,vt,t, xs)
        implicit real *8 (a-h,o-z)
        dimension xs(*),ws(*),vt(k,k)

c
c        convert expectations back from the diagonal coordinate (ws)
c        system to the original coordinate system (xs)
c
        do i=1,k
        z = 0
        do j=1,k
        z = z + vt(j,i)*ws(j)
        enddo
        xs(i) = z
        enddo

        do i=1,k1
        xs(i) = xs(i)*cos(t)
        enddo
        
        do i=1,k2
        xs(i+k1) = xs(i+k1)*sin(t)
        enddo

        return
        end
c
c
c
c
c
        subroutine entry_sqr(xin,xout,k)
        implicit real *8 (a-h,o-z)
        dimension xin(*),xout(*)
c
c        square the entries of an array
c
        do i=1,k
        xout(i) = xin(i)**2
        enddo

        return
        end
c
c
c
c
c
        subroutine y_proj(y,us,ss,n,k,ysmall,resid)
        implicit real *8 (a-h,o-z)
        dimension y(n),us(n,k),ss(*),ysmall(k)

c
c        compute the residual of the least squares solution
c

        rr = 0

        do i=1,k
        ys = 0
        do j=1,n
        ys = ys + y(j)*us(j,i)
        enddo
        ysmall(i) = ys
        rr = rr + ys*ys
        enddo

        yy = 0
        do i=1,n
        yy = yy + (y(i))**2
        enddo

        resid = yy-rr

        return
        end
c
c
c
c
c
        subroutine rescale_a(t,a,asca,n,k,k1,k2)
        implicit real *8 (a-h,o-z)
        dimension a(n,k),asca(n,k)

c
c        multiply a by diagonal matrix
c

        tc = cos(t)
        ts = sin(t)

        do i=1,n

        do j=1,k1
        asca(i,j) = a(i,j)*tc
        enddo

        do j=(k1+1),k
        asca(i,j) = a(i,j)*ts
        enddo

        enddo

        return
        end
c
c
c
c
c
        subroutine get_rho0(ci, n, rho0)
        implicit real *8 (a-h,o-z)
c
c        find the lower bound of integration of the inner integral
c        using bisection and the formula for the maximum of the 
c        integrand
c

c
c        compute maximum of integrand
c
        dn = n - 2
        rho_max = sqrt(1/2.0*(sqrt(4*ci + dn**2) - dn))
        call eval_logdens_rho(rho_max,n,ci,fmax)
        fmax_log10 = fmax / log(10.0d0)

c
c        bisection
c        
        xmin = 1.0d-16
        xmax = rho_max 
        do 110 i=1,10
        rho0 = (xmin + xmax)/2.0d0
        call eval_logdens_rho(rho0,n,ci,f)
        f_log10 = f / log(10.0d0)
        if (f_log10 .lt. fmax_log10-18) xmin = rho0
        if (f_log10 .ge. fmax_log10-18) xmax = rho0
 110    continue

ccc        call eval_logdens_rho(rho0, n, ci, f)
ccc        call prin2('f at left endpoint*', exp(f), 1)

        return
        end
c
c
c
c
c
        subroutine get_rho1(ci, n, rho1)
        implicit real *8 (a-h,o-z)

c
c        find the upper bound of integration of the inner integral
c        using bisection and the formula for the maximum of the 
c        integrand and the second derivative at the maximum
c

c
c        compute maximum of integrand and second derivative at maximum
c

        dn = n - 2
        rho_max = sqrt(1/2.0*(sqrt(4*ci + dn**2) - dn))
        call eval_logdens_rho(rho_max,n,ci,fmax)
        fmax_log10 = fmax / log(10.0d0)

        dder2 = -1-3*ci/rho_max**4 + (n-2)/rho_max**2
        sd = 1/sqrt(-dder2)

c
c        bisection
c        
        xmin = rho_max
        xmax = rho_max + 40*sd
        do 120 i=1,10
        rho1 = (xmin + xmax)/2.0d0
        call eval_logdens_rho(rho1,n,ci,f)

        f_log10 = f / log(10.0d0)
        if (f_log10 .lt. fmax_log10-18) xmax = rho1
        if (f_log10 .ge. fmax_log10-18) xmin = rho1
 120    continue

ccc        call prin2('rho1*', rho1, 1)
ccc        call eval_logdens_rho(rho1, n, ci, f)
ccc        if (f .gt. -20) then
ccc          call prin2('bad!*', 1.0d0, 1)
ccc          call prin2('f*', f, 1)
ccc        endif

        return
        end
c
c
c
c
c
        subroutine plot_heatmap(nn1, nn2, fs)
        implicit real *8 (a-h,o-z)
        real*8 fs(1)
        data in/11/

        tol = -20.0d0

        call max_vec(nn1*nn2, fs, f1, ind)
ccc        call prin2('max*', f1, 1)

        do 320 i=1,nn1*nn2
        fs(i) = fs(i) - f1
        if (fs(i) .lt. tol) fs(i) = tol
ccc        call prin2('fs*', fs, nnodes**2)
 320    continue

c
c        if this function gets called twice, nothing gets
c        plotted, so change iw
c
        iw = in
        call pyimage(iw,nn1,nn2,fs,'title* ')

        in = in + 1

        return
        end
c
c
c
c
c       
        subroutine cc_arr(n, x, y)
        implicit real*8 (a-h,o-z)
        real*8 x(1),y(1)
        
        do 100 i=1,n
        y(i) = x(i)
 100    continue

        return 
        end
c
c
c
c
c
        subroutine mat_vec_mult(a, n, m, x, y)
        implicit real*8 (a-h,o-z)
        real*8 a(n, 1), x(1), y(1)


        do 450 i = 1, n
        y(i) = 0
        do 425 j = 1, m
        y(i) = y(i) + x(j) * a(i, j)
 425    continue
 450    continue

        
        return
        end
c
c
c
c
c
        subroutine mat_mult(a, n, m, b, k, c)
        implicit real*8 (a-h,o-z)
        real*8 a(n, 1), x(1), y(1), b(m, 1), c(n, 1)

        do 1075 i=1,n
        do 1050 j=1,k
        c(i,j) = 0
        do 1025 kk = 1,m
        c(i,j) = c(i,j) + a(i, kk) * b(kk, j)
 1025   continue
 1050   continue
 1075   continue
        
        return
        end
c
c
c
c
c
        subroutine inner_prod(u, v, n, y)
        implicit real *8 (a-h,o-z)
        real*8 u(1), v(1)

        y = 0
        do 1050 i=1,n
        y = y + u(i) * v(i)
 1050   continue

        return
        end
c
c
c
c
c
        subroutine mat_trans(u, n, m, ut)
        implicit real *8 (a-h,o-z)
        real*8 u(n, 1), ut(m, 1)

        do 1075 i=1,n
        do 1070 j=1,m
        ut(j, i) = u(i, j)
 1070   continue
 1075   continue

        return
        end
c
c
c
c
c
        subroutine theta_lege_nodes_whts(nn, t0, t1, ts, whts)
        implicit real *8 (a-h,o-z)
        real*8 ts(*), whts(*), u(1), v(1), x(nn+10),
     1      xleg(1000),wleg(1000)
        data nloc/0/
        save

        if (nn .ne. nloc) then
        nloc   = nn
        ifsave = 0

        itype=1
        call legeexps(itype,nn,xleg,u,v,wleg)
        endif

        do 60 i=1,nn
        tmp = (xleg(i) + 1)/2.0d0
        ts(i) = t0 + (t1-t0)*tmp
        whts(i) = wleg(i)*(t1-t0)/2.0d0
 60     continue

ccc        call prin2('wleg*', wleg, nn)
ccc        call prin2('xleg*', xleg, nn)

        return
        end
c
c
c
c
c
        subroutine lege_nodes_whts(nn, t0, t1, ts, whts)
        implicit real *8 (a-h,o-z)
        real*8 ts(*), whts(*), u(1), v(1), x(nn+10),
     1      xleg(1000),wleg(1000)
        data nloc/0/
        save

        if (nn .ne. nloc) then
        nloc   = nn
        ifsave = 0

        itype=1
        call legeexps(itype,nn,xleg,u,v,wleg)
        endif

        do 60 i=1,nn
        tmp = (xleg(i) + 1)/2.0d0
        ts(i) = t0 + (t1-t0)*tmp
        whts(i) = wleg(i)*(t1-t0)/2.0d0
 60     continue

ccc        call prin2('wleg*', wleg, nn)
ccc        call prin2('xleg*', xleg, nn)

        return
        end
c
c
c
c
c
        subroutine max_vec(n, v, f1, ind)
        implicit real *8 (a-h,o-z)
        real*8 v(1)

        f1 = v(1)
        ind = 1
        do 100 i=1,n-1
        if (v(i+1) .gt. f1) then
          f1 = v(i+1)
          ind = i+1
        endif
 100    continue

        return
        end
c
c
c
c
c
        subroutine pad_mat(b,n,m,nn,c)
        implicit real *8 (a-h,o-z)
        dimension b(n,m),c(nn,m)

        do i=1,nn
        do j=1,m
        c(i,j) = 0
        enddo 
        enddo

        do i=1,n
        do j=1,m
        c(i,j) = b(i,j)
        enddo
        enddo
        
        return
        end
c
c
c
c
c
        subroutine read_params(csv_file, n, k1, k2, a, y)
        implicit real *8 (a-h,o-z)
        real*8 a(*), y(*)
        character(*) csv_file
c
c        read the csv file written in rstudio by test43.R
c
ccc        csv_file = 'dense_params.csv'
        open (2, file=csv_file)
c
c        read matrix size
c
 100    format(i12)
        read(2, 100) n
ccc        call prinf('n*', n, 1)
        read(2, 100) k1
ccc        call prinf('k*', k, 1)
        read(2, 100) k2
ccc        call prinf('k*', k, 1)
c
c        read a
c
        k = k1+k2
ccc        call prinf('nk*', n*k, 1)
 200    format(f12.6)
        do 250 i=1,n*k
        read(2, 200) a(i)
ccc        call prinf('i*', i, 1)
ccc        call prin2('a(i)*', a(i), 1)
 250    continue

ccc        call prin2('a*', a, n*k)
c
c        read y
c
        do 350 i=1,n
        read(2, 200) y(i)
 350    continue

ccc        call prin2('y*', y, n)
        return
        end
c
c
c
c
c
        subroutine max_vec_abs(n, v, f1, ind)
        implicit real *8 (a-h,o-z)
        real*8 v(1)

        f1 = abs(v(1))
        ind = 1
        do 100 i=1,n-1
        if (abs(v(i+1)) .gt. f1) then
          f1 = abs(v(i+1))
          ind = i+1
        endif
 100    continue

        return
        end
c
c
c
c
c
        subroutine mat_mult_dense_diag(a, n, k, s, c)
        implicit real*8 (a-h,o-z)
        real*8 a(n, *), s(*), c(n, *)

        do 1075 i=1,n
        do 1050 j=1,k
        c(i,j) = s(j)*a(i,j)
 1050   continue
 1075   continue
        
        return
        end
c
c
c
c
c
        subroutine mat_mult_diag_dense(s, n, a, c)
        implicit real*8 (a-h,o-z)
        real*8 a(n, 1), s(1), c(n, 1)

        do 1075 i=1,n
        do 1050 j=1,n
        c(i,j) = s(i)*a(i,j)
 1050   continue
 1075   continue
        
        return
        end
c
c
c
c
c
        subroutine compute_ata(a,k,aout)
        implicit real *8 (a-h,o-z)
        dimension a(k,k),aout(k,k)

        do i=1,k
        do j=1,k
        z = 0
        do ijk=1,k
        z = z + a(ijk,i)*a(ijk,j)
        enddo
        aout(i,j) = z
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine sym_rescale_a(t,a,asca,k,k1,k2)
        implicit real *8 (a-h,o-z)
        dimension a(k,k),asca(k,k)

        tc = cos(t)
        ts = sin(t)

        tcc= tc*tc
        tss= ts*ts
        tcs= tc*ts

        do i=1,k1
        do j=1,k1
        asca(i,j) = a(i,j)*tcc
        enddo
        enddo

        do i=k1+1,k
        do j=k1+1,k
        asca(i,j) = a(i,j)*tss
        enddo
        enddo

        do i=1,k1
        do j=k1+1,k
        tt = a(i,j)*tcs
        asca(i,j) = tt
        asca(j,i) = tt
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine read_means(csv_file, m, dexps)
        implicit real *8 (a-h,o-z)
        real*8 dexps(*)
        character(100) csv_file

c
c        read the posterior means computed by stan from 
c        text file
c

        open (2, file = csv_file)

c
c        read posterior expectations
c

 200    format(f12.6)

        do 250 i=1,m+3
        read(2, 200) dexps(i)
 250    continue

        return
        end
c
c
c
c
c
        subroutine dd_abs_max(v1, v2, n, dd_max)
        implicit real*8 (a-h,o-z)
        real*8 v1(*), v2(*)

        dd_max = 0.0d0
        do i=1,n
        dd = abs(v1(i) - v2(i))
        if (dd .gt. dd_max) dd_max = dd
        enddo

        return
        end
c
c
c
c
c
        subroutine dd_abs_max_rel(v1, v2, n, dd_max)
        implicit real*8 (a-h,o-z)
        real*8 v1(*), v2(*)

        dd_max = 0.0d0
        do i=1,n
        dd = abs((v1(i) - v2(i))/v1(i))
        if (dd .gt. dd_max) dd_max = dd
        enddo

        return
        end
c
c
c
c
c
        subroutine write_exps_stds(file_out, k, dsums, stds)
        implicit real*8 (a-h,o-z)
        real*8 dsums(*), stds(*)
        character(*) file_out

        open (2, file=file_out)

 210    format(f22.16,','f22.16)
        do i=1,k+3
        write(2, 210) dsums(i), stds(i)
        enddo

        return
        end
c
c
c
c
c
        subroutine read_exps_stds(file_in, k, dsums, stds)
        implicit real*8 (a-h,o-z)
        real*8 dsums(*), stds(*)
        character(*) file_in

        open (2, file=file_in)

 210    format(f22.16,X,f22.16)
        do i=1,k+3
        read(2, 210) dsums(i), stds(i)
        enddo

        return
        end
c
c
c
c
c
        subroutine set_zero(n, v)
        implicit real *8 (a-h,o-z)
        real*8 v(*)

        do i=1,n
        v(i) = 0.0d0
        enddo

        return
        end
c
c
c
c
c
        subroutine div_arr(n, alph, v)
        implicit real *8 (a-h,o-z)
        real*8 v(*)

        do i=1,n
        v(i) = v(i) / alph
        enddo

        return
        end
c
c
c
c
c
        subroutine write_dd_exps_stds(file_out, k, dd, dsums, stds)
        implicit real*8 (a-h,o-z)
        real*8 dsums(*), stds(*)
        character(*) file_out

ccc        csv_file = 'expectations.dat'
        open (2, file=file_out)

 210    format(f22.16,','f22.16)

        write(2, 210) dd

        do i=1,k+3
        write(2, 210) dsums(i), stds(i)
        enddo

        return
        end
