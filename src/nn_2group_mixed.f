c
c
c
c
c
        subroutine mixed_test_integration()
        implicit real *8 (a-h,o-z)
        real*8 a(150 000 000), y(500 000), ss(100 000)
        real *8, allocatable :: dsums(:), dexps(:),
     1     dsums2(:), dds(:), stds(:), stds2(:), cov(:, :)
        character(100) csv_file, filename

c
c        parameters
c
        filename = 'params.dat'
        call mixed_read_params(filename, n, k1, k2, a, y, ss)
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
        allocate(cov(k, k))

c
c        integrate
c
        nn = 80
        nnt = 40
        call cpu_time(t1)
        call mixed_2group(nnt, nn, n, k1, k2, k, a, y, ss, 
     1     dsums, dsum, stds, cov)
        call cpu_time(t2)
        call prin2('dsums*', dsums, k+2)
        call prin2('stds*', stds, k+2)
        call prin2('total time*', t2-t1, 1)

c
c        double nodes in all directions
c
        nn2 = 2*nn
        nnt2 = 2*nnt
        call mixed_2group(nnt2, nn2, n, k1, k2, k, a, y, ss,
     1     dsums2, dsum2, stds2, cov)
ccc        call prin2('dsums*', dsums2, k+2)
ccc        call prin2('stds*', stds2, k+2)

c
c        check error
c
        call mixed_dd_abs_max(dsums2, dsums, k+2, dd_max)
        call prin2('max posterior mean error*', dd_max, 1)

        call mixed_dd_abs_max(stds2, stds, k+2, dd_max)
        call prin2('max posterior std error*', dd_max, 1)

        filename = 'exps.dat'
        call mixed_write_exps_stds(filename, k, dsums, stds)
ccc        call prin2('dsums*', dsums,k+3)

c
c        print last k2 cols of covariance matrix
c
        filename = 'cov_cols.dat'
        call mixed_write_cov_cols(filename, k, k2, cov(1,k1+1))

c
c        print all dds
c
        do i=1,k+2
        dds(i) = dsums(i) - dsums2(i)
        enddo

        filename = 'dds_exps.dat'
        call mixed_write_dds(filename, k+2, dds)

        do i=1,k+2
        dds(i) = stds(i) - stds2(i)
        enddo

        filename = 'dds_stds.dat'
        call mixed_write_dds(filename, k+2, dds)

        stop

c
c        stan comparison
c
        csv_file = 'means.dat'
        call mixed_read_means(csv_file, k, dexps)
        call prin2('stan expectations*', dexps, k+2)

        call mixed_dd_abs_max(dexps, dsums, k+2, dd_max)
        call prin2('stan max error*', dd_max, 1)
        
        return
        end
c
c
c
c
c
        subroutine mixed_2group(nnt, nn, n, k1, k2, k, a, y, ss, 
     1     dsums, stds, dsums_cov)
        implicit real *8 (a-h,o-z)
        real*8 a(n,*), y(*), dsums(*), stds(*),
     1     ts(nnt),rhos(nn),phis(nn), s(k),s2(k), ss(k2),
     1     ysmall(k),ys2(k),vmoms(100 000), dsumsi(k+3+10),
     1     whts_ts(nnt+10),whts_phis(nn+10),whts_rhos(nn+10),
     1     dsum_xsi(k+2), dsum_vars(k+10), dsums_cov(k, k),
     1     s_inv(k), s0(k), x1(k), tmp_mat2(k**2), ax(n), tmp_arr(k)
        real *8, allocatable :: vt(:,:),v(:,:), b(:,:),ynew(:), 
     1     dsums_covi(:,:), xxti(:,:),tmp_mat(:,:), a2(:,:), a2t(:,:),
     1     ata(:,:), u(:,:), ut(:,:)
c
c        this subroutine computes posterior expectations of the two
c        group bayesian linear regression model
c
c          y ~ normal(A1*x1 + A2*x2, sig1)
c          x1 ~ normal(0, sig2)
c          x2 ~ normal(0, [ss(1), ss(2), ...,ss(k2)])
c          sig{1,2} ~ normal+(0, [d1, d2])
c
c        this subroutine computes E[x_i], var(x_i), E[sig{1,2}]. 
c

        done = 1.0d0
        pi = 4*atan(done)

        k = k1+k2

        allocate (ut(k,n))
        allocate (u(n,k))
        allocate (vt(k,k))
        allocate (v(k,k))
        allocate (b(k,k))
        allocate (ynew(k))
        allocate (dsums_covi(k,k))
        allocate (xxti(k,k))
        allocate (tmp_mat(k,k))
        allocate (a2(n,k))
        allocate (a2t(k,n))
        allocate (ata(k,k))

c
c        adjust for fixed priors on coefficients k1+1 to k1+k2
c
        call  mixed_cc_arr(n*k, a, a2)
        do j=1,k2
        do i=1,n
        a2(i, k1+j) = a2(i, k1+j) * ss(j)**2
        enddo
        enddo

c
c       hyperpriors -- variances, not standard deviations
c
        d1 = 1.0d0**2
        d2 = 1.0d0**2

c
c        before computing integral find least squares solution to
c        ax = y as well as projection of y onto the
c        columns of a
c

c
c        second method
c
        call mixed_mat_trans(a2, n, k, a2t)
        call dmatmat(k, n, a2t, k, a2, ata)
        call deigs(k, ata, v, s)
        call mixed_mat_trans(v, k, k, vt)

        do i=1,k
        s_inv(i) = 1/s(i)
        s0(i) = sqrt(s(i))
        enddo

c
c        solve linear system 
c        
        call mixed_mat_mult_diag_dense(s0, k, vt, b)
ccc        call prin2('b*', b, k**2)

        call mixed_mat_mult(a2t, k, n, y, 1, tmp_mat)
        call mixed_mat_mult(vt, k, k, tmp_mat, 1, tmp_mat2)
        call mixed_mat_mult_diag_dense(s_inv, k, tmp_mat2, tmp_mat)
        call mixed_mat_mult(v, k, k, tmp_mat, 1, x1)

        call mixed_mat_mult(vt, k, k, x1, 1, tmp_arr)
        do i=1,k
        ynew(i) = s0(i)*tmp_arr(i)
        enddo

c
c        residual 
c
        call mixed_mat_mult(a2, n, k, x1, 1, ax)
        resid = 0.0d0
        do i=1,n
        resid = resid + (ax(i) - y(i))**2
        enddo


c
c        construct grid for outer parameter of integration, theta
c
        t0 = 0.0
        t1 = pi/2.0d0
        call mixed_theta_lege_nodes_whts(nnt, t0, t1, ts, 
     1     whts_ts)

c
c        initialize sums (integrals) to be computed
c          dsums - expectations
c          dsums_cov - covariance
c          ss1 - E[sig1**2]
c          ss2 - E[sig2**2]
c          ss3 - E[sig3**2]
c          fm - scaling constant
c
        ktmp = k+2
        call mixed_set_zero(ktmp, dsums)
        dsum = 0
        call mixed_set_zero(k*k, dsums_cov)
        ss1 = 0
        ss2 = 0
        fm = -1.0d250
c
c        theta integral
c       
        do 130 i=1,nnt
ccc        t = ts(i)
ccc        wht_t = whts_ts(i)
        t = ts(nnt+1-i)
        wht_t = whts_ts(nnt+1-i)
ccc        call prin2('t*', t, 1)

c
c        compute phi integral
c
        call mixed_eval_inner(nn, n, k1, k2, k, d1, d2, b,
     1     t, resid, ynew, dsumsi,
     1     dsumi, ss1i, ss2i, dsums_covi, dsum_xsi, xxti, fmi)

c
c        due to underflow issues, integrate over theta by computing 
c        a sum of the form \sum_i exp(fi)*gi such that at the end
c        we have an expression exp(fm)*dsum 
c
        fi = fmi
        wt = wht_t
        if (fi .gt. fm) then
          dsum = dsum * exp(fm-fi) + dsumi*wt

          do ijk=1,k
          dsums(ijk) = dsums(ijk) * exp(fm-fi) + dsum_xsi(ijk)*wt
          enddo
          dsums(k+1) = dsums(k+1) * exp(fm-fi) + dsumsi(k+1)*wt
          dsums(k+2) = dsums(k+2) * exp(fm-fi) + dsumsi(k+2)*wt

          ss1 = ss1 * exp(fm-fi) + ss1i*wt
          ss2 = ss2 * exp(fm-fi) + ss2i*wt

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

          ss1 = ss1 + exp(fi-fm)*ss1i*wt
          ss2 = ss2 + exp(fi-fm)*ss2i*wt

          do i1=1,k
          do i2=1,k
          tmp = dsums_covi(i1,i2)+xxti(i1,i2)
          dsums_cov(i1,i2)=dsums_cov(i1,i2) + exp(fi-fm)*tmp*wt
          enddo
          enddo
  
        endif

ccccccc
        if (fi/log(10.0d0) .lt. fm/log(10.0d0) - 20.0d0) goto 140
ccccccc

 130    continue
 140    continue

c
c        scale first and second moments by normalizing constant
c
        do i=1,k+2
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
        do j=1,k
        dsums_cov(i, j)=dsums_cov(i, j) - dsums(i)*dsums(j)
        enddo
        enddo

c
c        get stds
c
        do i=1,k
        stds(i) = sqrt(dsums_cov(i, i))
        enddo

c
c        get variances of sig1, sig2, sig3
c
        ss1 = ss1/dsum
        ss2 = ss2/dsum
        stds(k+1) = sqrt(ss1 - dsums(k+1)**2)
        stds(k+2) = sqrt(ss2 - dsums(k+2)**2)

c
c        readjust for scale parameter priors 
c        
        do j=1,k2
        dsums(k1+j) = dsums(k1+j) * ss(j)**2
        stds(k1+j) = stds(k1+j) * ss(j)**2
        enddo

c
c        plot
c        
ccc        call prin2('fs*', fs, nnt*nn)
ccc        call mixed_plot_heatmap(nnt, nn, fs)


        return
        end
c
c
c
c
c
        subroutine mixed_eval_inner(nn, n, k1, k2, k, d1, d2, b,
     1     t, resid, ynew, dsumsi, 
     1     dsumi, ss1i, ss2i, dsums_covi, dsum_xsi, xxti, fmi)
        implicit real *8 (a-h,o-z)
        real*8 ysmall(k), ys2(k), ynew(*),phis(nn), s(k),s2(k), 
     1     xs(k), vmoms(k), dsumsi(k+2), dsums_covi(k,k),
     1     whts_phis(nn), dsum_xsi(k+2), dsum_vars(k+2),
     1     wwti(k,k), b(k,k), xxti(k,k), u(k, k)
        real *8, allocatable :: v(:,:), vt(:,:), c(:,:), ct(:,:), 
     1     tmp_mat(:,:), asca(:,:)


        done = 1.0d0
        pi = 4*atan(done)

        k = k1+k2

        allocate (vt(k,k))
        allocate (v(k,k))
        allocate (c(k,k))
        allocate (ct(k,k))
        allocate (tmp_mat(k,k))
        allocate (asca(k,k))

c
c        compute the diagonal form of a^t*a for each theta which
c        will be used to convert the integrand to a diagonal 
c        gaussian for each (theta, phi, rho) 
c
        call mixed_rescale_a(t,b,asca,k,k,k1,k2)
        call dsvd(k,k, asca, u, s, vt)
        call mixed_y_proj(ynew,u,s,k,k,ysmall,res2)
        call mixed_entry_sqr(ysmall,ys2,k)
        call mixed_entry_sqr(s,s2,k)

c
c        initialize sums taken in inner loop
c          wwti is E[w*w^t] where w is in the coordinate 
c          system that depends on theta
c
        dsumi = 0
        ss1i = 0
        ss2i = 0
        ktmp = k+3
        call mixed_set_zero(ktmp, dsumsi)
        call mixed_set_zero(k**2, wwti)
        call mixed_set_zero(k, dsum_vars)
        fmi = -1.0d250

c
c        for each theta, get upper integration bound for phi
c
        call mixed_get_phi(nn, n, k, t, ysmall, ys2, s, s2, 
     1     resid, d1, d2, phi0i, phi1i, fmax_theta)
        call mixed_lege_nodes_whts(nn, phi0i, phi1i, phis, whts_phis)
ccc        call prin2('phis new*', phis, nn)

c
c          phi integral 
c
        do 120 j=1,nn
        phi = phis(j)
        wht_phi = whts_phis(j)
        call mixed_get_beta_alpha_prefact(phi,ysmall,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)
        call mixed_get_mjs(k,s2,ysmall,phi,vmoms)

c
c        compute rho as a function of theta and phi
c
        sig22 = 1/tan(t)**2
        sig12 = -cos(phi)**2*(sig22+1)/(cos(phi)**2 - 1)
        rho = sqrt(sig12 + sig22 + 1)

        a1 = cos(phi)**2/d1 + sin(phi)**2*cos(t)**2/d2
        call mixed_eval_logdens_rho(n, rho, a1, exp_fact, f)

c
c        compute determinant of jacobian
c
        call mixed_eval_jac_det(rho, phi, t, djac)

c
c        for fixed theta, compute integral over phi by a sum of 
c        the form \sum_i exp(fi)*gi so that
c        we have an expression exp(fmi)*dsum, the reason we do 
c        this is that exp(fi) is often smaller than 10^{-250}
c       

        fi = prefact + f - log(abs(djac))
        wt = wht_phi
        gi = 1.0d0
        gi1 = rho
        gi2 = rho**2
        wt = wht_phi
        if (fi .gt. fmi) then
          dsumi = dsumi * exp(fmi-fi) + gi*wt

          do ijk=1,k
          dsumsi(ijk) = dsumsi(ijk) * exp(fmi - fi) + gi*wt*vmoms(ijk) 
          enddo        

          dsumsi(k+1) = dsumsi(k+1)*exp(fmi-fi) + gi1*cos(phi)*wt
          dsumsi(k+2) = dsumsi(k+2)*exp(fmi-fi) + gi1*sin(phi)*cos(t)*wt

          ss1i = ss1i*exp(fmi-fi) + gi2*cos(phi)**2*wt
          ss2i = ss2i*exp(fmi-fi) + gi2*sin(phi)**2*cos(t)**2*wt

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

          ss1i = ss1i + exp(fi-fmi)*gi2*cos(phi)**2*wt
          ss2i = ss2i + exp(fi-fmi)*gi2*sin(phi)**2*cos(t)**2*wt

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
        call mixed_mat_trans(vt, k, k, v)
        call mixed_get_xs_to_ws_matrix(k,k1,k2,v,t,c,ct)
        call mixed_mat_mult_dense_diag(c, k, k, dsum_vars, tmp_mat)
        call mixed_mat_mult(tmp_mat, k, k, ct, k, dsums_covi)

c
c        recover matrix E[x*x^t]
c
        call mixed_mat_mult(c, k, k, wwti, k, tmp_mat)
        call mixed_mat_mult(tmp_mat, k, k, ct, k, xxti)
        call mixed_get_xs_from_ws(dsumsi,k,k1,k2,vt,t, dsum_xsi)

        return
        end
c
c
c
c
c
        subroutine mixed_get_phi(nn, n, k, t, ysmall, ys2, s, s2, 
     1     resid, d1, d2, phi0i, phi1i, fmax)
        implicit real *8 (a-h,o-z)
        real*8 phis(nn),whts_phis(nn), ysmall(*),ys2(*),s(*),s2(*),
     1     fs(nn), tmp_mat(nn)

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
        call mixed_lege_nodes_whts(nn, phi0, phi1, phis, whts_phis)

c
c        tabulate function
c
        fmax = -1.0d250
        do 119 j=1,nn
        phi = phis(j)
        call mixed_get_beta_alpha_prefact(phi,ysmall,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)
        a1 = cos(phi)**2/d1
        a1 = a1 + sin(phi)**2*cos(t)**2/d2
c
c        compute rho as a function of theta and phi
c
        sig22 = 1/tan(t)**2
        sig12 = -cos(phi)**2*(sig22+1)/(cos(phi)**2 - 1)
        rho = sqrt(sig12 + sig22 + 1)

        call mixed_eval_jac_det(rho, phi, t, djac)
        call mixed_eval_logdens_rho(n, rho, a1, exp_fact, f)

        fi = prefact + f - log(djac)
ccc        call prin2('fi*', fi, 1)
        fs(j) = fi

 119    continue

ccc        call prin2('fs*', fs, nn)


c       
c        get max 
c
        call mixed_max_vec(nn, fs, fmax, ind)

        do i=1,nn
        tmp_mat(i) = 1
        if (fs(i) .gt. (fmax - 16)) tmp_mat(i) = 0
        enddo
        call mixed_get_int_bds(nn, tmp_mat, i0, i1)

        phi0i = phis(i0)
        phi1i = phis(i1)

c
c        when using legendre nodes, choose endpoint instead of first 
c         node
c       
        if (i0 .eq. 1) phi0i = phi0
        if (i1 .eq. nn) phi0i = phi1

        return
        end
c
c
c
c
c
        subroutine mixed_eval_jac_det(rho, phi, t, f)
        implicit real *8 (a-h,o-z)
c
c        compute the determinant of the change of variables
c        jacobian for the mixed effects model
c

        sig1 = rho*cos(phi)
        sig2 = rho*sin(phi)*cos(t)
        x1 = sig1
        y1 = sig2

        alph = x1**2+y1**2+1
        djac = 1/(1+y1**2)
        djac = djac*(1/sqrt(alph)-x1**2/alph**(1.5))/sqrt(1-x1**2/alph)
        f = abs(djac)

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
c   mixed_2group - evaluates two-group normal normal model with 
c         normal+(0, 1) hyperpriors on the scale parameters 
c 
c
c
c
c
c
c
c
c
c
        subroutine mixed_get_xs_to_ws_matrix(k,k1,k2,v,t,a,at)
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

        call mixed_mat_mult_diag_dense(s1, k, v, a)
        call mixed_mat_trans(a, k, k, at)

        return
        end
c
c
c
c
c
        subroutine mixed_get_mjs(k,s2,ys,phi,vmoms)
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
        subroutine mixed_eval_logdens_rho(n, rho, a1, c1, f)
        implicit real *8 (a-h,o-z)

        f = -n*log(rho) - rho**2/2.0d0*a1 - c1/(2*rho**2)

        return
        end
c
c
c
c
c
        subroutine mixed_get_beta_alpha_prefact(phi,ys,ys2,s,s2,
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
        prefact = -(n-k)/2.0d0*log(cp)+alpha
        exp_fact= resid/cp+beta

        return
        end
c
c
c
c
c
        subroutine mixed_get_xs_from_ws(ws,k,k1,k2,vt,t, xs)
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
        subroutine mixed_entry_sqr(xin,xout,k)
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
        subroutine mixed_y_proj(y,us,ss,n,k,ysmall,resid)
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
        subroutine mixed_rescale_a(t,a,asca,n,k,k1,k2)
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
        subroutine mixed_plot_heatmap(nn1, nn2, fs)
        implicit real *8 (a-h,o-z)
        real*8 fs(1)
        data in/11/

        tol = -20.0d0

        call mixed_max_vec(nn1*nn2, fs, f1, ind)
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
        subroutine mixed_cc_arr(n, x, y)
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
        subroutine mixed_mat_vec_mult(a, n, m, x, y)
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
        subroutine mixed_mat_mult(a, n, m, b, k, c)
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
        subroutine mixed_inner_prod(u, v, n, y)
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
        subroutine mixed_mat_trans(u, n, m, ut)
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
        subroutine mixed_theta_lege_nodes_whts(nn, t0, t1, ts, whts)
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
        subroutine mixed_lege_nodes_whts(nn, t0, t1, ts, whts)
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
        subroutine mixed_max_vec(n, v, f1, ind)
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
        subroutine mixed_read_params(csv_file, n, k1, k2, a, y, ss)
        implicit real *8 (a-h,o-z)
        real*8 a(*), y(*), ss(*)
        character(*) csv_file
c
c        read the csv file written in rstudio by test43.R
c
        open (2, file=csv_file)
c
c        read matrix size
c
 100    format(i12)
        read(2, 100) n
        read(2, 100) k1
        read(2, 100) k2
c
c        read a
c
        k = k1+k2
 200    format(f12.6)
        do 250 i=1,n*k
        read(2, 200) a(i)
 250    continue

c
c        read y
c
        do 350 i=1,n
        read(2, 200) y(i)
 350    continue

c
c        read ss
c
        do i=1,k2
        read(2, 200) ss(i)
        enddo

c
c        scale hyperpriors
c
ccc        read(2, 200) sigy
ccc        read(2, 200) sig2


        return
        end
c
c
c
c
c
        subroutine mixed_max_vec_abs(n, v, f1, ind)
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
        subroutine mixed_mat_mult_dense_diag(a, n, k, s, c)
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
        subroutine mixed_mat_mult_diag_dense(s, n, a, c)
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
        subroutine mixed_read_means(csv_file, m, dexps)
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

        do 250 i=1,m+2
        read(2, 200) dexps(i)
 250    continue

        return
        end
c
c
c
c
c
        subroutine mixed_dd_abs_max(v1, v2, n, dd_max)
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
        subroutine mixed_write_exps_stds(file_out, k, dsums, stds)
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
        subroutine mixed_read_exps_stds(file_in, k, dsums, stds)
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
        subroutine mixed_set_zero(n, v)
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
        subroutine mixed_get_int_bds(n, v, i0, i1)
        implicit real *8 (a-h,o-z)
        real*8 v(n)

c
c        this subroutine takes as input an array of the form
c        (1, 1, ..., 1, 0, ..., 0, 1, 1, ..., 1) and returns
c        the indices of the last 1 before the 0s and the first
c        1 after the 0s
c

c
c        if first and last element is 0, return endpoints
c
        if ((v(1) .lt. 0.5d0) .and. (v(n) .lt. 0.5d0)) then
          i0 = 1
          i1 = n
          goto 900
        endif
c
c        if first element is 0, go until we get a 1
c
        i0 = 1
        i1 = n
        if (v(1) .lt. 0.5d0) then
          i0 = 1
          do i=2,n
          if (v(i) .gt. 0.5d0) then
            i1 = i
            goto 900
          endif
          enddo
        endif

c
c        if last element is 0, go from end until we get 1
c
        if (v(n) .lt. 0.5d0) then
          i1 = n
          do i=1,n-1
          if (v(n-i) .gt. 0.5d0) then
            i0 = n-i
            goto 900
          endif
          enddo
        endif

c
c        otherwise find the last 1 before the 0's and first 
c        1 after 0's
c

        do i=1,n
        if (v(i) .lt. 0.5d0) then
          i0 = i-1
          goto 800
        endif
        enddo

 800    continue

        do i=1,n
        ind = n+1-i
ccc        call prin2('v(ind)*', v(ind), 1)
ccc        call prinf('ind*', ind, 1)
        if (v(ind) .lt. 0.5d0) then
          i1 = ind+1
          goto 900
        endif
        enddo


 900    continue

        return
        end
c
c
c
c
c
        subroutine mixed_write_cov_cols(file_out, k, k1, cov)
        implicit real*8 (a-h,o-z)
        real*8 cov(*)
        character(*) file_out

        open (2, file=file_out)
 210    format(f22.16)

        m = k*k1
        do i=1,m
        write(2, 210) cov(i)
        enddo

        return
        end
c
c
c
c
c
        subroutine mixed_write_dds(file_out, k, dds)
        implicit real*8 (a-h,o-z)
        real*8 dds(*)
        character(*) file_out

        open (2, file=file_out)
 210    format(f22.16)

        do i=1,k
        write(2, 210) dds(i)
        enddo

        return
        end
