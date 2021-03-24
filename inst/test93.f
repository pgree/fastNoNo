c
c
c
c
c
        subroutine testit(n)
        implicit real*8 (a-h,o-z)

        call test_integration()
ccc        call utest()
ccc        call utest2()

        return
        end
c
c
c
c
c
        subroutine test_integration()
        implicit real *8 (a-h,o-z)
        real*8 a(10 000 000), x(100 000), y(1000 000), dexps(10 000),
     1     dsums2(10 000), dsums(100 000), dds(10 000), stds(10 000),
     1     stds2(10 000)
        character(100) csv_file, filename

c
c        parameters
c
        csv_file = 'test_dense_params2.dat'
        csv_file = 'test_dense_params.dat'
        csv_file = 'abort_params.dat'
        call read_params(csv_file, n, k1, k2, a, y)
ccc        call prin2('a = *',a,n*k)
ccc        call prin2('y = *',y,n)
        k = k1+k2
        call prinf('n*', n, 1)
        call prinf('k1*', k1, 1)
        call prinf('k2*', k2, 1)

c
c        integrate
c
        nn = 60
        nn_theta = 20
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

        stop

c
c        stan comparison
c
        csv_file = 'means.dat'
        csv_file = 'test_means_dense_stan2.dat'
        csv_file = 'test_means_dense_stan.dat'
        call read_means(csv_file, k, dexps)
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
        subroutine utest2()
        implicit real *8 (a-h,o-z)
        real*8 a(10 000 000), x(100 000), y(1000 000), dexps(10 000),
     1     dsums2(10 000), dsums(100 000), dds(10 000), stds(10 000),
     1     stds2(10 000)
        character(100) csv_file, filename

c
c        parameters
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
c        parameters
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

        filename = 'test_means_dense_fast.dat'
        call read_exps_stds(filename, k, dsums2, stds2)

        call dd_abs_max(dsums2, dsums, k+3, dd_max)
        call prin2('max posterior mean error*', dd_max, 1)

        call dd_abs_max(stds2, stds, k+3, dd_max)
        call prin2('max posterior std error*', dd_max, 1)

c
c        stan comparison
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
c
c
        subroutine dense_eval(nn_theta, nn, n, k1, k2, a, y, dsums, 
     1     dsum, stds)
        implicit real *8 (a-h,o-z)
        real*8 a(*), y(*), dsums(*),
     1     thetas(nn_theta+10),rhos(nn+10),phis(nn+10),phis0(nn+10),
     1     s(100 000),s2(100 000), whts_phis0(nn+10),
     1     ysmall(100 000),ys2(100 000),xs(k1+k2),
     1     vmoms(100 000),ys(100 000),dsumsi(k1+k2+3+10),
     1     whts_thetas(nn_theta+10),whts_phis(nn+10),whts_rhos(nn+10),
     1     dsum_xsi(k1+k2+10), dsum_vars(10 000),
     1     ztilde(10 000),xtilde(10 000), prefacts(10 000),
     1     stds(10 000), c(10 000), c_inv(10 000), tmp_mat(10 000)
        real *8, allocatable :: u(:,:),vt(:,:),asca(:,:),ut(:,:),
     1     b(:,:),ynew(:),v(:,:),bb(:,:), fmaxs(:), fs(:,:),
     1     dsums_cov(:,:), dsums_covi(:,:), wwti(:,:), xxti(:,:),
     1     gs(:,:)

c
c        this subroutine computes E[x_i], var(x_i), a
c

        done = 1.0d0
        pi = 4*atan(done)

        k = k1+k2

        allocate (u(n,k))
        allocate (ut(n,k))
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



        call dsvd(n,k, a, u, s, vt)
        call mat_trans(u, n, k, ut)
        call mat_mult(ut, k, n, a, k, b)
        call y_proj(y,u,s,n,k,ysmall,resid)
        call mat_vec_mult(ut, k, n, y, ynew)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc         .   .   .   for eigenvalue version

        call mat_trans(u, n, k, ut)
        call compute_ata(b,k,bb)
        do ii=1,k
        ztilde(ii) = ynew(ii)/s(ii)
        enddo
        call mat_trans(vt, k, k, v)
        call mat_vec_mult(v, k, k, ztilde, xtilde)

cccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c        construct grid
c
        theta0 = 0.0
        theta1 = pi/2.0d0
        call theta_lege_nodes_whts(nn_theta, theta0, theta1, thetas, 
     1     whts_thetas)

c
c        initialize to 0
c

        do ijk=1,k+3
        dsums(ijk)=0
        enddo
        dsum = 0

        do ijk=1,k
        do ijk2=1,k
        dsums_cov(ijk,ijk2)=0
        enddo
        enddo

        fm = -1.0d250

c
c        theta integral
c       
        do 130 i=1,nn_theta
        t = thetas(i)
        wht_theta = whts_thetas(i)


        call rescale_a(t,b,asca,k,k,k1,k2)
        call dsvd(k,k, asca, u, s, vt)
        call y_proj(ynew,u,s,k,k,ysmall,res2)
        call entry_sqr(ysmall,ys2,k)
        call entry_sqr(s,s2,k)

        dsumi = 0
        do ijk=1,k+3
        dsumsi(ijk)=0
        enddo

        do i1=1,k
        do i2=1,k
        dsums_covi(i1,i2)=0.0d0
        wwti(i1, i2) = 0.0d0
        enddo
        enddo

        do ijk=1,k
        dsum_vars(ijk) = 0.0d0
        enddo

        fmi = -1.0d250
c
c          phi integral 
c
        call get_phi1_fmax(nn, n, k, t, ysmall, ys2, s, s2, resid,
     1     phi1i, fmax_theta)
        fmaxs(i) = fmax_theta
ccc        call prin2('phi1i*', phi1i, 1)
ccc        call prin2('fmax_theta*', fmax_theta, 1)

        phi0 = 0.0
        call lege_nodes_whts(nn, phi0, phi1i, phis, whts_phis)
ccc        call prin2('phis new*', phis, nn)

        do 120 j=1,nn
        phi = phis(j)
        call compute_beta_alpha_prefact2(phi,ysmall,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)
        call compute_mjs(k,s2,ysmall,phi,vmoms)

        wht_phi = whts_phis(j)

c
c          rho integral 
c
        call eval_rho_int2(nn, n, k, exp_fact, 
     1     dsum_rho, dsum_rho1, dsum_rho2, fmax)

c
c        compute sum of the form exp(fi)*gi such that at the end
c        we have an expression exp(fmi)*dsum 
c

ccc        call prin2('dsumi*', dsumi, 1)
ccc        call prin2('dsum_rho*', dsum_rho, 1)
        fi = prefact + fmax - fmaxs(1)
        gi = dsum_rho
        gi1 = dsum_rho1
        gi2 = dsum_rho2
        wt = wht_phi
ccc        fs(i, j) = log10(gi * exp(fi))
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

        call mat_trans(vt, k, k, v)
        call get_xs_to_ws_matrix(k,k1,k2,v,vt,t,c,c_inv)
        call mat_mult_dense_diag(c, k, k, dsum_vars, tmp_mat)
        call mat_mult(tmp_mat, k, k, c_inv, k, dsums_covi)

        call mat_mult(c, k, k, wwti, k, tmp_mat)
        call mat_mult(tmp_mat, k, k, c_inv, k, xxti)

        call get_xs_from_ws(dsum_xsi,dsumsi,k,k1,k2,vt,t)
ccc        call prin2('xs 2*', dsums_xsi, k)
ccc        call prin2('dsumi*', dsumi, 1)

c
c        compute sum of the form exp(fi)*gi such that at the end
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
c        compute error with 
c
ccc        call get_err(nn_theta, nn, fs, gs)

c
c        scale by normalizing constant
c
        call scaling(k, dsum, dsums, dsums_cov, stds)
ccc        call prin2('dsums*', dsums, k)
ccc        call prin2('stds*', stds, k)

c
c        plot
c        
ccc        call prin2('fs*', fs, nn_theta*nn)
        call plot_heatmap(nn_theta, nn, fs)


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
        subroutine scaling(k, dsum, dsums, dsums_cov, stds)
        implicit real *8 (a-h,o-z)
        dimension dsums(*), dsums_cov(k, *), stds(*)

        do i=1,k+3
        dsums(i)=dsums(i)/dsum
        enddo

        do i=1,k
        do j=1,k
        dsums_cov(i, j)=dsums_cov(i, j)/dsum
        enddo
        enddo

        do i=1,k
        stds(i) = sqrt(dsums_cov(i, i) - dsums(i)**2)
        enddo

        return
        end
c
c
c
c
c
        subroutine get_xs_to_ws_matrix(k,k1,k2,v,vt,t,a,a_inv)
        implicit real *8 (a-h,o-z)
        dimension s1(10 000),s1_inv(10 000),v(k,k),a(*),a_inv(*)

c
c        construct diagonal 
c
        do i=1,k1
        s1(i) = cos(t)
ccc        s1(i) = sqrt(cos(t))
        s1_inv(i) = 1/s1(i)
        enddo
        
        do i=1,k2
        s1(i+k1) = sin(t)
ccc        s1(i+k1) = sqrt(sin(t))
        s1_inv(i+k1) = 1/s1(i+k1)
        enddo


        call mat_mult_diag_dense(s1, k, v, a)
        call mat_mult_dense_diag(vt, k, k, s1_inv, a_inv)
ccc        call mat_mult_diag_dense(s1_inv, k, vt, a)
ccc        call mat_mult_dense_diag(v, k, k, s1, a_inv)

ccc        call prin2('s1*', s1, k)
ccc        call prin2('a*', a, k*k)

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
        call mat_trans(a, k, k, a_inv)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc


        return
        end
c
c
c
c
c
        subroutine theta_int1_eigs(t, nn, n, k1, k2, bb, y, 
     1     resid, xtilde, dsumi, dsum_xsi, dsum1i, 
     1      dsum2i, dsum3i)
        implicit real *8 (a-h,o-z)
        real*8 y(*), bb(*),dsum_xsi(*),
     1     rhos(nn+10),phis(nn+10), s(100 000),s2(100 000), 
     1     ys2(100 000), vmoms(100 000),ys(100 000),dsumsi(k1+k2+10),
     1     whts_thetas(nn+10),whts_phis(nn+10),whts_rhos(nn+10),
     1     xtilde(*),wtilde(10 000),ztilde(10 000),ysmall(10 000)
        real *8, allocatable :: u(:,:),vt(:,:),asca(:,:),ut(:,:),
     1     fmaxs(:,:),v(:,:)


        done = 1.0d0
        pi = 4*atan(done)

        k = k1+k2

        allocate (u(n,k))
        allocate (ut(k,n))
        allocate (vt(k,k))
        allocate (v(k,k))
        allocate (asca(n,k))
        allocate (fmaxs(nn,nn))

        do ii=1,k1
        wtilde(ii) = xtilde(ii)/cos(t)
        enddo
        do ii=k1+1,k
        wtilde(ii) = xtilde(ii)/sin(t)
        enddo

        call sym_rescale_a(t,bb,asca,k,k1,k2)
        call deigs(k, asca, u, s2)
        do ii=1,k
        s(ii) = sqrt(s2(ii))
        enddo
        call mat_trans(u, k, k, vt)
        call mat_vec_mult(vt, k, k, wtilde, ztilde)
        do ii=1,k
        wtilde(ii) = ztilde(ii)*s(ii)
        enddo
        do ii=1,k
        ysmall(ii) = wtilde(ii)
        enddo

        call entry_sqr(ysmall,ys2,k)

ccc        call prin2('s2 new*', s2, k)
ccc        call prin2('s new*', s, k)
ccc        call prin2('ysmall*', ysmall, k)
ccc        call prin2('ys2 new*', ys2, k)
ccc        call prin2('resid new*', resid, 1)

        dsumi = 0
        dsum1i = 0
        dsum2i = 0
        dsum3i = 0
        do ijk=1,k+3
        dsumsi(ijk)=0
        enddo

c
c        phi integral 
c
        call get_phi1(nn, n, k, t, ysmall, ys2, s, s2, resid, phi1i)
ccc        call prin2('phi1i new*', phi1i, 1)
        phi0 = 0.0
        call lege_nodes_whts(nn, phi0, phi1i, phis, whts_phis)
ccc        call prin2('phis new*', phis, nn)

        do 120 j=1,nn
        phi = phis(j)
        call compute_beta_alpha_prefact(phi,ysmall,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)
        call compute_mjs(k,s2,ysmall,phi,vmoms)

        wht_phi = whts_phis(j)

c
c        rho integral 
c
        call eval_rho_int(nn, n, k, prefact, exp_fact, 
     1     t, phi, rhos, whts_rhos, 
     1     dsum1_rho, dsum2_rho, dsum3_rho, dsum_rho, dsum_rho2, fmax)
        fmaxs(1, j) = fmax
        fmaxmax = fmaxs(1, 1)
        dsum_rho = exp(fmax - fmaxmax) * dsum_rho
        dsum1 = exp(fmax - fmaxmax) * dsum1_rho
        dsum2 = exp(fmax - fmaxmax) * dsum2_rho
        dsum3 = exp(fmax - fmaxmax) * dsum3_rho


        dsum_rho = dsum_rho * wht_phi
        dsumsi(k+1) = dsumsi(k+1) + dsum1 * wht_phi
        dsumsi(k+2) = dsumsi(k+2) + dsum2 * wht_phi
        dsumsi(k+3) = dsumsi(k+3) + dsum3 * wht_phi

        dsum1i = dsum1i + dsum1 * wht_phi
        dsum2i = dsum2i + dsum2 * wht_phi
        dsum3i = dsum3i + dsum3 * wht_phi

        do ijk=1,k
        dsumsi(ijk) = dsumsi(ijk) + vmoms(ijk) * dsum_rho
        enddo        
        dsumi = dsumi + dsum_rho

 120    continue

        call get_xs_from_ws(dsum_xsi,dsumsi,k,k1,k2,vt,t)

        return
        end
c
c
c
c
c
        subroutine theta_int1(t, nn, n, k1, k2, b, y, ynew, ysmall,
     1     resid, dsumi, dsum_xsi, dsum1i, dsum2i, dsum3i)
        implicit real *8 (a-h,o-z)
        real*8 y(*), b(*),ysmall(*),ynew(*),dsum_xsi(*),
     1     rhos(nn+10),phis(nn+10), s(100 000),s2(100 000), 
     1     ys2(100 000), vmoms(100 000),ys(100 000),dsumsi(10000),
     1     whts_thetas(1000),whts_phis(1000),whts_rhos(1000)
        real *8, allocatable :: u(:,:),vt(:,:),asca(:,:),ut(:,:),
     1     fmaxs(:,:),v(:,:)


        done = 1.0d0
        pi = 4*atan(done)

        k = k1+k2

        allocate (u(n,k))
        allocate (ut(k,n))
        allocate (vt(k,k))
        allocate (v(k,k))
        allocate (asca(n,k))
        allocate (fmaxs(nn,nn))


        call rescale_a(t,b,asca,k,k,k1,k2)
        call dsvd(k,k, asca, u, s, vt)
        call y_proj(ynew,u,s,k,k,ysmall,res2)
        call entry_sqr(ysmall,ys2,k)
        call entry_sqr(s,s2,k)

ccc        call prin2('s2 new*', s2, k)
ccc        call prin2('s new*', s, k)
ccc        call prin2('ysmall*', ysmall, k)
ccc        call prin2('ys2 new*', ys2, k)
ccc        call prin2('resid new*', resid, 1)

        dsumi = 0
        dsum1i = 0
        dsum2i = 0
        dsum3i = 0
        do ijk=1,k+3
        dsumsi(ijk)=0
        enddo

c
c        phi integral 
c
        call get_phi1(nn, n, k, t, ysmall, ys2, s, s2, resid, phi1i)
ccc        call prin2('phi1i new*', phi1i, 1)
        phi0 = 0.0
        call lege_nodes_whts(nn, phi0, phi1i, phis, whts_phis)
ccc        call prin2('phis new*', phis, nn)


        do 120 j=1,nn
        phi = phis(j)
        call compute_beta_alpha_prefact(phi,ysmall,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)
        call compute_mjs(k,s2,ysmall,phi,vmoms)

        wht_phi = whts_phis(j)

c
c        rho integral 
c
        call eval_rho_int(nn, n, k, prefact, exp_fact, 
     1     t, phi, rhos, whts_rhos, 
     1     dsum1_rho, dsum2_rho, dsum3_rho, dsum_rho, dsum_rho2, fmax)
        fmaxs(1, j) = fmax
        fmaxmax = fmaxs(1, 1)
ccc        dsum_rho = exp(fmax - fmaxmax) * dsum_rho
ccc        dsum1 = exp(fmax - fmaxmax) * dsum1_rho
ccc        dsum2 = exp(fmax - fmaxmax) * dsum2_rho
ccc        dsum3 = exp(fmax - fmaxmax) * dsum3_rho

        dsum_rho = exp(fmax) * dsum_rho
        dsum1 = exp(fmax) * dsum1_rho
        dsum2 = exp(fmax) * dsum2_rho
        dsum3 = exp(fmax) * dsum3_rho


        dsum_rho = dsum_rho * wht_phi
        dsumsi(k+1) = dsumsi(k+1) + dsum1 * wht_phi
        dsumsi(k+2) = dsumsi(k+2) + dsum2 * wht_phi
        dsumsi(k+3) = dsumsi(k+3) + dsum3 * wht_phi

        dsum1i = dsum1i + dsum1 * wht_phi
        dsum2i = dsum2i + dsum2 * wht_phi
        dsum3i = dsum3i + dsum3 * wht_phi

        do ijk=1,k
        dsumsi(ijk) = dsumsi(ijk) + vmoms(ijk) * dsum_rho
        enddo        
        dsumi = dsumi + dsum_rho


 120    continue

        call get_xs_from_ws(dsum_xsi,dsumsi,k,k1,k2,vt,t)

        return
        end
c
c
c
c
c
        subroutine get_phi1(nn, n, k, t, ysmall, ys2, s, s2, resid,
     1     phi1i)
        implicit real *8 (a-h,o-z)
        real*8 rhos(nn+10),whts_rhos(nn+10),phis0(nn+10),
     1     whts_phis0(nn+10), ysmall(*),ys2(*),s(*),s2(*),
     1     fmaxs(nn+10)


        done = 1.0d0
        pi = 4*atan(done)

        phi0 = 0.0d0
        phi1 = pi/2.0d0
        call lege_nodes_whts(nn, phi0, phi1, phis0, whts_phis0)

        fmax = 0.0d0
        do 119 j=1,nn
ccc        call prin2('fmax*', fmax, 1)
        phi = phis0(j)
        call compute_beta_alpha_prefact(phi,ysmall,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)

c
c        get maximum
c

        dn = n - 2
        rho_max = sqrt(1.0d0/2.0d0*(sqrt(4.0d0*exp_fact + dn*dn) - dn))
        call eval_logdens_rho(rho_max,n,prefact,exp_fact,fmaxmax)
        fmaxs(j) = fmaxmax

c
c        compute integral in rho
c
        call get_rho0(exp_fact, n, rho0)
        call get_rho1(exp_fact, n, rho1)
        call lege_nodes_whts(nn, rho0, rho1, rhos, whts_rhos)
        dsum_rho = 0
        do 110 ik=1,nn
        rho = rhos(ik)
        call eval_logdens_rho(rho,n,prefact,exp_fact,val_sphere)
ccc        call prin2('val_sphere*', val_sphere, 1)
        dsum_rho=dsum_rho+exp(val_sphere-fmaxs(1))*whts_rhos(ik)*rho**2
 110    continue
ccc        call prin2('dsum_rho*', dsum_rho, 1)

        if (dsum_rho .gt. fmax) fmax = dsum_rho
c
c        check if current point is significantly less than the max
c
        thresh = 1.0d-22
        if (dsum_rho/fmax .lt. thresh) then
ccc          call prin2('fmax*', dsum_rho/fmax, 1)
          phi1i = phi
          goto 121
        endif
 119    continue

c
c        must use full interval 
c
        phi1i = pi/2.0d0


 121    continue
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
        real*8 rhos(nn+10),whts_rhos(nn+10),phis0(nn+10),
     1     whts_phis0(nn+10), ysmall(*),ys2(*),s(*),s2(*),
     1     fmaxs(nn+10)


        done = 1.0d0
        pi = 4*atan(done)

        phi0 = 0.0d0
        phi1 = pi/2.0d0
        call lege_nodes_whts(nn, phi0, phi1, phis0, whts_phis0)

        fmax = 0.0d0
        do 119 j=1,nn
        phi = phis0(j)
        call compute_beta_alpha_prefact(phi,ysmall,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)

c
c        get maximum of integrand and use that maximum for the 
c        first value of phi for scaling
c
        dn = n - 2
        rho_max = sqrt(1.0d0/2.0d0*(sqrt(4.0d0*exp_fact + dn*dn) - dn))
        call eval_logdens_rho(rho_max,n,prefact,exp_fact,fmaxmax)
        fmaxs(j) = fmaxmax

c
c        compute integral in rho
c
        call eval_rho_int2(nn, n, k, exp_fact, 
     1     dsum_rho, dsum_rho1, dsum_rho2, fmaxi)
        fi = dsum_rho * exp(prefact + fmaxi - fmaxs(1))
        fi = log(dsum_rho) + prefact + fmaxi - fmaxs(1)
ccc        call prin2('dsum_rho2*', dsum_rho, 1)

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
c        so return that value, readjusting for the scaling
c
        fmax = fmax + fmaxs(1)

        return
        end
c
c
c
c
c
        subroutine eval_rho_int2(nn, n, k, exp_fact, 
     1     dsum_rho, dsum_rho1, dsum_rho2, fmax)
        implicit real *8 (a-h,o-z)
        real*8 whts_rhos(nn),rhos(nn)

        call get_rho0(exp_fact, n, rho0)
        call get_rho1(exp_fact, n, rho1)
        call lege_nodes_whts(nn, rho0, rho1, rhos, whts_rhos)

c
c        get maximum
c
        dn = n - 2
        rho_max = sqrt(1.0d0/2.0d0*(sqrt(4.0d0*exp_fact + dn*dn) - dn))
        call eval_logdens_rho2(rho_max,n,exp_fact,fmax)
ccc        fmax = fmax + prefact
        
c
c        compute integral with respect to rho
c
        dsum_rho = 0
        dsum_rho1 = 0
        dsum_rho2 = 0
        do 110 i=1,nn
        rho = rhos(i)
        call eval_logdens_rho2(rho,n,exp_fact,f)
        wht_rho = whts_rhos(i)

        dsum_rho = dsum_rho + exp(f-fmax)*wht_rho
        dsum_rho1 = dsum_rho1 + rho*exp(f-fmax)*wht_rho
        dsum_rho2 = dsum_rho2 + rho**2*exp(f-fmax)*wht_rho

c        call prin2('dsum_rho*', dsum_rho, 1)
c        call prin2('f*', f, 1)
c        call prin2('fmax*', fmax, 1)

 110    continue


        return
        end
c
c
c
c
c
        subroutine eval_rho_int(nn, n, k, prefact, exp_fact, 
     1     t, phi, rhos, whts_rhos,  
     1     dsum1, dsum2, dsum3, dsum_rho, dsum_rho2, fmax)
        implicit real *8 (a-h,o-z)
        real*8 whts_rhos(*),rhos(*)

        call get_rho0(exp_fact, n, rho0)
        call get_rho1(exp_fact, n, rho1)
        call lege_nodes_whts(nn, rho0, rho1, rhos, whts_rhos)

        dsum1 = 0.0d0
        dsum2 = 0.0d0
        dsum3 = 0.0d0
c
c        get maximum
c
        dn = n - 2
        rho_max = sqrt(1.0d0/2.0d0*(sqrt(4.0d0*exp_fact + dn*dn) - dn))
        call eval_logdens_rho2(rho_max,n,prefact,exp_fact,fmax)
        
c
c        compute integral with respect to rho
c
        dsum_rho = 0
        dsum_rho2 = 0
        do 110 ik=1,nn
        rho = rhos(ik)
        call eval_logdens_rho2(rho,n,prefact,exp_fact,f)
        wht = whts_rhos(ik)*sin(phi)

        dsum_rho = dsum_rho + exp(f-fmax)*wht
        dsum_rho2 = dsum_rho2 + rho**2*exp(f-fmax)*wht

        sig1 = rho*cos(phi)
        sig2 = rho*sin(phi)*cos(t)
        sig3 = rho*sin(phi)*sin(t)

        dsum1 = dsum1 + sig1*exp(f-fmax)*wht
        dsum2 = dsum2 + sig2*exp(f-fmax)*wht
        dsum3 = dsum3 + sig3*exp(f-fmax)*wht
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

ccc        call prin2('phi = *',phi,1)
ccc        call prinf('k = *',k,1)
ccc        call prin2('ys = *',ys,k)

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
        subroutine eval_logdens_rho2(rho,n,exp_fact,val)
        implicit real *8 (a-h,o-z)


ccc        val = prefact - exp_fact/2.0d0/(rho*rho)-rho*rho/2.0d0
        val = -exp_fact/2.0d0/(rho*rho)
c        call prin2('rho*', rho, 1)
c        call prin2('val*', val, 1)
        val = val -rho*rho/2.0d0
c        call prin2('val*', val, 1)
c        call prin2('exp_fact*', exp_fact, 1)

        val = val - (n-2)*log(rho)
c        call prin2('val*', val, 1)

        return
        end
c
c
c
c
c
        subroutine eval_logdens_rho(rho,n,prefact,exp_fact,val)
        implicit real *8 (a-h,o-z)

        val = prefact - exp_fact/2.0d0/(rho*rho)-rho*rho/2.0d0
        val = val - n*log(rho)

        return
        end
c
c
c
c
c
        subroutine compute_beta_alpha_prefact2(phi,ys,ys2,s,s2,
     1      n,k,resid,alpha,beta,prefact,exp_fact)
        implicit real *8 (a-h,o-z)
        dimension ys(*),ys2(*),s(*),s2(*)

        done = 1
        pi   = atan(done)*4

        sp = sin(phi)**2
        cp = cos(phi)**2

        alpha = 0

        do i=1,k
        alpha = alpha - log(cp+sp*s2(i))
        enddo

        alpha = alpha/2 + k/2.0d0*log(2*pi)

        beta = 0

        do i=1,k
        beta = beta + ys2(i)/(s2(i)*sp+cp)
        enddo

        prefact = -(n-k)/2.0d0*log(cp)+alpha + log(sin(phi))
        exp_fact= resid/cp+beta


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

        done = 1
        pi   = atan(done)*4

        sp = sin(phi)**2
        cp = cos(phi)**2

        alpha = 0

        do i=1,k
        alpha = alpha - log(cp+sp*s2(i))
        enddo

        alpha = alpha/2 + k/2.0d0*log(2*pi)

        beta = 0

        do i=1,k
        beta = beta + ys2(i)/(s2(i)*sp+cp)
        enddo

        prefact = -(n-k)/2.0d0*log(cp)+alpha
        exp_fact= resid/cp+beta


        return
        end
c
c
c
c
c
        subroutine get_xs_from_ws(xs,ws,k,k1,k2,vt,t)
        implicit real *8 (a-h,o-z)
        dimension xs(*),ws(*),vt(k,k)

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
        subroutine eval_logdens_with_m(sig,r,n,k,ys,s,
     1      ys2,s2,resid,val,ws)
        implicit real *8 (a-h,o-z)
        dimension ws(*),ys(*),s(*),ys2(*),s2(*)

        call eval_logdens(sig,r,n,k,ys2,s2,resid,val)

        sig2 = sig*sig
        r2   = r*r

        do i=1,k
        ws(i) = (ys(i)*s(i))/(s2(i)+sig2/r2)
        enddo

        return
        end
c
c
c
c
c
        subroutine eval_logdens(sig,r,n,k,ys2,s2,resid,val)
        implicit real *8 (a-h,o-z)
        dimension ys2(*),s2(*)

        sig2 = sig*sig
        r2   = r*r

        call prin2('resid = *',resid,1)

        call eval_prefac(sig2,r2,n,k,vprefac)
        call eval_expons(sig2,r2,s2,k,ys2,resid,vexp)

        call prin2('vprefac = *',vprefac,1)
        call prin2('vexp = *',vexp,1)
        val = vprefac + vexp 

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
        subroutine eval_prefac(sig2,r2,n,k,val)
        implicit real *8 (a-h,o-z)

        val = 0
        val = -sig2 - r2
        val = val - (n-k)*log(sig2)-k*log(r2)
        val = val/2.0d0

        return
        end
c
c
c
c
c
        subroutine eval_expons(sig2,r2,s2,k,ys2,resid,val)
        implicit real *8 (a-h,o-z)
        dimension s2(*),ys2(*)

        val = 0
        vale =0

        do i=1,k
ccc        fac = -1.0d0/(s2(i))/(r2*s2(i)+sig2)
        fac = -1.0d0/(r2*s2(i)+sig2)
        vale = vale + (ys2(i))*fac
        call prin2('fac = *',s2(i)+sig2/r2,1)
        enddo
        call prin2('ys2 = *',ys2,k)

        vale = vale/2.0d0 - resid/2.0d0/sig2
        call prin2('resid = *',resid,1)
        call prin2('resid val = *',resid/2/sig2,1)

        val = 0
        do i=1,k
        val = val + log(s2(i)+sig2/r2)
        enddo

        val = -val/2.0d0 + k*log(atan(1.0d0)*8.0d0)/2.0d0
        call prin2('atan = *',atan(1.0d0)*8.0d0,1)
        val = val + vale

        return
        end
c
c
c
c
c
        subroutine y_proj_with_d(y,us,ss,n,k,ysmall,resid,
     1      yd)
        implicit real *8 (a-h,o-z)
        dimension y(n),us(n,k),ss(*),ysmall(k),yd(n)

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

ccc        call prin2('yy = *',yy,1)
ccc        call prin2('rr = *',rr,1)

        do i=1,n
        yt = y(i)
        do j=1,k
        yt = yt - us(i,j)*ysmall(j)
        enddo
        yd(i) = yt
        enddo
        

        resid = yy-rr

cccccc
cccccc         .   .   .   to delete
ccc        
ccc        dinprod = 0
ccc        dinprod2= 0
ccc        do i=1,n
ccc        dinprod = dinprod + us(i,1)*us(i,1)
ccc        dinprod2= dinprod2+ us(i,1)*us(i,2)
ccc        enddo
ccc        call prin2('dinprod = *',dinprod,1)
ccc        call prin2('dinprod2= *',dinprod2,1)

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

ccc        call prin2('yy = *',yy,1)
ccc        call prin2('rr = *',rr,1)

        resid = yy-rr

cccccc
cccccc         .   .   .   to delete
ccc        
ccc        dinprod = 0
ccc        dinprod2= 0
ccc        do i=1,n
ccc        dinprod = dinprod + us(i,1)*us(i,1)
ccc        dinprod2= dinprod2+ us(i,1)*us(i,2)
ccc        enddo
ccc        call prin2('dinprod = *',dinprod,1)
ccc        call prin2('dinprod2= *',dinprod2,1)

        return
        end
c
c
c
c
c
        subroutine diagonalize_subs(a,n,k,k1,k2,u1,s1,vt1,
     1      u2,s2,vt2,a1,a2)
        implicit real *8 (a-h,o-z)
        dimension a(n,k),a1(n,k1),a2(n,k2),u1(*),s1(*),vt1(*),
     1      u2(*),s2(*),vt2(*)

        do i=1,n
        do j=1,k1
        a1(i,j) = a(i,j)
        enddo
        enddo

        do i=1,n
        do j=1,k2
        a2(i,j) = a(i,j+k1)
        enddo
        enddo

        call dsvd(n,k1, a1, u1, s1, vt1)
        call dsvd(n,k2, a2, u2, s2, vt2)

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

        dn = n - 2
        rho_max = sqrt(1/2.0*(sqrt(4*ci + dn**2) - dn))
        call eval_log_psi_rho(rho_max, ci, n, f_max)
        f_max_log10 = f_max / log(10.0d0)

        dder2 = -1-3*ci/rho_max**4 + (n-2)/rho_max**2
        sd = 1/sqrt(-dder2)

c
c        find left endpoint
c        
        xmin = 1.0d-16
        xmax = rho_max 
        do 110 i=1,10
        rho0 = (xmin + xmax)/2.0d0
        call eval_log_psi_rho(rho0, ci, n, f)
        f_log10 = f / log(10.0d0)
        if (f_log10 .lt. f_max_log10-18) xmin = rho0
        if (f_log10 .ge. f_max_log10-18) xmax = rho0
 110    continue

ccc        call eval_log_psi_rho(rho0, ci, n, f)
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

        dn = n - 2
        rho_max = sqrt(1/2.0*(sqrt(4*ci + dn**2) - dn))
        call eval_log_psi_rho(rho_max, ci, n, f_max)
        f_max_log10 = f_max / log(10.0d0)

        dder2 = -1-3*ci/rho_max**4 + (n-2)/rho_max**2
        sd = 1/sqrt(-dder2)

c
c        find right endpoint
c        
        xmin = rho_max
        xmax = rho_max + 40*sd
        do 120 i=1,10
        rho1 = (xmin + xmax)/2.0d0
        call eval_log_psi_rho(rho1, ci, n, f)

        f_log10 = f / log(10.0d0)
        if (f_log10 .lt. f_max_log10-18) xmax = rho1
        if (f_log10 .ge. f_max_log10-18) xmin = rho1
 120    continue

ccc        call prin2('rho1*', rho1, 1)
ccc        call eval_log_psi_rho(rho1, ci, n, f)
ccc        call prin2('f at right endpoint*', exp(f), 1)
        if (f .gt. -20) then
ccc          call prin2('bad!*', 1.0d0, 1)
ccc          call prin2('f*', f, 1)
        endif

        return
        end
c
c
c
c
c
        subroutine eval_log_psi_rho(rho, c, n, f)
        implicit real *8 (a-h,o-z)

        f = -rho**2/2.0d0 - c/(2*rho**2)
        f = f - (n-2)*log(rho)

        return
        end
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
c
c
c
        subroutine gaus3d(nn, n, k1, k2, a, y, dints, dsum)
        implicit real *8 (a-h,o-z)
        real*8 a(1 000 000), x(100 000), y(100 000), b(100 000), 
     1     at(100 000), c_inv(100 000), d2(100 000), c(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     xexps(10 000), xexps2(10 000), ts(10 000), dexps_x(10 000),
     1     dsums(10 000), dints(10 000), fs(nn, nn), u(1), v(1),
     1     whts(10 000), dints_sig(1), yta(10 000), ata(10 000)

        k = k1+k2
c
c        find expectations dumb
c
c
c        legendre nodes
c        
        t0 = 0.0
        t1 = 15.0
        call lege_nodes_whts(nn, t0, t1, ts, whts)
ccc        call prin2('ts*', ts, nn)
ccc        call prin2('whts*', whts, nn)

c
c        compute ata 
c
        call mat_trans(a, n, k, at)
        call mat_mult(at, k, n, a, k, ata)
        call mat_mult(y, 1, n, a, k, yta)

c
c        initialize sums to 0
c

        dsum = 0
        do 80 i=1,k+3
        dsums(i) = 0
        dexps_x(i) = 0
 80     continue

c
c        take sum
c
        do 130 i=1,nn
        do 120 j=1,nn
        do 110 ik=1,nn
        sig1 = ts(i)
        sig2 = ts(j)
        sig3 = ts(ik)

        call q_tilde(n, k1, k2, a, y, sig1, sig2, sig3, f, dexps_x)
ccc        call prin2('f*', f, 1)
        f = exp(f)

        wht = whts(i)*whts(j)*whts(ik)
        dsum = dsum + f * wht
        fs(i, j) = f

        do 90 ijk=1,k
        term1 = dexps_x(ijk) * f * wht
        dsums(ijk) = dsums(ijk) + term1
 90     continue
        dsums(k+1) = dsums(k+1) + sig1 * f * wht
        dsums(k+2) = dsums(k+2) + sig2 * f * wht
        dsums(k+3) = dsums(k+3) + sig3 * f * wht


ccc        call corner_ind_3d(nn, i, j, ik, coef)
ccc        if ((i .eq. nn) .or. ((j .eq. nn) .or. (ik .eq. nn))) then
ccc          call prin2('f*', f, 1)
ccc        endif

 110    continue
 120    continue
 130    continue

c
c        x expectations
c
        do 140 i=1,k+3
        dints(i) = dsums(i) / dsum
 140    continue

ccc        call prin2('dints*', dints, k+3)
ccc        call prin2('dconst*', dconst, 1)

        return
        end
c
c
c
c
c
        subroutine q_tilde2(n, k1, k2, a, y, sig1, sig2, sig3,
     1     f_const, dexps_x)
        implicit real *8 (a-h,o-z)
        real*8 a(1), x(1), y(1),
     1     at(100 000), d2(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     ax(10 000), s1(10 000), 
     1     b(10 000), b2(10 000), ata(10 000), ux(10 000),
     1     au(10 000), yta(10 000), u(10 000), ut(10 000),
     1     s2(10 000), uxt(10 000), z(10 000), w(10 000),
     1     dexps_z(10 000)


        k = k1+k2
        done = 1
        pi = atan(done)*4

        call get_s_s2(k1, k2, sig2, sig3, s, s1)
c
c        term1
c
        call inner_prod(y, y, n, yty)
        a0 = yty/(2*sig1**2) 
c
c        construct b
c

        call mat_trans(a, n, k, at)
        call mat_mult(at, k, n, a, k, ata)
c
c        scale ata 
c

        do 90 i=1,k*k
ccc        b(i) = ata(i)/(2*sig1**2)
        b(i) = ata(i)
 90     continue

c
c        add diagonal 
c
        do 100 i=1,k
        ind = i + (i-1)*k
ccc        b(ind) = b(ind) + s1(i)
        b(ind) = b(ind) + s1(i)*2*sig1**2
 100    continue

ccc        call prin2('ata*', ata, k*k)
ccc        call prin2('b*', b, k*k)

c
c        eigendecomposition of b
c
        ifvect=1
        eps=1.0d-14
        call cc_arr(k*k, b, s2)
        call jaceig(s2,k,eps,ifvect,u,nelems)
        call prinf('nelems= *',nelems,1)

        call mat_trans(u, k, k, ut)
        call mat_diag(s2, k, d2)
ccc        call prin2('s2*', s2, k)
c
c        construct z 
c

        call mat_vec_mult(ut, k, k, x, z)

c
c        evaluate xt * b * x
c
        call mat_mult(a, n, k, u, k, au)
        call mat_mult(y, 1, n, au, k, w)

        f_const = 0
        do 200 i=1,k
ccc        a1 = w(i)/sig1**2
        a1 = 2*w(i)
        a2 = s2(i)

        dexps_z(i) = a1/(2*a2)
ccc        f_const = f_const + a1**2/(4*a2) + log(sqrt(2*pi/(2*a2)))
        f_const = f_const - yty/(k*2*sig1**2)
        f_const = f_const + a1**2/(4*a2*2*sig1**2)
        f_const = f_const + log(sqrt(2*pi/(2*(a2/(2*(sig1**2))))))
 200    continue

ccc        f_const = f_const - a0
        f_const = f_const - k1*log(sig2) - k2*log(sig3) - n*log(sig1)
        f_const = f_const - sig1**2/2.0 - sig2**2/2.0 - sig3**2/2.0 
ccc        call prin2('f_const*', f_const, 1)

c
c        take sum of terms
c

        call mat_mult(u, k, k, dexps_z, k, dexps_x)

        return
        end
c
c
c
c
c
        subroutine q_tilde2_old(n, k1, k2, a, y, sig1, sig2, sig3,
     1     f_const, dexps_x)
        implicit real *8 (a-h,o-z)
        real*8 a(1), x(1), y(1),
     1     at(100 000), d2(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     ax(10 000), s1(10 000), 
     1     b(10 000), b2(10 000), ata(10 000), ux(10 000),
     1     au(10 000), yta(10 000), u(10 000), ut(10 000),
     1     s2(10 000), uxt(10 000), z(10 000), w(10 000),
     1     dexps_z(10 000)


        k = k1+k2
        done = 1
        pi = atan(done)*4

        call get_s_s2(k1, k2, sig2, sig3, s, s1)
c
c        term1
c
        call inner_prod(y, y, n, yty)
ccc        a0 = yty/(2*sig1**2) 
        a0 = yty
c
c        construct b
c

        call mat_trans(a, n, k, at)
        call mat_mult(at, k, n, a, k, ata)
c
c        scale ata 
c

        do 90 i=1,k*k
ccc        b(i) = ata(i)/(2*sig1**2)
        b(i) = ata(i)
 90     continue

c
c        add diagonal 
c
        do 100 i=1,k
        ind = i + (i-1)*k
        b(ind) = b(ind) + 2*sig1**2*s1(i)
 100    continue


ccc        call prin2('ata*', ata, k*k)
ccc        call prin2('b*', b, k*k)

c
c        eigendecomposition of b
c
        ifvect=1
        eps=1.0d-14
        call cc_arr(k*k, b, s2)
        call jaceig(s2,k,eps,ifvect,u,nelems)
        call mat_trans(u, k, k, ut)
        call mat_diag(s2, k, d2)
ccc        call prin2('s2*', s2, k)
c
c        construct z 
c

        call mat_vec_mult(ut, k, k, x, z)

c
c        evaluate xt * b * x
c
        call mat_mult(a, n, k, u, k, au)
        call mat_mult(y, 1, n, au, k, w)

        f_const = 0
        do 200 i=1,k
ccc        a1 = w(i)/sig1**2
        a1 = 2*w(i)
        a2 = s2(i)

        dexps_z(i) = a1/(2*a2)
ccc        f_const = f_const + a1**2/(4*a2) + log(sqrt(2*pi/(2*a2)))
        
        f_const = f_const + a1**2/(4*a2)
        f_const = f_const + log(sqrt(2*pi/(2*(a2/(2*(sig1**2))))))
 200    continue

        f_const = f_const - a0
        f_const = f_const - k1*log(sig2) - k2*log(sig3) - n*log(sig1)
        f_const = f_const - sig1**2/2.0 - sig2**2/2.0 - sig3**2/2.0 

c
c        take sum of terms
c

        call mat_mult(u, k, k, dexps_z, k, dexps_x)

        return
        end
c
c
c
c
c
        subroutine test_coefs()
        implicit real *8 (a-h,o-z)
        real*8 a(1 000 000), x(100 000), y(100 000), b(100 000), 
     1     at(100 000), c_inv(100 000), d2(100 000), c(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     xexps(10 000), xexps2(10 000), ts(10 000), dexps_x(10 000),
     1     dsums(10 000), dints(10 000), dds(10 000), dints2(10 000),
     1     fs(10 000), u(10 000), v(10 000), whts(10 000),
     1     coefs(10 000), xs(10 000)
c
c        parameters
c
        n = 3
        k1 = 2
        k2 = 1
        k = k1+k2

        a(1) = 1.0
        a(2) = 1.0
        a(3) = 1.0
        a(4) = 1.0
        a(5) = 1.0
        a(6) = 0.0
        a(7) = 1.0
        a(8) = 0.0
        a(9) = 0.0

        a(1) = 1.0
        a(2) = 1-1.0d-5
        a(3) = 1-1.0d-5
        a(4) = 1-1.0d-5
        a(5) = 1.0
        a(6) = 1-1.0d-5
        a(7) = 1-1.0d-5
        a(8) = 1-1.0d-5
        a(9) = 1.0

        y(1) = 1.0
        y(2) = 0.1
        y(3) = 0.5


        sig1 = 0.1
        sig2 = 0.1
c
c        construct legendre nodes
c
        nn = 100
        t0 = 1e-10
ccc        t0 = 1e-1
        t1 = 5.0
        itype=2
        call legeexps(itype,nn,x,u,v,whts)
c
c        shift to t0 to t1
c
        do 60 i=1,nn
        tmp = (x(i) + 1)/2.0d0
        x(i) = t0 + (t1-t0)*tmp
 60     continue
ccc        call prin2('x*', x, nn)

        do 100 i=1,nn
        sig3 = x(i)
        call q_tilde(n, k1, k2, a, y, sig1, sig2, sig3, f, dexps_x)
        fs(i) = exp(f) * dexps_x(1)
 100    continue

        call prin2('fs*', fs, nn)

        call mat_vec_mult(u, nn, nn, fs, coefs)

        call prin2('coefs*', coefs, nn)
        do 110 i=1,nn
        xs(i) = i
        coefs(i) = log(abs(coefs(i)))/log(10.0d0)
 110    continue
        call prin2('coefs*', coefs, nn)

        do 120 i=1,nn-1
        tmp = coefs(i)/coefs(i+1)
ccc        call prin2('dd*', tmp, 1)
 120    continue

        iw = 12
        itype1 = 3
        call quagraph(iw,x,fs,nn,itype1,'title*')

        iw = 21
        itype1 = 3
        call quagraph(iw,xs,coefs,nn,itype1,'title*')


        return
        end
c
c
c
c
c
        subroutine test4_2d()
        implicit real *8 (a-h,o-z)
        real*8 a(1 000 000), x(100 000), y(100 000), b(100 000), 
     1     at(100 000), c_inv(100 000), d2(100 000), c(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     xexps(10 000), xexps2(10 000), ts(10 000), dexps_x(10 000),
     1     dsums(10 000), dints(10 000), dds(10 000), dints2(10 000)
c
c        parameters
c
        n = 3
        k1 = 2
        k2 = 1
        k = k1+k2
        call corrand_norm(n*k,a,y)
        call corrand_norm(k,x,tmp_mat)

        a(1) = 1.0
        a(2) = 1.0
        a(3) = 1.0
        a(4) = 1.0
        a(5) = 1.0
        a(6) = 0.0
        a(7) = 1.0
        a(8) = 0.0
        a(9) = 0.0

        y(1) = 1.0
        y(2) = 0.1
        y(3) = 0.5

        nn = 80
        call trap_2d(nn, n, k1, k2, a, x, y, dints)
ccc        call gaus_2d(nn, n, k1, k2, a, x, y, dints)
        call prin2('dints*', dints, k)

        nn2 = 2*nn
        call trap_2d(nn2, n, k1, k2, a, x, y, dints2)
ccc        call gaus_2d(nn2, n, k1, k2, a, x, y, dints2)
        call prin2('dints2*', dints2, k)

        do 100 i=1,k
        dds(i) = (dints(i) - dints2(i))/(dints(i)+dints2(i))*2
 100    continue

        call prin2('dds*', dds, k)

        return
        end
c
c
c
c
c
        subroutine gaus_2d(nn, n, k1, k2, a, x, y, dints)
        implicit real *8 (a-h,o-z)
        real*8 a(1 000 000), x(100 000), y(100 000), b(100 000), 
     1     at(100 000), c_inv(100 000), d2(100 000), c(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     xexps(10 000), xexps2(10 000), ts(10 000), dexps_x(10 000),
     1     dsums(10 000), dints(10 000), fs(nn, nn), u(1), v(1),
     1     whts(10 000)

        k = k1+k2
c
c        find expectations dumb
c
        t0 = 0.0
        t1 = 5.0
        call equispaced_nodes(nn, t0, t1, ts)

c
c        legendre nodes
c        
        itype=1
        call legeexps(itype,nn,x,u,v,whts)
ccc        call prin2('x*', x, nn)
ccc        call prin2('whts*', whts, nn)

c
c        shift to t0 to t1
c
        do 60 i=1,nn
        tmp = (x(i) + 1)/2.0d0
        x(i) = t0 + (t1-t0)*tmp
 60     continue
ccc        call prin2('x*', x, nn)

c
c        initialize sums to 0
c

        dsum = 0
        do 80 i=1,k
        dsums(i) = 0
        dexps_x(i) = 0
 80     continue

c
c        take sum
c
ccc        sig3 = 0.01
        sig2 = 0.01
        sig1 = 0.5
        do 130 i=1,nn
        do 120 j=1,nn
ccc        sig1 = x(i)
        sig2 = x(j)
        sig3 = x(j)

ccc        call q2_log(n, k1, k2, a, y, sig1, sig2, sig3, f, dexps_x)
ccc        call prin2('f*', f, 1)

        call q_tilde(n, k1, k2, a, y, sig1, sig2, sig3, f, dexps_x)
ccc        call prin2('f2*', f, 1)
        f = exp(f)

        dsum = dsum + f * whts(i)*whts(j)
        fs(i, j) = f

        do 90 ijk=1,k
        dsums(ijk) = dsums(ijk) + dexps_x(ijk)*f*whts(i)*whts(j)
 90     continue

 120    continue
 130    continue

c
c        integrate
c
        dsum = dsum*((t1 - t0) / 2.0)**2
ccc        call prin2('dsum*', dsum, 1)

        do 140 i=1,k
        dints(i) = ((t1 - t0) / 2.0)**2 * dsums(i) / dsum
ccc        dints(i) = dsums(i)
 140    continue

ccc        call prin2('dconst*', dconst, 1)
ccc        call prin2('dsum1*', dsum1, 1)

ccc        call prin2('dints*', dints, k)

c
c        plot
c
        tol = 40
        iw = 11
        call plot_heatmap(iw, nn, tol, fs)

        return
        end
c
c
c
c
c
        subroutine trap_2d(nn, n, k1, k2, a, x, y, dints)
        implicit real *8 (a-h,o-z)
        real*8 a(1 000 000), x(100 000), y(100 000), b(100 000), 
     1     at(100 000), c_inv(100 000), d2(100 000), c(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     xexps(10 000), xexps2(10 000), ts(10 000), dexps_x(10 000),
     1     dsums(10 000), dints(10 000), fs(nn, nn), u(1), v(1),
     1     whts(10 000)

        k = k1+k2
c
c        find expectations dumb
c
        t0 = 1e-8
        t0 = 1e-6
        t1 = 5.0
        call equispaced_nodes(nn, t0, t1, ts)

c
c        take sum 
c
        dsum = 0
        do 80 i=1,k
        dsums(i) = 0
        dexps_x(i) = 0
 80     continue


        sig1 = 0.1
        do 130 i=1,nn
        do 120 j=1,nn
        sig2 = ts(i)
        sig3 = ts(j)

        call q2_log(n, k1, k2, a, y, sig1, sig2, sig3, f, dexps_x)
        f = exp(f)
        fs(i, j) = f

        call corner_ind_2d(nn, i, j, coef)
        dsum = dsum + f*coef
        do 96 ijk=1,k
        dsums(ijk) = dsums(ijk) + dexps_x(ijk)*f*coef
 96     continue

ccc        call prin2('dexps_x(2)*', dexps_x(2), 1)
ccc        call prin2('exp(f)*', exp(f), 1)
ccc        fs(i, j) = dexps_x(2)*exp(f)
ccc        call prin2('fs*', log(abs(fs(i, j))), 1)

 120    continue
 130    continue

ccc        call prin2('fs*', fs, nn**2)
ccc        call max_vec(nn*nn, fs, f1, ind)
ccc        call prin2('max fs*', f1, 1)

c
c        integrate
c
ccc        call prin2('dsum*', dsum, 1)

        do 140 i=1,k
        dints(i) = dsums(i) / dsum
 140    continue
ccc        call prin2('dints*', dints, k)
c
c        plot
c
        tol = 40
        iw = 11
        call plot_heatmap(iw, nn, tol, fs)

        return
        end
c
c
c
c
c
        subroutine test3()
        implicit real *8 (a-h,o-z)
        real*8 a(1 000 000), x(100 000), y(100 000), b(100 000), 
     1     at(100 000), c_inv(100 000), d2(100 000), c(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     xexps(10 000), xexps2(10 000), ts(10 000), dexps_x(10 000),
     1     dsums(10 000), dints(10 000), dds(10 000)
c
c        parameters
c
        n = 3
        k1 = 2
        k2 = 1
        k = k1+k2
        sig1 = 0.5
        sig2 = 2.0
        sig3 = 4.0

        call corrand_norm(n*k,a,y)
        call corrand_norm(k,x,tmp_mat)
        a(1) = 1.0
        a(2) = 1.0
        a(3) = 1.0
        a(4) = 1.0
        a(5) = 1.0
        a(6) = 0.0
        a(7) = 1.0
        a(8) = 0.0
        a(9) = 0.0

        y(1) = 1.0
        y(2) = 0.1
        y(3) = 0.5

c
c        find expectations dumb
c

        nn = 40
        t0 = -5.0
        t1 = 6.0
        call equispaced_nodes(nn, t0, t1, ts)

        dsum = 0
        do 80 i=1,k
        dsums(i) = 0
 80     continue
        do 130 i=1,nn
        do 120 j=1,nn
        do 110 ik=1,nn
        x(1) = ts(i)
        x(2) = ts(j)
        x(3) = ts(ik)
        call q2_log_old(n, k1, k2, a, x, y, sig1, sig2, sig3, 
     1     f, f_const, dexps_x)
        call q_tilde(n, k1, k2, a, x, y, sig1, sig2, sig3, tot,
     1     f_log_const, dexps_x)

        f = exp(f)
        dsum = dsum + f

        do 90 ijk=1,k
        dsums(ijk) = dsums(ijk) + x(ijk)*f
 90     continue

 110    continue
 120    continue
 130    continue
        

        box_vol = (ts(2) - ts(1))**3

        dconst = dsum * box_vol
        call prin2('dconst*', dconst, 1)

c
c        formulas for constant and expectations
c
        f_const = exp(f_log_const)
        call prin2('f_const*', f_const, 1)

        dd = (f_const - dconst) / f_const
        call prin2('dd*', dd, 1)

        do 140 i=1,k
        dints(i) = dsums(i) * box_vol / dconst
 140    continue

ccc        call prin2('dsum1*', dsum1, 1)

        call prin2('dints*', dints, k)
        call prin2('dexps_x*', dexps_x, k)

        do 160 i=1,k
        dds(i) = (dints(i) - dexps_x(i)) / dints(i)
 160    continue

        call prin2('dds*', dds, k)

        return
        end
c
c
c
c
c
        subroutine q_tilde(n, k1, k2, a, y, sig1, sig2, sig3,
     1     f_const, dexps_x)
        implicit real *8 (a-h,o-z)
        real*8 a(1), x(1), y(1),
     1     at(100 000), d2(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     ax(10 000), s1(10 000), 
     1     b(10 000), b2(10 000), ata(10 000), ux(10 000),
     1     au(10 000), yta(10 000), u(10 000), ut(10 000),
     1     s2(10 000), uxt(10 000), z(10 000), w(10 000),
     1     dexps_z(10 000)


        k = k1+k2
        done = 1
        pi = atan(done)*4

        call get_s_s2(k1, k2, sig2, sig3, s, s1)
ccc        call prin2('s1 = *',s1,3)
c
c        term1
c
        call inner_prod(y, y, n, yty)
        a0 = yty/(2*sig1**2) 
c
c        construct b
c

        call mat_trans(a, n, k, at)
        call mat_mult(at, k, n, a, k, ata)
c
c        scale ata 
c

        do 90 i=1,k*k
        b(i) = ata(i)/(2*sig1**2)
 90     continue

c
c        add diagonal 
c
        do 100 i=1,k
        ind = i + (i-1)*k
        b(ind) = b(ind) + s1(i)
 100    continue

ccc        call prin2('ata*', ata, k*k)
ccc        call prin2('b*', b, k*k)

c
c        eigendecomposition of b
c
        ifvect=1
        eps=1.0d-14
        call cc_arr(k*k, b, s2)
ccc        call prin2('b = *',b,k*k)
ccc        call prin2('s2 = *',s2,k*k)
        call jaceig(s2,k,eps,ifvect,u,nelems)
ccc        call prin2('s2 = *',s2,k*k)
ccc        call prinf('nelems=  *',nelems,1)
        call mat_trans(u, k, k, ut)
        call mat_diag(s2, k, d2)
ccc        call prin2('s2*', s2, k)
ccc        call prin2('here *',0,0)
c
c        construct z 
c

        call mat_vec_mult(ut, k, k, x, z)

c
c        evaluate xt * b * x
c
        call mat_mult(a, n, k, u, k, au)
        call mat_mult(y, 1, n, au, k, w)

        f_const = 0
        do 200 i=1,k
        a1 = w(i)/sig1**2
        a2 = s2(i)

        dexps_z(i) = a1/(2*a2)
        f_const = f_const + a1**2/(4*a2) + log(sqrt(2*pi/(2*a2)))
 200    continue

        f_const = f_const - a0
        f_const = f_const - k1*log(sig2) - k2*log(sig3) - n*log(sig1)
        f_const = f_const - sig1**2/2.0 - sig2**2/2.0 - sig3**2/2.0 
ccc        call prin2('f_const*', f_const, 1)

c
c        take sum of terms
c

        call mat_mult(u, k, k, dexps_z, k, dexps_x)

        return
        end
c
c
c
c
c
        subroutine q2_log_old(n, k1, k2, a, x, y, sig1, sig2, sig3, f,
     1     f_const, dexps_x)
        implicit real *8 (a-h,o-z)
        real*8 a(1), x(1), y(1), c(10 000), c_inv(10 000),
     1     w(10 000), d1(10 000), z(10 000), dexps_x(1), dexps_z(1000)

        done = 1
        pi = atan(done)*4
ccc        call prin2('pi*', pi, 1)

c
c        precomp
c
        k = k1+k2
        call precomp(n, k1, k2, sig2, sig3, a, y, c, c_inv, w, d1)
        call inner_prod(y, y, n, yty)
        call mat_vec_mult(c_inv, k, k, x, z)
        a0 = yty/(2*sig1**2) 

c
c        compute term2
c
        dsum = 0
        f_const = 0
        do 200 i=1,k
        a1 = w(i)/sig1**2
        a2 = 1/(2*sig1**2) + d1(i)

        dsum = dsum + a2*(z(i)-a1/(2*a2))**2 - a1**2/(4*a2)
        f_const = f_const + a1**2/(4*a2) + log(sqrt(2*pi/(2*a2)))
        dexps_z(i) = a1/(2*a2)
 200    continue

c
c        magnitude of pdf, depends on x
c
        f = -a0 - dsum
        f = f - k1*log(sig2) - k2*log(sig3) - n*log(sig1)
        f = f - sig1**2/2.0 - sig2**2/2.0 - sig3**2/2.0 

c
c        normalizing constant
c
        f_const = f_const - a0 
        f_const = f_const - k1*log(sig2) - k2*log(sig3) - n*log(sig1)
        f_const = f_const - sig1**2/2.0 - sig2**2/2.0 - sig3**2/2.0 

c
c        back to x expectations
c
        call mat_mult(c, k, k, dexps_z, k, dexps_x)

        return
        end
c
c
c
c
c
        subroutine q2_log(n, k1, k2, a, y, sig1, sig2, sig3,
     1     f_const, dexps_x)
        implicit real *8 (a-h,o-z)
        real*8 a(1), x(1), y(1), c(10 000), c_inv(10 000),
     1     w(10 000), d1(10 000), z(10 000), dexps_x(1), dexps_z(1000)

        done = 1
        pi = atan(done)*4
ccc        call prin2('pi*', pi, 1)

c
c        precomp
c
        k = k1+k2
        call precomp(n, k1, k2, sig2, sig3, a, y, c, c_inv, w, d1)
        call inner_prod(y, y, n, yty)
        a0 = yty/(2*sig1**2) 
ccc        call prin2('d1*', d1, k)

c
c        compute term2
c
        f_const = 0
        do 200 i=1,k
        a1 = w(i)/sig1**2
        a2 = 1/(2*sig1**2) + d1(i)
        f_const = f_const + a1**2/(4*a2) + log(sqrt(2*pi/(2*a2)))
ccc        call prin2('f_const*', f_const, 1)
        dexps_z(i) = a1/(2*a2)
 200    continue

        f_const = f_const - a0 
        f_const = f_const - k1*log(sig2) - k2*log(sig3) - n*log(sig1)
        f_const = f_const - sig1**2/2.0 - sig2**2/2.0 - sig3**2/2.0 

        call mat_mult(c, k, k, dexps_z, k, dexps_x)

        return
        end
c
c
c
c
c
        subroutine corner_ind_3d(nn, i, j, k, coef)
        implicit real *8 (a-h,o-z)

        ii=0
        ij=0
        ik=0
        coef = 1.0d0

        if ((i .eq. 1) .or. (i .eq. nn)) ii=1
        if ((j .eq. 1) .or. (j .eq. nn)) ij=1
        if ((k .eq. 1) .or. (k .eq. nn)) ik=1

        ijk = ii + ij + ik

        if (ijk .eq. 1) coef=1/2.0d0
        if (ijk .eq. 2) coef=1/4.0d0
        if (ijk .eq. 3) coef=1/8.0d0
        

        return
        end
c
c
c
c
c
        subroutine corner_ind_2d(nn, i, j, coef)
        implicit real *8 (a-h,o-z)

        ii=0
        ij=0
        coef = 1

        if ((i .eq. 1) .or. (i .eq. nn)) ii=1
        if ((j .eq. 1) .or. (j .eq. nn)) ij=1

        ijk = ii + ij

        if (ijk .eq. 1) coef=1/2.0d0
        if (ijk .eq. 2) coef=1/4.0d0

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
        subroutine plot_heatmap_old(iw, nn, tol, fs)
        implicit real *8 (a-h,o-z)
        real*8 fs(1)
        data in/11/

        call max_vec(nn*nn, fs, f1, ind)
ccc        call prin2('max*', f1, 1)

        do 320 i=1,nn**2
        fs(i) = fs(i) - f1
        if (fs(i) .lt. f1-tol) fs(i) = -tol
ccc        call prin2('fs*', fs, nnodes**2)
 320    continue

ccc        iw = 11
c
c        if this function gets called twice, nothing gets
c        plotted, so change iw
c
        iw = in
        call pyimage(iw,nn,nn,fs,'title* ')

        in = in + 1

        return
        end
c
c
c
c
c
        subroutine test_pdf()
        implicit real *8 (a-h,o-z)
        real*8 a(1 000 000), x(100 000), y(100 000), b(100 000), 
     1     at(100 000), c_inv(100 000), d2(100 000), c(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     xexps(10 000), xexps2(10 000), ts(10 000), dexps_x(10 000),
     1     dsums(10 000), dints(10 000), dds(10 000)
c
c        parameters
c
        n = 3
        k1 = 2
        k2 = 1
        k = k1+k2
        call corrand_norm(n*k,a,y)
        call corrand_norm(k,x,tmp_mat)
        a(1) = 1.0
        a(2) = 1.0
        a(3) = 1.0
        a(4) = 1.0
        a(5) = 1.0
        a(6) = 0.0
        a(7) = 1.0
        a(8) = 0.0
        a(9) = 0.0

        y(1) = 1.0
        y(2) = 0.1
        y(3) = 0.5

        sig1 = 1.0
        sig2 = 2.0
        sig3 = 1.5
        x(1) = 1.0
        x(2) = 1.5
        x(3) = 0.0
        call q2_log_old(n, k1, k2, a, x, y, sig1, sig2, sig3, 
     1     f1, f_const0, dexps_x)
        call prin2('f1*', f1, 1)

        sig1 = 0.5
        sig2 = 1.0
        sig3 = 0.5
        x(1) = 1
        x(2) = 3.0
        x(3) = 0.0
        call q2_log_old(n, k1, k2, a, x, y, sig1, sig2, sig3, 
     1     f2, f_const0, dexps_x)
        call prin2('f2*', f2, 1)

        call prin2('dd*', f2 - f1, 1)

        return
        end
c
c
c
c
c
        subroutine test2()
        implicit real *8 (a-h,o-z)
        real*8 a(1 000 000), x(100 000), y(100 000), b(100 000), 
     1     at(100 000), c_inv(100 000), d2(100 000), c(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000),
     1     xexps(10 000), xexps2(10 000), ts(10 000), dexps_x(10 000)

c
c        parameters
c
        n = 3
        k1 = 2
        k2 = 1
        k = k1+k2
        sig1 = 1.0
        sig2 = 1.0
        sig3 = 1.0

        call corrand_norm(n*k,a,y)
        call corrand_norm(k,x,tmp_mat)
        a(1) = 1.0
        a(2) = 1.0
        a(3) = 1.0
        a(4) = 1.0
        a(5) = 1.0
        a(6) = 0.0
        a(7) = 1.0
        a(8) = 0.0
        a(9) = 0.0
        a(1) = 1.0
        a(2) = 0.0
        a(3) = 0.0
        a(4) = 0.0
        a(5) = 1.0
        a(6) = 0.0
        a(7) = 0.0
        a(8) = 0.0
        a(9) = 1.0

        y(1) = 1.0
        y(2) = 0.1
        y(3) = 0.5

c
c        find expectations dumb
c

        nn = 37
        t0 = -3.0
        t1 = 5.0
        call equispaced_nodes(nn, t0, t1, ts)

        dsum = 0
        dsum1 = 0
        dsum2 = 0
        do 130 i=1,nn
        do 120 j=1,nn
        do 110 ik=1,nn
        x(1) = ts(i)
        x(2) = ts(j)
        x(3) = ts(ik)
        call q2_log_old(n, k1, k2, a, x, y, sig1, sig2, sig3, 
     1     f, dexps_x)
ccc        call prin2('dexps_x*', dexps_x, k)
        f = exp(f)
        dsum = dsum + f
        dsum1 = dsum1 + x(1)*f
        dsum2 = dsum2 + x(2)*f
        dsum3 = dsum3 + x(3)*f
ccc        call prin2('dsum2*', dsum2, 1)
 110    continue
 120    continue
 130    continue
        
        box_vol = (ts(2) - ts(1))**2

        dconst = dsum * box_vol
        dint1 = dsum1 * box_vol / dconst
        dint2 = dsum2 * box_vol / dconst
        dint3 = dsum3 * box_vol / dconst

ccc        call prin2('dconst*', dconst, 1)
ccc        call prin2('dsum1*', dsum1, 1)
        call prin2('dint1*', dint1, 1)
        call prin2('dint2*', dint2, 1)
        call prin2('dint3*', dint3, 1)
ccc        call prin2('box_vol*', box_vol, 1)

        return
        end
c
c
c
c
c
        subroutine q2_log_old_old(n, k1, k2, a, x, y, sig1, sig2, 
     1     sig3, f)
        implicit real *8 (a-h,o-z)
        real*8 a(1), x(1), y(1), c(10 000), c_inv(10 000),
     1     w(10 000), d1(10 000), z(10 000)

c
c        precomp
c
        k = k1+k2
        call precomp(n, k1, k2, sig2, sig3, a, y, c, c_inv, w, d1)
        call inner_prod(y, y, n, yty)
        a0 = yty/(2*sig1**2) 
        call mat_vec_mult(c_inv, k, k, x, z)

c
c        compute term2
c
        dsum = 0
        do 200 i=1,k
        dsum = dsum + z(i)**2*(1/(2*sig1**2) + d1(i))
        dsum = dsum - z(i)*(w(i)/sig1**2)
 200    continue

        f = -a0 - dsum
        f = f - k1*log(sig2) - k2*log(sig3) - n*log(sig1)
        f = f - sig1**2/2.0 - sig2**2/2.0 - sig3**2/2.0

ccc        call prin2('term1*', term1, 1)
ccc        call prin2('term2*', term2, 1)
ccc        call prin2('f*', f, 1)

        return
        end
c
c
c
c
c
        subroutine precomp(n, k1, k2, sig2, sig3, a, y, c, c_inv, w, d1)
        implicit real *8 (a-h,o-z)
        real*8 a(1), x(100 000), y(1), utd0(100 000),
     1     d0(100 000), ax(100 000), u(100 000), b(100 000), 
     1     a2(100 000), ata(100 000), at(100 000), ut(100 000), 
     1     a_star(100 000), a_start(100 000), d0_sqrt(100 000), 
     1     c_inv(1), z(100 000), u2(100 000), s2(100 000), 
     1     d2(100 000), b2(10 000), a_star_inv(100 000), act(10 000),
     1     d0_sqrt_inv(100 000), a_star_invt(100 000), ac(10 000),
     1     tmp_mat(100 000), d1(1), u2t(100 000), d1_mat(10 000), 
     1     d(100 000), s(10 000), c(1), ct(100 000), w(1)


        k = k1 + k2

        call mat_trans(a, n, k, at)
        call mat_mult(at, k, n, a, k, ata)
ccc        call prin2('ata*', ata, k*k)
ccc        call prin2('y*', y, n)

c
c        eigendecomp of ata
c
        call eigendecomp_ata(k, ata, u, ut, d0)
        call mat_mult(ut, k, k, d0, k, utd0)
        call mat_mult(utd0, k, k, u, k, a2)
ccc        call prin2('u*', u, n*n)
ccc        call prin2('a2*', a2, k*k)
ccc        call prin2('ata*', ata, k*k)

c
c        get a_star and related mats
c

        call get_astars(k, u, ut, d0, a_star, a_start, a_star_inv,
     1     a_star_invt)
ccc        call prin2('a_star*', a_star, k*k)
ccc        call prin2('a_star_inv*', a_star_inv, k*k)

c       
c        construct full diagonal matrices d and d2
c        
        call get_s_s2(k1, k2, sig2, sig3, s, s2)
        call mat_diag(s, k, d)
        call mat_diag(s2, k, d2)
ccc        call prin2('s*', s, k)
ccc        call prin2('s2*', s2, k)

c
c        find eigendecomposition of b = a_star_invt d2 a_star_inv
c
        call mat_mult(a_star_invt, k, k, d2, k, tmp_mat)
        call mat_mult(tmp_mat, k, k, a_star_inv, k, b)
        call eigendecomp_b(k, b, d1, u2, u2t)
ccc        call prin2('b*', b, k*k)
ccc        call prin2('u2*', u2, k*k)
ccc        call prin2('d1*', d1, k)

        call mat_diag(d1, k, d1_mat)
        call mat_mult(u2, k, k, d1_mat, k, utd0)
        call mat_mult(utd0, k, k, u2t, k, b2)
ccc        call prin2('b*', b, k*k)
ccc        call prin2('b2*', b2, k*k)

c
c        construct c, c_inv, ct
c
        call mat_mult(a_star_inv, k, k, u2, k, c)
        call mat_mult(u2t, k, k, a_star, k, c_inv)
ccc        call mat_mult(c, k, k, c_inv, k, tmp_mat)
ccc        call prin2('c*', c, k*k)

c
c        test c^t d2 c = d1
c
ccc        call mat_trans(c, k, k, ct)
ccc        call mat_mult(ct, k, k, d2, k, tmp_mat)
ccc        call mat_mult(tmp_mat, k, k, c, k, b)
ccc        call prin2('b*', b, k*k)
ccc        call prin2('d1*', d1, k)
ccc        call prin2('d2*', d2, k*k)

c
c        precomp w
c
        call mat_mult(a, n, k, c, k, ac)
        call mat_trans(ac, n, k, act)
        call mat_vec_mult(act, k, n, y, w)
ccc        call prin2('ac*', ac, k*k)
ccc        call prin2('w*', w, k)

        return
        end
c
c
c
c
c
        subroutine test1()
        implicit real *8 (a-h,o-z)
        real*8 a(1 000 000), x(100 000), y(100 000), b(100 000), 
     1     at(100 000), c_inv(100 000), d2(100 000), c(100 000), 
     1     tmp_mat(100 000), d1(100 000), d(100 000), s(10 000)

c
c        parameters
c
        n = 3
        k1 = 2
        k2 = 1
        n = 10
        k1 = 6
        k2 = 4
        k = k1 + k2

        n = 3
        k1 = 2
        k2 = 1
        k = k1+k2
        sig1 = 1.0
        sig2 = 2.0
        sig3 = 1.0


        call corrand_norm(n*k,a,y)
        call corrand_norm(k,x,tmp_mat)
        a(1) = 1.0
        a(2) = 1.0
        a(3) = 1.0
        a(4) = 1.0
        a(5) = 1.0
        a(6) = 0.0
        a(7) = 1.0
        a(8) = 0.0
        a(9) = 0.0

        x(1) = 1.0
        x(2) = 2.0
        x(3) = 3.0

        y(1) = 1.0
        y(2) = 1.5
        y(3) = 1.7

c
c        new evaluation
c

        call q2_log_old(n, k1, k2, a, x, y, sig1, sig2, sig3, f)
        call prin2('quad form new*', f, 1)

c
c        compare to slow way
c

        call quad_form_dumb(n, k1, k2, sig1, sig2, 
     1     sig3, a, x, y, tot_dumb)
        call prin2('quad form dumb*', tot_dumb, 1)

        dd = (f - tot_dumb)/f
        call prin2('dd*', dd, 1)

        return
        end
c
c
c
c
c
        subroutine eigendecomp_ata(n, ata, u, ut, d0)
        implicit real *8 (a-h,o-z)
        real*8 ata(1), a_star(1), a_start(1), u(1), ut(1), d0(1)
        real*8 a_star_inv(100 000), a_star_invt(100 000)
        real*8 b(10 000), d0_sqrt(10 000)
c
c        eigendecomposition of ata
c
        ifvect=1
        eps=1.0d-14
        call cc_arr(n*n, ata, b)
        call jaceig(b,n,eps,ifvect,ut,nelems)
  
ccc        call prin2('after jackeig, eigenvalues are*',b,n)
ccc        call prin2('after jackeig, eigenvectors are*',ut,k*k)
ccc        call prinf('after jackeig, nelems=*',nelems,1)

        call mat_trans(ut, n, n, u)
        call mat_diag(b, n, d0)

        
        return
        end
c
c
c
c
c
        subroutine get_astars(n, u, ut, d0, a_star, a_start, a_star_inv,
     1     a_star_invt)
        implicit real *8 (a-h,o-z)
        real*8 d0(1), u(1), a2(100 000), ut(1), a_star(1)
        real*8 a_start(1), d0_sqrt(100 000), a_star_inv(1)
        real*8 a_star_invt(1), tmp_mat(100 000), d0_sqrt_inv(100 000)

c
c        make large diagonal matrices
c
        do 140 i=1,n*n
        d0_sqrt(i) = 0
        d0_sqrt_inv(i) = 0
        if (d0(i) .ne. 0) then
          d0_sqrt_inv(i) = 1/sqrt(d0(i))
          d0_sqrt(i) = sqrt(d0(i))
        endif
 140    continue

c
c        construct a_star and a_start
c

ccc        call prin2('d0_sqrt*', d0_sqrt, n*n)
ccc        call prin2('d0_sqrt_inv*', d0_sqrt_inv, n*n)
ccc        call prin2('d0*', d0, n*n)

        call mat_mult(ut, n, n, d0_sqrt, n, a_start)
        call mat_trans(a_start, n, n, a_star)
ccc        call prin2('a_start*', a_start, n*n)
ccc        call prin2('a_star*', a_star, n*n)
ccc        call mat_mult(a_start, n, n, a_star, n, a2)
ccc        call prin2('a2*', a2, n*n)

c
c        compute a_star_inv and a_star_invt
c
        call mat_mult(d0_sqrt_inv, n, n, u, n, a_star_invt)
        call mat_trans(a_star_invt, n, n, a_star_inv)

ccc        call mat_mult(a_star, n, n, a_star_inv, n, tmp_mat)
ccc        call prin2('tmp_mat*', tmp_mat, n*n)
ccc        call prin2('a_star*', a_star, n*n)
ccc        call prin2('a_start*', a_start, n*n)
ccc        call prin2('a_star_inv*', a_star_inv, n*n)

        return
        end
c
c
c
c
c
        subroutine get_s_s2(k1, k2, sig2, sig3, s, s2)
        implicit real *8 (a-h,o-z)
        real*8 s(1), s2(1)

ccc        call prin2('sig2 = *',sig2,1)
ccc        call prin2('sig3 = *',sig3,1)
c
c        construct the diagonal of the second quadratic form
c

        do 160 i=1,k1
        s2(i) = 1/(2*sig2**2)
        s(i) = 1/sqrt(2*sig2**2)
 160    continue

        do 180 i=1,k2
        s2(k1+i) = 1/(2*sig3**2)
        s(k1+i) = 1/sqrt(2*sig3**2)
 180    continue
ccc        call prin2('s2*', s2, n)


        return
        end
c
c
c
c
c
        subroutine eigendecomp_b(n, b, d1, u2, u2t)
        implicit real *8 (a-h,o-z)
        real*8 b(1), d1(10 000), u2(1), u2t(1), d1_mat(1)

        ifvect=1
        eps=1.0d-14
        call cc_arr(n*n, b, d1)

        call jaceig(d1,n,eps,ifvect,u2,nelems)
ccc        call prin2('d1*', d1, n)
ccc        call prin2('u2*', u2, n*n)

        call mat_trans(u2, n, n, u2t)

        return
        end
c
c
c
c
c
        subroutine quad_form_dumb(n, k1, k2, sig1, sig2, 
     1     sig3, a, x, y, tot)
        implicit real *8 (a-h,o-z)
        real*8 a(1), x(1), y(1), ax(100 000), v(100 000)

c
c        compute first term, dumb
c
        k = k1 + k2
        call mat_vec_mult(a, n, k, x, ax)
ccc        call prin2('ax*', ax, n)

        do 100 i=1,n
        v(i) = ax(i) - y(i)
 100    continue
ccc        call prin2('v*', v, n)

        call inner_prod(v, v, n, vtv)
ccc        call prin2('vtv*', vtv, 1)

        term1 = vtv/(2.0*sig1**2)
ccc        call prin2('term1*', term1, 1)

c
c        compute second term, dumb
c
        term2 = 0
        do 120 i=1,k1
        term2 = term2 + x(i)**2/(2*sig2**2)
 120    continue

        do 130 i=1,k2
        term2 = term2 + x(k1+i)**2/(2*sig3**2)
 130    continue
        
ccc        call prin2('term2*', term2, 1)
        tot = term1 + term2
ccc        call prin2('total dumb*', tot, 1)

c
c        include priors on sigs and prefactors
c
        tot = -tot
        tot = tot - k1*log(sig2) - k2*log(sig3) - n*log(sig1)
        tot = tot - sig1**2/2.0 - sig2**2/2.0 - sig3**2/2.0

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
        subroutine symmetrize_mat(a, n, a2)
        implicit real *8 (a-h,o-z)
        real*8 a(n, n), a2(n, n)

        do 100 i=1,n
        do 120 j=1,n
        a2(i, j) = a(i, j)/2.0 + a(j, i)/2.0
 120    continue
 100    continue

        return
        end
c
c
c
c
c
        subroutine mat_diag(s, n, d_mat)
        implicit real *8 (a-h,o-z)
        real*8 s(1), d_mat(n, n)

        do 1075 i=1,n
        do 1070 j=1,n
        d_mat(i, j) = 0
        if (i .eq. j) then
          d_mat(i, i) = s(i)
        endif
 1070   continue
 1075   continue

        return
        end
c
c
c
c
c
        subroutine equispaced_nodes(n, t0, t1, ts)
        implicit real *8 (a-h,o-z)
        real*8 ts(1)

        h = (t1 - t0)/(n-1)
        do 100 i=1,n
        ts(i) = t0 + (i-1)*h
 100    continue

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

ccc        csv_file = 'means.csv'
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

ccc        csv_file = 'expectations.dat'
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
        subroutine read_dd_exps_stds(file_in, k, dd, dsums, stds)
        implicit real*8 (a-h,o-z)
        real*8 dsums(*), stds(*)
        character(*) file_in

        open (2, file=file_in)

 210    format(f22.16,X,f22.16)

c
c        read error
c
        read(2, 210) dd

c
c        read stds and means
c
        do i=1,k+3
        read(2, 210) dsums(i), stds(i)
        enddo

        return
        end
