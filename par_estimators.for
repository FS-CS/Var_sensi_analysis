      module par_estimators
c   One can use the OpenMP library to parallelize loops containing heavy computations. Loop parallelization is useful to the repeated computation of outputs of the tested model that is usually complex and can therefore be slow.
c      use omp_lib
      double precision :: ap,amin,amax,bp,bmin,bmax,cp,cmin,cmax
c      double precision cp1, cm, cv
      character prefix*50

      contains

c   Subroutines for estimator's computations : as long as ind, ind2 or ind3 > acc (accuracy goal), the code iteratively adds Nsamp points to the sample and re-computes the estimator (because of the estimator's additive nature, there is no need to re-run the model for all points in the sample, only for the new ones). The estimator's value is considered to have converged when its variation in value relative to the previous iteration is found smaller than the accuracy goal (acc) three consecutive times

c   Subroutine for computation of the estimators Si (i=1,2,3 here stored in the variable iz)
      subroutine convEstim(iz,Nsamp,dps,acc,model,stresti,Vires,Sires,
     &meanres,varres)
      implicit none
      integer, intent(in) :: iz,Nsamp,dps
      double precision, intent(in) :: acc
      double precision ind,ind2,ind3,Vip1,Vip2,Vip3,Vip4,Sip1,Sip2,Sip3,
     &Sip4,lz,model
      double precision, dimension(dps,Nsamp,4) :: Asob,Bsob,AB1,AB2,AB3
      double precision, dimension(dps+2,dps,Nsamp,4) :: amaster
      double precision, dimension(dps+2,Nsamp,4), target :: fmaster
      double precision, pointer :: fABi(:,:)
      double precision, dimension(4) :: mean,var
      integer niter,iSi,jSi,kSi,iloop,jloop,kloop
      double precision, intent(out) :: Vires,Sires,meanres,varres
      character striz*1
      character stresti*3

c   Si-estimator computation requires arrays fA,fB and fABi ; depending which i=1,2,3 is chosen in the subroutine call, fABi will change consequently ; we use a pointer for fABi to point towards the corresponding array, which is computed later and then used in the Vi function
      if (iz==1) then
         fABi=>fmaster(3,:,:)
      else if (iz==2) then
         fABi=>fmaster(4,:,:)
      else if (iz==3) then
         fABi=>fmaster(5,:,:)
      else
         write(*,*) char(10),'argument "iz" must be 1, 2 or 3'
         stop
      end if
      write(striz,'(I1)') iz
      if (stresti.eq.'Vi ') then
         write(6,*) 'Estimator = S'//striz
      else if (stresti.eq.'Vti') then
         write(6,*) 'Estimator = St'//striz
      else
         write(6,*) 'argument "stresti" must be "Vi " or "Vti"'
         stop
      end if
      write(6,*) '   mean',char(9),'sigma',char(9),'    S(t)'//striz,
     &char(9),'ind1',char(9),'    ind2',char(9),'ind3'

      niter=0
      ind=1d0
      ind2=1d0
      ind3=1d0
      lz=5E-3
      mean(1)=0d0
      var(1)=0d0
      Vip4=0d0

c   Beginning of the do while loop for convergence check : if at least one of the indicators is greater than the accuracy goal, the subroutine extends the sample in the parameters space by Nsamp points, applies the model to them and (re)computes the mean, variance and other estimators for the extended sample 
      do while (ind.gt.acc .or. ind2.gt.acc .or. ind3.gt.acc)

c   Generation of the arrays Asob,Bsob,AB1,AB2,AB3 for the Sobol sequence : Asob & Bsob are two arrays containing a sample of Nsamp points of the de-normalised phase-space (of dimension dps=3 in this case), while the ABi are a different mix of Asob & Bsob's columns (respectively BsobAsobAsob,AsobBsobAsob,AsobAsobBsob), the whole of them forming a low-discrepancy sequence
c   The different arrays are then stored in a master array : amaster, in order to apply the model to all elements of amaster in a (possibly parallelized) loop
         do iSi=1,4
            do jSi=1,Nsamp
               do kSi=1,dps
c   Random normalised numbers are created for each elements of the arrays Asob & Bsob, then rescaled by revnorm to become numbers between their lower and upper bounds
                  Asob(kSi,jSi,iSi)=dble(rand())
                  Bsob(kSi,jSi,iSi)=dble(rand())
               end do
            end do
         end do

         AB1(1,:,:)=Bsob(1,:,:)
         AB1(2,:,:)=Asob(2,:,:)
         AB1(3,:,:)=Asob(3,:,:)

         AB2(1,:,:)=Asob(1,:,:)
         AB2(2,:,:)=Bsob(2,:,:)
         AB2(3,:,:)=Asob(3,:,:)

         AB3(1,:,:)=Asob(1,:,:)
         AB3(2,:,:)=Asob(2,:,:)
         AB3(3,:,:)=Bsob(3,:,:)

c   amaster is of dimension (dps+2,dps,Nsamp,4) : d+2 corresponding to the dps+2 differents arrays used to generate the Sobol sequence : A,B,AB1,AB2,AB3 for dps=3 ; (dps,Nsamp) is the size of each of these arrays containing Nsamp points in the dps-dimensional parameter phase-space ; 4 is the number of additionnal samples we need to compute the estimators 4 different times and make our triple convergence test
         do iSi=1,4
            do jSi=1,Nsamp
               do kSi=1,dps
                   amaster(1,kSi,jSi,iSi)=Asob(kSi,jSi,iSi)
                   amaster(2,kSi,jSi,iSi)=Bsob(kSi,jSi,iSi)
                   amaster(3,kSi,jSi,iSi)=AB1(kSi,jSi,iSi)
                   amaster(4,kSi,jSi,iSi)=AB2(kSi,jSi,iSi)
                   amaster(5,kSi,jSi,iSi)=AB3(kSi,jSi,iSi)
               end do
            end do
         end do

c   Computation of the model for each value of amaster, stored in fmaster
c   One can parallelize the loop computing the model's ouputs if the latter is slow using the OpenMP library (other loops in the program should be fast)
c      call omp_set_num_threads(2)
c!$omp parallel do
         do iloop=1,4
            do jloop=1,Nsamp
               do kloop=1,dps+2
                  ap=revnorm(amaster(kloop,1,jloop,iloop),amin,amax)
                  bp=revnorm(amaster(kloop,2,jloop,iloop),bmin,bmax)
                  cp=revnorm(amaster(kloop,3,jloop,iloop),cmin,cmax)
                  fmaster(kloop,jloop,iloop)=model()
               end do
            end do
         end do
c!$omp end parallel do

c   Computation of the mean and variance of the model for the 4 different steps of the evaluation, each using a sample Nsamp points-larger than the previous one
         do iSi=1,4
            mean(iSi)=(mean(iSi)+sum(fmaster(:,:,iSi))/((dps+2)*Nsamp))
     &      /2
         end do
         do iSi=1,4
            var(iSi)=(var(iSi)+sum((fmaster(:,:,iSi)-mean(iSi))**2)/(
     &      (dps+2)*Nsamp))/2
         end do

c   Computation of the Vi (with a dedicated funtion) and Si estimators of the model for the 4 different steps of the evaluation, each using a sample that is Nsamp-points larger than the previous one
         if (stresti.eq.'Vi ') then
            Vip1=(Vip4+Vi(fmaster(1,:,1),fmaster(2,:,1),fABi(:,1),
     &      Nsamp))/2
            Sip1=Vip1/var(1)
            Vip2=(Vip1+Vi(fmaster(1,:,2),fmaster(2,:,2),fABi(:,2),
     &      Nsamp))/2
            Sip2=Vip2/var(2)
            Vip3=(Vip2+Vi(fmaster(1,:,3),fmaster(2,:,3),fABi(:,3),
     &      Nsamp))/2
            Sip3=Vip3/var(3)
            Vip4=(Vip3+Vi(fmaster(1,:,4),fmaster(2,:,4),fABi(:,4),
     &      Nsamp))/2
            Sip4=Vip4/var(4)
         else if (stresti.eq.'Vti') then
            Vip1=(Vip4+Vti(fmaster(1,:,1),fmaster(2,:,1),fABi(:,1),
     &      Nsamp))/2
            Sip1=Vip1/var(1)
            Vip2=(Vip1+Vti(fmaster(1,:,2),fmaster(2,:,2),fABi(:,2),
     &      Nsamp))/2
            Sip2=Vip2/var(2)
            Vip3=(Vip2+Vti(fmaster(1,:,3),fmaster(2,:,3),fABi(:,3),
     &      Nsamp))/2
            Sip3=Vip3/var(3)
            Vip4=(Vip3+Vti(fmaster(1,:,4),fmaster(2,:,4),fABi(:,4),
     &      Nsamp))/2
            Sip4=Vip4/var(4)
         end if

c   Computation of the convergence indicators
         ind=abs((Sip2-Sip1)/Sip1)
         ind2=abs((Sip3-Sip2)/Sip2)
         ind3=abs((Sip4-Sip3)/Sip3)
c   Writing the resulting mean, variance and onvergence indicators on screen at the end of the do while loop
         niter=niter+1
         write(*,16) mean(4),sqrt(var(4)),Sip4,ind,ind2,ind3

c   If the estimator is 0 or close to 0, its relative variations will never converge under the accuracy goal ; therefore if the estimator is found to be smaller than lz=0.5% four times in a row, its value is returned as such and can be considered to be 0 or negligible 
         if (abs(Sip1).lt.lz .and. abs(Sip2).lt.lz .and. 
     &   abs(Sip3).lt.lz .and. abs(Sip4).lt.lz) then
            write(6,*) 'Estimator is less than',real(lz),', possibly 0'
            goto 19
         end if
      end do

c   Once convergence is checked, the subroutine exits the loop and returns the last computed values of the estimators (thus the ones coming from the biggest samples) and writes them on screen
 19   Vires=Vip4
      Sires=Sip4
      meanres=mean(4)
      varres=var(4)
      write(*,*) ' niter',char(9),'      mean',char(9),
     &'  sigma',char(9),'      Vi/Vti',char(9),'  Si/Sti'
      write(*,16) niter*1d0,meanres,sqrt(varres),Vires,Sires
      write(*,*) char(10)

 16   format(8(1pe10.2,2x))

      end subroutine convEstim


c This subroutine computes the mean and standard deviation only
      subroutine convMeanVar(Nsamp,dps,acc,model,meanres,varres)
      implicit none
      integer, intent(in) :: Nsamp,dps
      double precision, intent(in) :: acc
      double precision indm,indm2,indm3,indv,indv2,indv3,lz,model
      double precision, dimension(dps,Nsamp,4) :: Asob,Bsob,AB1,AB2,AB3
      double precision, dimension(dps+2,dps,Nsamp,4) :: amaster
      double precision, dimension(dps+2,Nsamp,4), target :: fmaster
      double precision, dimension(4,Nsamp) :: initlist
      double precision, pointer :: fABi(:,:)
      double precision, dimension(4) :: mean,var
      integer niter,iSi,jSi,kSi,iloop,jloop,kloop
      double precision, intent(out) :: meanres,varres

      niter=0
      indm=1d0
      indm2=1d0
      indm3=1d0
      indv=1d0
      indv2=1d0
      indv3=1d0
      lz=5E-3

      write(6,*) 'Computing initial values...',char(10)
      do iSi=1,4
         do jSi=1,Nsamp
         ap=revnorm(rand()*1d0,amin,amax)
         bp=revnorm(rand()*1d0,bmin,bmax)
         cp=revnorm(rand()*1d0,cmin,cmax)
c   Example using rand_normal for a parameter cp1 with a gaussian-based distribution
c         cp1=rand_normal(amaster(kloop,3,jloop,iloop))
         initlist(iSi,jSi)=model()
         end do
      end do
      do iSi=1,4
         mean(iSi)=sum(initlist(iSi,:))/Nsamp
         var(iSi)=sum((initlist(iSi,:)-mean(iSi))**2)/Nsamp
      end do

      write(6,*) '   mean',char(9),'sigma',char(9),'   indm1',
     &char(9),'indm2',char(9),'   indm3',char(9),'indv',char(9)
     &,'    indv2',char(9),'indv3'


      do while (indm.gt.acc .or. indm2.gt.acc .or. indm3.gt.acc .or. 
     &indv.gt.acc .or. indv2.gt.acc .or. indv3.gt.acc)
         do iSi=1,4
            do jSi=1,Nsamp
               do kSi=1,dps
                  Asob(kSi,jSi,iSi)=dble(rand())
                  Bsob(kSi,jSi,iSi)=dble(rand())
               end do
            end do
         end do

         AB1(1,:,:)=Bsob(1,:,:)
         AB1(2,:,:)=Asob(2,:,:)
         AB1(3,:,:)=Asob(3,:,:)

         AB2(1,:,:)=Asob(1,:,:)
         AB2(2,:,:)=Bsob(2,:,:)
         AB2(3,:,:)=Asob(3,:,:)

         AB3(1,:,:)=Asob(1,:,:)
         AB3(2,:,:)=Asob(2,:,:)
         AB3(3,:,:)=Bsob(3,:,:)

         do iSi=1,4
            do jSi=1,Nsamp
               do kSi=1,dps
                   amaster(1,kSi,jSi,iSi)=Asob(kSi,jSi,iSi)
                   amaster(2,kSi,jSi,iSi)=Bsob(kSi,jSi,iSi)
                   amaster(3,kSi,jSi,iSi)=AB1(kSi,jSi,iSi)
                   amaster(4,kSi,jSi,iSi)=AB2(kSi,jSi,iSi)
                   amaster(5,kSi,jSi,iSi)=AB3(kSi,jSi,iSi)
               end do
            end do
         end do

c      call omp_set_num_threads(2)
c!$omp parallel do
         do iloop=1,4
            do jloop=1,Nsamp
               do kloop=1,dps+2
                  ap=revnorm(amaster(kloop,1,jloop,iloop),amin,amax)
                  bp=revnorm(amaster(kloop,2,jloop,iloop),bmin,bmax)
                  cp=revnorm(amaster(kloop,2,jloop,iloop),cmin,cmax)
c   Example using rand_normal for a parameter cp1 with a gaussian-based distribution of mean cm and variance cv
c                  cp=rand_normal(cm,cv,amaster(kloop,3,jloop,iloop))
                  fmaster(kloop,jloop,iloop)=model()
               end do
            end do
         end do
c!$omp end parallel do

         do iSi=1,4
            mean(iSi)=(mean(iSi)+sum(fmaster(:,:,iSi))/((dps+2)*Nsamp))
     &      /2
         end do
         do iSi=1,4
            var(iSi)=(var(iSi)+sum((fmaster(:,:,iSi)-mean(iSi))**2)/(
     &      (dps+2)*Nsamp))/2
         end do

         indm=abs((mean(2)-mean(1))/mean(1))
         indm2=abs((mean(3)-mean(2))/mean(2))
         indm3=abs((mean(4)-mean(3))/mean(3))

         indv=abs((var(2)-var(1))/var(1))
         indv2=abs((var(3)-var(2))/var(2))
         indv3=abs((var(4)-var(3))/var(3))

         niter=niter+1
         write(*,16) mean(4),sqrt(var(4)),indm,indm2,indm3,indv,indv2,
     &   indv3

         if (abs(mean(1)).lt.lz .and. abs(mean(2)).lt.lz .and. 
     &   abs(mean(3)).lt.lz .and. abs(mean(4)).lt.lz) then
            write(6,*) 'Mean is less than',real(lz),', possibly 0'
            goto 19
         end if
         if (abs(var(1)).lt.lz .and. abs(var(2)).lt.lz .and. 
     &   abs(var(3)).lt.lz .and. abs(var(4)).lt.lz) then
            write(6,*) 'Var is less than',real(lz),', possibly 0'
            goto 19
         end if

      end do

 19   meanres=mean(4)
      varres=var(4)

      write(*,*) char(10),'  niter       mean        sigma'
      write(*,16) niter*1d0,meanres,sqrt(varres)

 16   format(8(1pe10.2,2x))

      end subroutine convMeanVar


c**********   Helper functions **********


c   Function revnorm that transforms a normalised variable into a real one 0<xnorm<1 --> xmin<x<xmax
      double precision function revnorm(xnorm,xmin,xmax)
      implicit none
      double precision xnorm,xmin,xmax
      revnorm=(xmax-xmin)*xnorm+xmin
      end function revnorm

c   Function evaluating the 1st-order sensitivity index (or main effect index) Vi (not normalised by the variance of the whole model)
      double precision function Vi(fA,fB,fABi,Nsamp)
      implicit none
      integer Nsamp,iVi
      double precision resVi
      double precision, dimension (Nsamp) :: fA,fB,fABi
      resVi=0d0
      do iVi=1,Nsamp
         resVi=resVi+1d0/Nsamp*fB(iVi)*(fABi(iVi)-fA(iVi))
      end do
      Vi=resVi
      end function Vi

c   Function evaluating the total-order index (or total-effect index) Vti (not normalised by the variance of the whole model)
      double precision function Vti(fA,fB,fABi,Nsamp)
      implicit none
      integer Nsamp,iVti
      double precision resVti
      double precision, dimension (Nsamp) :: fA,fB,fABi
      resVti=0d0
      do iVti=1,Nsamp
         resVti=resVti+1d0/(2d0*Nsamp)*(fA(iVti)-fABi(iVti))**2
      end do
      Vti=resVti
      end function Vti

c   Function evaluating the 2nd-order sensitivity index (or interaction-effect index) Sij (which is normalised by the variance of the whole model)
      double precision function Sij(Si,Sj,Sk,Sti,Stj,Stk)
      implicit none
      double precision Si,Sj,Sk,Sti,Stj,Stk
      Sij=1d0/2d0*(-Si-Sj+Sk+Sti+Stj-Stk)
      end function Sij

c Function that generates random values following a normal distribution of given mean and standard deviation 
c It needs 2 random numbers ; ran1 comes from the Sobol sequence in amaster ; the second one is generated with rand()
      function rand_normal(mean,stdev,ran1) result(c)
      implicit none
      double precision :: mean,stdev,c,r,theta,ran1,pi
      pi=3.14d0
      r=(-2.0d0*log(ran1))**0.5
      theta = 2.0d0*pi*rand()
      c= mean+stdev*r*sin(theta)
      end function rand_normal


      end module par_estimators

